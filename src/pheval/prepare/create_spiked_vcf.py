import gzip
import random
import re
import time
import urllib.parse
from copy import copy
from dataclasses import dataclass
from pathlib import Path

from phenopackets import Family, File, Phenopacket

from pheval.prepare.custom_exceptions import InputError
from pheval.utils.file_utils import all_files, files_with_suffix, is_gzipped
from pheval.utils.logger import get_logger
from pheval.utils.phenopacket_utils import (
    IncompatibleGenomeAssemblyError,
    PhenopacketRebuilder,
    PhenopacketUtil,
    ProbandCausativeVariant,
    phenopacket_reader,
    write_phenopacket,
)

logger = get_logger()
genome_assemblies = {
    "GRCh38": {
        "1": 248956422,
        "2": 242193529,
        "3": 198295559,
        "4": 190214555,
        "5": 181538259,
        "6": 170805979,
        "7": 159345973,
        "8": 145138636,
        "9": 138394717,
        "10": 133797422,
        "11": 135086622,
        "12": 133275309,
        "13": 114364328,
        "14": 107043718,
        "15": 101991189,
        "16": 90338345,
        "17": 83257441,
        "18": 80373285,
        "19": 58617616,
        "20": 64444167,
        "21": 46709983,
        "22": 50818468,
    },
    "GRCh37": {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
    },
}


@dataclass
class VcfHeader:
    """Data obtained from VCF header.

    Args:
        sample_id (str): The sample identifier from the VCF header.
        assembly (str): The assembly information obtained from the VCF header.
        chr_status (bool): A boolean indicating whether the VCF denotes chromosomes as chr or not.
    """

    sample_id: str
    assembly: str
    chr_status: bool


def read_vcf(vcf_file: Path) -> list[str]:
    """
    Read the contents of a VCF file into memory, handling both uncompressed and gzipped files.

    Args:
        vcf_file (Path): The path to the VCF file to be read.

    Returns:
        List[str]: A list containing the lines of the VCF file.
    """
    open_fn = gzip.open if is_gzipped(vcf_file) else open
    vcf = open_fn(vcf_file)
    vcf_contents = [line.decode() for line in vcf.readlines()] if is_gzipped(vcf_file) else vcf.readlines()
    vcf.close()
    return vcf_contents


class VcfHeaderParser:
    """Class for parsing the header of a VCF file."""

    def __init__(self, vcf_contents: list[str]):
        """
        Initialise the VcfHeaderParser.

        Args:
            vcf_contents (list[str]): The contents of the VCF file as a list of strings.
        """
        self.vcf_contents = vcf_contents

    def parse_assembly(self) -> tuple[str, bool]:
        """
        Parse the genome assembly and format of vcf_records.

        Returns:
            Tuple[str, bool]: A tuple containing the assembly and chromosome status (True/False).
        """
        vcf_assembly = {}
        chr_status = False
        for line in self.vcf_contents:
            if line.startswith("##contig=<ID"):
                tokens = line.split(",")
                chromosome = re.sub(r"^.*?ID=", "", next(token for token in tokens if "ID=" in token))
                if "chr" in chromosome:
                    chr_status = True
                    chromosome = chromosome.replace("chr", "")
                contig_length = re.sub("[^0-9]+", "", next(token for token in tokens if "length=" in token))
                vcf_assembly[chromosome] = int(contig_length)
                vcf_assembly = {i: vcf_assembly[i] for i in vcf_assembly if i.isdigit()}
        assembly = next(k for k, v in genome_assemblies.items() if v == vcf_assembly)
        return assembly, chr_status

    def parse_sample_id(self) -> str:
        """
        Parse the sample ID of the VCF.

        Returns:
            str: The sample ID extracted from the VCF header.
        """
        for line in self.vcf_contents:
            if line.startswith("#CHROM"):
                return line.split("\t")[9].rstrip()

    def parse_vcf_header(self) -> VcfHeader:
        """
        Parse the header of the VCF.

        Returns:
            VcfHeader: An instance of VcfHeader containing sample ID, assembly, and chromosome status.
        """
        assembly, chr_status = self.parse_assembly()
        sample_id = self.parse_sample_id()
        return VcfHeader(sample_id, assembly, chr_status)


@dataclass
class VcfFile:
    """
    Represents a VCF file with its name, contents, and header information.

    Attributes:
        vcf_file_name (str): The name of the VCF file.
        vcf_contents (List[str]): The contents of the VCF file.
        vcf_header (VcfHeader): The parsed header information of the VCF file.
    """

    vcf_file_name: str = None
    vcf_contents: list[str] = None
    vcf_header: VcfHeader = None

    @staticmethod
    def populate_fields(template_vcf: Path):
        """
        Populate the fields of the VcfFile instance using the contents of a template VCF file.

        Args:
            template_vcf (Path): The path to the template VCF file.

        Returns:
            VcfFile: An instance of VcfFile with populated fields.

        """
        contents = read_vcf(template_vcf)
        return VcfFile(template_vcf.name, contents, VcfHeaderParser(contents).parse_vcf_header())


def select_vcf_template(
    phenopacket_path: Path,
    proband_causative_variants: list[ProbandCausativeVariant],
    hg19_vcf_info: VcfFile,
    hg38_vcf_info: VcfFile,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> VcfFile:
    """
    Select the appropriate VCF template based on the assembly information of the proband causative variants.

    Args:
        phenopacket_path (Path): The path to the Phenopacket file.
        proband_causative_variants (List[ProbandCausativeVariant]): A list of causative variants from the proband.
        hg19_vcf_info (VcfFile): VCF file info for hg19 template vcf.
        hg38_vcf_info (VcfFile): CF file info for hg38 template vcf.
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files.
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files.

    Returns:
        VcfFile: The selected VCF template file based on the assembly information of the proband causative variants.

    """
    if proband_causative_variants[0].assembly in ["hg19", "GRCh37"]:
        if hg19_vcf_info:
            return hg19_vcf_info
        elif hg19_vcf_dir:
            return VcfFile.populate_fields(random.choice(all_files(hg19_vcf_dir)))
        else:
            raise InputError("Must specify hg19 template VCF!")
    elif proband_causative_variants[0].assembly in ["hg38", "GRCh38"]:
        if hg38_vcf_info:
            return hg38_vcf_info
        elif hg38_vcf_dir:
            return VcfFile.populate_fields(random.choice(all_files(hg38_vcf_dir)))
        else:
            raise InputError("Must specify hg38 template VCF!")
    else:
        raise IncompatibleGenomeAssemblyError(proband_causative_variants[0].assembly, phenopacket_path)


def check_variant_assembly(
    proband_causative_variants: list[ProbandCausativeVariant],
    vcf_header: VcfHeader,
    phenopacket_path: Path,
) -> None:
    """
    Check the assembly of the variant assembly against the VCF.

    Args:
        proband_causative_variants (List[ProbandCausativeVariant]): A list of causative variants from the proband.
        vcf_header (VcfHeader): An instance of VcfHeader representing the VCF file's header.
        phenopacket_path (Path): The path to the Phenopacket file.

    Raises:
        ValueError: If there are too many or incompatible genome assemblies found.
        IncompatibleGenomeAssemblyError: If the assembly in the Phenopacket does not match the VCF assembly.
    """
    compatible_genome_assembly = {"GRCh37", "hg19", "GRCh38", "hg38"}
    phenopacket_assembly = list({variant.assembly for variant in proband_causative_variants})
    if len(phenopacket_assembly) > 1:
        raise ValueError("Too many genome assemblies!")
    if phenopacket_assembly[0] not in compatible_genome_assembly:
        raise IncompatibleGenomeAssemblyError(phenopacket_assembly, phenopacket_path)
    if (phenopacket_assembly[0] in {"hg19", "GRCh37"} and vcf_header.assembly not in {"hg19", "GRCh37"}) or (
        phenopacket_assembly[0] in {"hg38", "GRCh38"} and vcf_header.assembly not in {"hg38", "GRCh38"}
    ):
        raise IncompatibleGenomeAssemblyError(assembly=phenopacket_assembly, phenopacket=phenopacket_path)


class VcfSpiker:
    """Class for spiking proband variants into template VCF file contents."""

    def __init__(
        self,
        vcf_contents: list[str],
        proband_causative_variants: list[ProbandCausativeVariant],
        vcf_header: VcfHeader,
    ):
        """
        Initialise the VcfSpiker.

        Args:
            vcf_contents (List[str]): Contents of the template VCF file.
            proband_causative_variants (List[ProbandCausativeVariant]): List of proband causative variants.
            vcf_header (VcfHeader): The VCF header information.
        """
        self.vcf_contents = vcf_contents
        self.proband_causative_variants = proband_causative_variants
        self.vcf_header = vcf_header

    def construct_variant_entry(self, proband_variant_data: ProbandCausativeVariant) -> list[str]:
        """
        Construct variant entries.

        Args:
            proband_variant_data (ProbandCausativeVariant): Data for the proband variant.

        Returns:
            List[str]: Constructed variant entry as a list of strings.
        """
        genotype_codes = {
            "hemizygous": "0/1",
            "homozygous": "1/1",
            "heterozygous": "0/1",
            "compound heterozygous": "0/1",
        }
        if self.vcf_header.chr_status is True and "chr" not in proband_variant_data.variant.chrom:
            proband_variant_data.variant.chrom = "chr" + proband_variant_data.variant.chrom
        return [
            proband_variant_data.variant.chrom,
            str(proband_variant_data.variant.pos),
            ".",
            proband_variant_data.variant.ref,
            (
                f"<{proband_variant_data.variant.alt}>"
                if proband_variant_data.variant.ref == "N"
                else proband_variant_data.variant.alt
            ),
            "100",
            "PASS",
            proband_variant_data.info if proband_variant_data.info else ".",
            "GT",
            genotype_codes[proband_variant_data.genotype.lower()] + "\n",
        ]

    def construct_vcf_records(self, template_vcf_name: str) -> list[str]:
        """
        Construct updated VCF records by inserting spiked variants into the correct positions within the VCF.

        Args:
            template_vcf_name (str): Name of the template VCF file.

        Returns:
            List[str]: Updated VCF records containing the spiked variants.
        """
        updated_vcf_records = copy(self.vcf_contents)
        for variant in self.proband_causative_variants:
            variant_entry = self.construct_variant_entry(variant)
            matching_indices = [
                i
                for i, val in enumerate(updated_vcf_records)
                if val.split("\t")[0] == variant_entry[0] and int(val.split("\t")[1]) < int(variant_entry[1])
            ]
            if matching_indices:
                logger.info(
                    f"Successfully spiked variant {variant.variant.chrom}-{variant.variant.pos}-"
                    f"{variant.variant.ref}-{variant.variant.alt} in {template_vcf_name}"
                )
                variant_entry_position = matching_indices[-1] + 1
            else:
                logger.warning(
                    f"Could not find entry position for {variant.variant.chrom}-{variant.variant.pos}-"
                    f"{variant.variant.ref}-{variant.variant.alt} in {template_vcf_name}, "
                    "inserting at end of VCF contents."
                )
                variant_entry_position = len(updated_vcf_records)
            updated_vcf_records.insert(variant_entry_position, "\t".join(variant_entry))
        return updated_vcf_records

    def construct_header(self, updated_vcf_records: list[str]) -> list[str]:
        """
        Construct the header of the VCF.

        Args:
            updated_vcf_records (List[str]): Updated VCF records.

        Returns:
            List[str]: Constructed header as a list of strings.
        """
        updated_vcf_file = []
        for line in updated_vcf_records:
            if line.startswith("#"):
                text = line.replace(
                    self.vcf_header.sample_id,
                    self.proband_causative_variants[0].proband_id,
                )
            else:
                text = line
            updated_vcf_file.append(text)
        return updated_vcf_file

    def construct_vcf(self, template_vcf_name: str) -> list[str]:
        """
        Construct the entire spiked VCF file by incorporating the spiked variants into the VCF.

        Args:
            template_vcf_name (str): Name of the template VCF file.

        Returns:
            List[str]: The complete spiked VCF file content as a list of strings.
        """
        return self.construct_header(self.construct_vcf_records(template_vcf_name))


class VcfWriter:
    """Class for writing VCF file."""

    def __init__(
        self,
        vcf_contents: list[str],
        spiked_vcf_file_path: Path,
    ):
        """
        Initialise the VcfWriter class.

        Args:
            vcf_contents (List[str]): Contents of the VCF file to be written.
            spiked_vcf_file_path (Path): Path to the spiked VCF file to be created.
        """
        self.vcf_contents = vcf_contents
        self.spiked_vcf_file_path = spiked_vcf_file_path

    def write_gzip(self) -> None:
        """
        Write the VCF contents to a gzipped VCF file.
        """
        encoded_contents = [line.encode() for line in self.vcf_contents]
        with gzip.open(self.spiked_vcf_file_path, "wb") as f:
            for line in encoded_contents:
                f.write(line)
        f.close()

    def write_uncompressed(self) -> None:
        """
        Write the VCF contents to an uncompressed VCF file.
        """
        with open(self.spiked_vcf_file_path, "w") as file:
            file.writelines(self.vcf_contents)
        file.close()

    def write_vcf_file(self) -> None:
        """
        Write the VCF file based on compression type.

        Determines the file writing method based on the compression type of the spiked VCF file path.
        Writes the VCF contents to the corresponding file format (gzip or uncompressed).
        """
        self.write_gzip() if is_gzipped(self.spiked_vcf_file_path) else self.write_uncompressed()


def spike_vcf_contents(
    phenopacket: Phenopacket | Family,
    phenopacket_path: Path,
    hg19_vcf_info: VcfFile,
    hg38_vcf_info: VcfFile,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> tuple[str, list[str]]:
    """
    Spike VCF records with variants obtained from a Phenopacket or Family.

    Args:
        phenopacket (Union[Phenopacket, Family]): Phenopacket or Family containing causative variants.
        phenopacket_path (Path): Path to the Phenopacket file.
        hg19_vcf_info (VcfFile): VCF file info for hg19 template vcf.
        hg38_vcf_info (VcfFile): VCF file info for hg38 template vcf.
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files.
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files.

    Returns:
        A tuple containing:
            assembly (str): The genome assembly information extracted from VCF header.
            modified_vcf_contents (List[str]): Modified VCF records with spiked variants.
    """
    phenopacket_causative_variants = PhenopacketUtil(phenopacket).causative_variants()
    chosen_template_vcf = select_vcf_template(
        phenopacket_path,
        phenopacket_causative_variants,
        hg19_vcf_info,
        hg38_vcf_info,
        hg19_vcf_dir,
        hg38_vcf_dir,
    )
    check_variant_assembly(phenopacket_causative_variants, chosen_template_vcf.vcf_header, phenopacket_path)
    return (
        chosen_template_vcf.vcf_header.assembly,
        VcfSpiker(
            chosen_template_vcf.vcf_contents,
            phenopacket_causative_variants,
            chosen_template_vcf.vcf_header,
        ).construct_vcf(chosen_template_vcf.vcf_file_name),
    )


def generate_spiked_vcf_file(
    output_dir: Path,
    phenopacket: Phenopacket | Family,
    phenopacket_path: Path,
    hg19_vcf_info: VcfFile,
    hg38_vcf_info: VcfFile,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> File:
    """
    Write spiked VCF contents to a new file.

    Args:
        output_dir (Path): Path to the directory to store the generated file.
        phenopacket (Union[Phenopacket, Family]): Phenopacket or Family containing causative variants.
        phenopacket_path (Path): Path to the Phenopacket file.
        hg19_vcf_info (VcfFile): VCF file info for hg19 template vcf.
        hg38_vcf_info (VcfFile): VCF file info for hg38 template vcf.
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files.
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files.
    Returns:
        File: The generated File object representing the newly created spiked VCF file.
    """
    vcf_assembly, spiked_vcf = spike_vcf_contents(
        phenopacket, phenopacket_path, hg19_vcf_info, hg38_vcf_info, hg19_vcf_dir, hg38_vcf_dir
    )
    spiked_vcf_path = output_dir.joinpath(phenopacket_path.name.replace(".json", ".vcf.gz"))
    VcfWriter(spiked_vcf, spiked_vcf_path).write_vcf_file()
    return File(
        uri=urllib.parse.unquote(spiked_vcf_path.as_uri()),
        file_attributes={"fileFormat": "vcf", "genomeAssembly": vcf_assembly},
    )


def spike_and_update_phenopacket(
    hg19_vcf_info: VcfFile,
    hg38_vcf_info: VcfFile,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
    output_dir: Path,
    phenopacket_path: Path,
) -> None:
    """
    Spike the VCF files with genetic variants relevant to the provided Phenopacket, update the Phenopacket
    accordingly, and write the updated Phenopacket to the specified output directory.

    Args:
        hg19_vcf_info (VcfFile): VCF file info for hg19 template vcf.
        hg38_vcf_info (VcfFile): VCF file info for hg38 template vcf.
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files.
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files.
        output_dir (Path): Directory where the updated Phenopacket will be saved.
        phenopacket_path (Path): Path to the original Phenopacket file.

    Returns:
        None
    """
    phenopacket = phenopacket_reader(phenopacket_path)
    spiked_vcf_file_message = generate_spiked_vcf_file(
        output_dir,
        phenopacket,
        phenopacket_path,
        hg19_vcf_info,
        hg38_vcf_info,
        hg19_vcf_dir,
        hg38_vcf_dir,
    )
    updated_phenopacket = PhenopacketRebuilder(phenopacket).add_spiked_vcf_path(spiked_vcf_file_message)
    write_phenopacket(updated_phenopacket, phenopacket_path)


def create_spiked_vcf(
    output_dir: Path,
    phenopacket_path: Path,
    hg19_template_vcf: Path,
    hg38_template_vcf: Path,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> None:
    """
    Create a spiked VCF for a Phenopacket.

    Args:
        output_dir (Path): The directory to store the generated spiked VCF file.
        phenopacket_path (Path): Path to the Phenopacket file.
        hg19_template_vcf (Path): Path to the hg19 template VCF file (optional).
        hg38_template_vcf (Path): Path to the hg38 template VCF file (optional).
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files (optional).
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files (optional).

    Raises:
        InputError: If both hg19_template_vcf and hg38_template_vcf are None.
    """
    if hg19_template_vcf is None and hg38_template_vcf is None:
        raise InputError("Either a hg19 template vcf or hg38 template vcf must be specified")
    hg19_vcf_info = VcfFile.populate_fields(hg19_template_vcf) if hg19_template_vcf else None
    hg38_vcf_info = VcfFile.populate_fields(hg38_template_vcf) if hg38_template_vcf else None
    spike_and_update_phenopacket(hg19_vcf_info, hg38_vcf_info, hg19_vcf_dir, hg38_vcf_dir, output_dir, phenopacket_path)


def create_spiked_vcfs(
    output_dir: Path,
    phenopacket_dir: Path,
    hg19_template_vcf: Path,
    hg38_template_vcf: Path,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> None:
    """
    Create a spiked VCF for a directory of Phenopackets.

    Args:
        output_dir (Path): The directory to store the generated spiked VCF file.
        phenopacket_dir (Path): Path to the Phenopacket directory.
        hg19_template_vcf (Path): Path to the template hg19 VCF file (optional).
        hg38_template_vcf (Path): Path to the template hg19 VCF file (optional).
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files (optional).
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files (optional).

    Raises:
        InputError: If both hg19_template_vcf and hg38_template_vcf are None.
    """
    if hg19_template_vcf is None and hg38_template_vcf is None and hg19_vcf_dir is None and hg38_vcf_dir is None:
        raise InputError("Need to specify a VCF!")
    hg19_vcf_info = VcfFile.populate_fields(hg19_template_vcf) if hg19_template_vcf else None
    hg38_vcf_info = VcfFile.populate_fields(hg38_template_vcf) if hg38_template_vcf else None
    for phenopacket_path in files_with_suffix(phenopacket_dir, ".json"):
        logger.info(f"Creating spiked VCF for: {phenopacket_path.name}")
        spike_and_update_phenopacket(
            hg19_vcf_info, hg38_vcf_info, hg19_vcf_dir, hg38_vcf_dir, output_dir, phenopacket_path
        )


def spike_vcfs(
    output_dir: Path,
    phenopacket_path: Path,
    phenopacket_dir: Path,
    hg19_template_vcf: Path,
    hg38_template_vcf: Path,
    hg19_vcf_dir: Path,
    hg38_vcf_dir: Path,
) -> None:
    """
    Create spiked VCF from either a Phenopacket or a Phenopacket directory.

    Args:
        output_dir (Path): The directory to store the generated spiked VCF file(s).
        phenopacket_path (Path): Path to a single Phenopacket file (optional).
        phenopacket_dir (Path): Path to a directory containing Phenopacket files (optional).
        hg19_template_vcf (Path): Path to the hg19 template VCF file (optional).
        hg38_template_vcf (Path): Path to the hg38 template VCF file (optional).
        hg19_vcf_dir (Path): The directory containing the hg19 VCF files (optional).
        hg38_vcf_dir (Path): The directory containing the hg38 VCF files (optional).
    """
    start_time = time.perf_counter()
    logger.info("Creating spiked VCFs.")
    output_dir.mkdir(exist_ok=True)
    logger.info(f" Created output directory: {output_dir}")
    if phenopacket_path is not None:
        logger.info(f"Spiking variants from {phenopacket_path}.")
        create_spiked_vcf(
            output_dir,
            phenopacket_path,
            hg19_template_vcf,
            hg38_template_vcf,
            hg19_vcf_dir,
            hg38_vcf_dir,
        )
    elif phenopacket_dir is not None:
        logger.info(f"Spiking variants from {len(all_files(phenopacket_dir))} phenopackets in {phenopacket_dir}.")
        create_spiked_vcfs(
            output_dir,
            phenopacket_dir,
            hg19_template_vcf,
            hg38_template_vcf,
            hg19_vcf_dir,
            hg38_vcf_dir,
        )
    logger.info(f"Finished spiking! Total time: {time.perf_counter() - start_time:.2f} seconds.")
