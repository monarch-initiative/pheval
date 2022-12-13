#!/usr/bin/python
import gzip
import logging
import secrets
import shutil
from dataclasses import dataclass
from pathlib import Path

import click

from pheval.utils.phenopacket_utils import (
    IncompatibleGenomeAssemblyError,
    PhenopacketRebuilder,
    PhenopacketUtil,
    ProbandCausativeVariant,
    phenopacket_reader,
)

from .custom_exceptions import InputError, MutuallyExclusiveOptionError
from ..utils.file_utils import all_files, files_with_suffix, is_gzipped

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
fh = logging.FileHandler(r"pheval.log")
logger.addHandler(ch)
logger.addHandler(fh)

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
    """Data obtained from VCF header"""

    sample_id: str
    assembly: str
    chr_status: bool


class VcfPicker:
    """Chooses a VCF file from random for a directory if provided, otherwise selects the single template."""

    def __init__(self, template_vcf: Path or None, vcf_dir: Path or None):
        self.template_vcf = template_vcf
        self.vcf_dir = vcf_dir

    def pick_file_from_dir(self) -> Path:
        return secrets.choice(all_files(self.vcf_dir))

    def pick_file(self) -> Path:
        """Selects a VCF file from random when given a directory, if not, template vcf is assigned."""
        return self.pick_file_from_dir() if self.vcf_dir is not None else self.template_vcf


class VcfParser:
    """Parses a VCF file."""

    def __init__(self, template_vcf_file: Path, file_extension_gz: bool):
        self.file_extension_gz = file_extension_gz
        open_fn = gzip.open if file_extension_gz else open
        self.vcf_file = open_fn(template_vcf_file)

    def close_vcf(self) -> None:
        self.vcf_file.close()

    def parse_assembly(self) -> tuple[str, bool]:
        """Parses the genome assembly and format of vcf_records."""
        vcf_assembly = {}
        chr_status = False
        for line in self.vcf_file:
            if self.file_extension_gz:
                line = line.decode().strip()
            if line.startswith("##contig=<ID"):
                line_split = line.split(",")
                chromosome = line_split[0].split("=")[2]
                if "chr" in chromosome:
                    chr_status = True
                    chromosome = chromosome.replace("chr", "")
                contig_length = line_split[1].split("=")[1]
                vcf_assembly[chromosome] = int(contig_length)
                vcf_assembly = {i: vcf_assembly[i] for i in vcf_assembly if i.isdigit()}
        self.vcf_file.seek(0)
        assembly = [k for k, v in genome_assemblies.items() if v == vcf_assembly][0]
        return assembly, chr_status

    def parse_sample_id(self) -> str:
        """Parses the sample ID of the VCF."""
        self.vcf_file.seek(0)
        for line in self.vcf_file:
            if self.file_extension_gz:
                line = line.decode().strip()
            if line.startswith("#CHROM"):
                return line.split("\t")[9].rstrip()

    def parse_vcf_header(self) -> VcfHeader:
        """Parses the header of the VCF."""
        assembly, chr_status = self.parse_assembly()
        sample_id = self.parse_sample_id()
        self.close_vcf()
        return VcfHeader(sample_id, assembly, chr_status)


class ProbandVariantChecker:
    """Performs checks on the proband variant data."""

    def __init__(
        self,
        proband_causative_variants: list[ProbandCausativeVariant],
        vcf_header: VcfHeader,
    ):
        self.proband_causative_variants = proband_causative_variants
        self.vcf_header = vcf_header
        self.compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]

    def check_variant_assembly(self) -> list[ProbandCausativeVariant]:
        """Checks that proband variant data is compatible with the software and VCF file.
        Returns a list of incorrect variants."""
        incompatible_variants = []
        for variant in self.proband_causative_variants:
            if variant.assembly not in self.compatible_genome_assembly:
                raise IncompatibleGenomeAssemblyError(variant.assembly, variant.phenopacket)
            if variant.assembly != self.vcf_header.assembly:
                incompatible_variants.append(variant)
        return incompatible_variants


class VcfSpiker:
    """Spikes proband variants into template VCF file contents."""

    def __init__(
        self,
        phenopacket: Path,
        vcf_file: Path,
        output_dir: Path,
        proband_causative_variants: list[ProbandCausativeVariant],
        vcf_header: VcfHeader,
    ):
        self.phenopacket = phenopacket
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.proband_causative_variants = proband_causative_variants
        self.vcf_header = vcf_header
        self.proband_vcf_file_name = self.output_dir.joinpath(
            self.phenopacket.name.replace(".json", ".vcf")
        )

    def construct_variant(self, proband_variant_data: ProbandCausativeVariant) -> list[str]:
        """Constructs variant entries."""
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
            proband_variant_data.variant.alt,
            "100",
            "PASS",
            "SPIKED_VARIANT_" + proband_variant_data.genotype.upper(),
            "GT:AD:DP:GQ:PL",
            genotype_codes[proband_variant_data.genotype.lower()] + ":0,2:2:12:180,12,0" + "\n",
        ]

    def construct_vcf_records(self, file_extension_gz: bool) -> list[str]:
        """Inserts spiked variant into correct position within VCF."""
        open_fn = gzip.open if file_extension_gz else open
        vcf_file = open_fn(self.vcf_file)
        vcf_contents = (
            [line.decode() for line in vcf_file.readlines()]
            if file_extension_gz
            else vcf_file.readlines()
        )
        vcf_file.close()
        for variant in self.proband_causative_variants:
            variant = self.construct_variant(variant)
            variant_entry_position = [
                i
                for i, val in enumerate(vcf_contents)
                if val.split("\t")[0] == variant[0] and int(val.split("\t")[1]) < int(variant[1])
            ][-1] + 1

            vcf_contents.insert(variant_entry_position, "\t".join(variant))
        return vcf_contents

    def construct_header(self, file_extension_gz: bool) -> list[str]:
        """Constructs the header of the VCF."""
        updated_vcf = []
        updated_records = self.construct_vcf_records(file_extension_gz)
        for line in updated_records:
            text = line.replace(
                self.vcf_header.sample_id,
                self.proband_causative_variants[0].proband_id,
            )
            updated_vcf.append(text)
        return updated_vcf


class VcfWriter:
    def __init__(
        self,
        template_vcf_name: Path,
        vcf_contents: list[str],
        spiked_vcf_file_name: Path,
        file_extension_gz: bool,
    ):
        self.template_vcf_name = template_vcf_name
        self.vcf_contents = vcf_contents
        self.spiked_vcf_file_name = spiked_vcf_file_name
        self.file_extension_gz = file_extension_gz

    def copy_template_to_new_file(self):
        shutil.copyfile(self.template_vcf_name, self.spiked_vcf_file_name)

    def write_gzip(self):
        encoded_contents = [line.encode() for line in self.vcf_contents]
        with gzip.open(self.spiked_vcf_file_name, "wb") as f:
            for line in encoded_contents:
                f.write(line)
        f.close()

    def write_uncompressed(self):
        with open(self.spiked_vcf_file_name, "w") as file:
            file.writelines(self.vcf_contents)
        file.close()

    def write_vcf_file(self):
        self.copy_template_to_new_file()
        self.write_gzip() if self.file_extension_gz else self.write_uncompressed()


def spike_vcf(phenopacket: Path, output_dir: Path, template_vcf: Path = None, vcf_dir: Path = None):
    """Creates a new VCF file with spiked variants specified in a phenopacket."""
    # this is a separate function to a click command as it will fail if annotated with click annotations
    # and referenced from another click command
    if template_vcf is None and vcf_dir is None:
        raise InputError("Either a template_vcf or vcf_dir must be specified")
    try:
        output_dir.mkdir()
        logger.info(f" Created a directory {output_dir}")
    except FileExistsError:
        pass
    phenopacket_contents = phenopacket_reader(phenopacket)
    phenopacket_causative_variants = PhenopacketUtil(phenopacket_contents).causative_variants(
        phenopacket
    )
    chosen_template_vcf = VcfPicker(template_vcf, vcf_dir).pick_file()
    file_extension_gz = is_gzipped(chosen_template_vcf)
    vcf_header = VcfParser(chosen_template_vcf, file_extension_gz).parse_vcf_header()
    incompatible_variants = ProbandVariantChecker(
        phenopacket_causative_variants, vcf_header
    ).check_variant_assembly()
    if len(incompatible_variants) != 0:
        for incompatible_variant in incompatible_variants:
            logger.error(
                f" Skipping... Proband variant does not match Human Genome Build of VCF: {phenopacket.absolute().name}."
                f" Variant Assembly -> {incompatible_variant.assembly} Expected: {vcf_header.assembly}"
            )
            raise IncompatibleGenomeAssemblyError(
                assembly=incompatible_variant.assembly,
                phenopacket=phenopacket.absolute().name,
            )
    spiked_vcf_file_name = (
        output_dir.joinpath(phenopacket.name.replace(".json", ".vcf.gz"))
        if file_extension_gz
        else output_dir.joinpath(phenopacket.name.replace(".json", ".vcf"))
    )
    vcf_contents = VcfSpiker(
        phenopacket,
        chosen_template_vcf,
        output_dir,
        phenopacket_causative_variants,
        vcf_header,
    ).construct_header(file_extension_gz)

    VcfWriter(
        chosen_template_vcf, vcf_contents, spiked_vcf_file_name, file_extension_gz
    ).write_vcf_file()
    phenopacket_rebuilder = PhenopacketRebuilder(phenopacket_contents)
    phenopacket_rebuilder.add_created_vcf_path(
        spiked_vcf_file_name,
        vcf_header.assembly,
    )
    phenopacket_rebuilder.write_phenopacket(phenopacket)


@click.command()
@click.option(
    "--phenopacket",
    "-p",
    metavar="FILE",
    required=True,
    help="Path to phenopacket file",
    type=Path,
)
@click.option(
    "--template-vcf",
    "-t",
    cls=MutuallyExclusiveOptionError,
    metavar="FILE",
    required=False,
    help="Template VCF file",
    mutually_exclusive=["vcf_dir"],
    type=Path,
)
@click.option(
    "--vcf-dir",
    "-v",
    cls=MutuallyExclusiveOptionError,
    metavar="PATH",
    help="Directory containing template VCF files",
    mutually_exclusive=["template_vcf"],
    type=Path,
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="vcf",
    type=Path,
)
def create_spiked_vcf(
    phenopacket: Path, output_dir: Path, template_vcf: Path = None, vcf_dir: Path = None
):
    """Spikes variants into a template VCF file for a singular phenopacket."""
    spike_vcf(phenopacket, output_dir, template_vcf, vcf_dir)


@click.command()
@click.option(
    "--phenopacket-dir",
    "-p",
    metavar="PATH",
    required=True,
    help="Path to phenopackets directory",
    type=Path,
)
@click.option(
    "--template-vcf",
    "-t",
    cls=MutuallyExclusiveOptionError,
    metavar="PATH",
    required=False,
    help="Template VCF file",
    mutually_exclusive=["vcf_dir"],
    type=Path,
)
@click.option(
    "--vcf-dir",
    "-v",
    cls=MutuallyExclusiveOptionError,
    metavar="PATH",
    help="Directory containing template VCF files",
    mutually_exclusive=["template_vcf"],
    type=Path,
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="vcf",
    type=Path,
)
def create_spiked_vcfs(
    phenopacket_dir: Path,
    output_dir: Path,
    template_vcf: Path = None,
    vcf_dir: Path = None,
):
    """Spikes variants into a template VCF file for a directory of phenopackets."""
    for phenopacket in files_with_suffix(phenopacket_dir, ".json"):
        spike_vcf(phenopacket, output_dir, template_vcf, vcf_dir)

    # or made a lambda one-liner for maximum wtf...
    # [spike_vcf(path, output_dir, template_vcf, vcf_dir) for path in phenopacket_dir.iterdir() if path.suffix ==
    # ".json"]
