#!/usr/bin/python
import logging
import os
import secrets
import shutil
from dataclasses import dataclass
from pathlib import Path

import click

from pheval.utils.phenopacket_utils import (
    CausativeVariant,
    IncompatibleGenomeAssemblyError,
    PhenopacketRebuilder,
    PhenopacketUtil,
    phenopacket_reader,
)

from .custom_exceptions import InputError, MutuallyExclusiveOptionError
from ..utils.file_utils import all_files, files_with_suffix

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

    def __init__(self, template_vcf: Path, vcf_dir: Path):
        self.template_vcf = template_vcf
        self.vcf_dir = vcf_dir

    def pick_file(self) -> Path:
        """Selects a VCF file from random when given a directory, if not, template vcf is assigned."""
        vcf_file = self.template_vcf
        if self.vcf_dir is not None:
            vcf_files = all_files(self.vcf_dir)
            vcf_file = secrets.choice(vcf_files)
        return vcf_file


class VcfParser:
    """Parses a VCF file."""

    def __init__(self, template_vcf_file: Path):
        self.vcf_file = open(template_vcf_file)

    def parse_assembly(self) -> tuple[str, bool]:
        """Parses the genome assembly and format of vcf_records."""
        assembly_dict = {}
        chr_status = False
        for line in self.vcf_file:
            if line.startswith("##contig=<ID"):
                line_split = line.split(",")
                chromosome = line_split[0].split("=")[2]
                if "chr" in chromosome:
                    chr_status = True
                    chromosome = chromosome.replace("chr", "")
                length = line_split[1].split("=")[1]
                assembly_dict[chromosome] = int(length)
                assembly_dict = {i: assembly_dict[i] for i in assembly_dict if i.isdigit()}
        self.vcf_file.seek(0)
        assembly = [k for k, v in genome_assemblies.items() if v == assembly_dict][0]
        return assembly, chr_status

    def parse_sample_id(self) -> str:
        """Parses the sample ID of the VCF."""
        for line in self.vcf_file:
            if line.startswith("#CHROM"):
                return line.split("\t")[9].rstrip()

    def parse_header(self) -> VcfHeader:
        """Parses the header of the VCF."""
        assembly, chr_status = self.parse_assembly()
        sample_id = self.parse_sample_id()
        self.vcf_file.close()
        return VcfHeader(sample_id, assembly, chr_status)


class ProbandVariantChecker:
    """Performs checks on the proband variant data."""

    def __init__(
        self,
        list_of_variant_proband_data: list[CausativeVariant],
        vcf_header: VcfHeader,
    ):
        self.list_of_variant_proband_data = list_of_variant_proband_data
        self.vcf_header = vcf_header
        self.compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]

    def check_variant_assembly(self) -> list[CausativeVariant]:
        """Checks that proband variant data is compatible with the software and VCF file.
        Returns a list of incorrect variants."""
        incorrect_variants = []
        for variant in self.list_of_variant_proband_data:
            if variant.assembly not in self.compatible_genome_assembly:
                raise IncompatibleGenomeAssemblyError(variant.assembly, variant.phenopacket)
            if variant.assembly != self.vcf_header.assembly:
                incorrect_variants.append(variant)
        return incorrect_variants


class VcfSpiker:
    """Spikes proband variants into template VCF file contents."""

    def __init__(
        self,
        phenopacket: Path,
        vcf_file: Path,
        output_dir: Path,
        list_of_variant_proband_data: list[CausativeVariant],
        vcf_header: VcfHeader,
    ):
        self.phenopacket = phenopacket
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.list_of_proband_variant_data = list_of_variant_proband_data
        self.vcf_header = vcf_header
        self.proband_vcf_file_name = os.path.join(
            self.output_dir, self.phenopacket.name.replace(".json", ".vcf")
        )

    def construct_variant(self, variant: CausativeVariant) -> list[str]:
        """Constructs variant entries."""
        genotype_codes = {
            "hemizygous": "0/1",
            "homozygous": "1/1",
            "heterozygous": "0/1",
            "compound heterozygous": "0/1",
        }
        if self.vcf_header.chr_status is True and "chr" not in variant.chrom:
            variant.chrom = "chr" + variant.chrom
        variant_data = [
            variant.chrom,
            str(variant.pos),
            ".",
            variant.ref,
            variant.alt,
            "100",
            "PASS",
            "SPIKED_VARIANT_" + variant.genotype.upper(),
            "GT:AD:DP:GQ:PL",
            genotype_codes[variant.genotype.lower()] + ":0,2:2:12:180,12,0" + "\n",
        ]
        return variant_data

    def construct_vcf_records(self) -> list[str]:
        """Inserts spiked variant into correct position within VCF."""
        vcf_file = open(self.vcf_file)
        vcf_contents = vcf_file.readlines()
        vcf_file.close()
        for variant in self.list_of_proband_variant_data:
            variant = self.construct_variant(variant)
            locs = [
                i
                for i, val in enumerate(vcf_contents)
                if val.split("\t")[0] == variant[0] and int(val.split("\t")[1]) < int(variant[1])
            ][-1] + 1
            vcf_contents.insert(locs, "\t".join(variant))
        return vcf_contents

    def construct_header(self) -> list[str]:
        """Constructs the header of the VCF."""
        updated_vcf = []
        updated_records = self.construct_vcf_records()
        for line in updated_records:
            text = line.replace(
                self.vcf_header.sample_id,
                self.list_of_proband_variant_data[0].proband_id,
            )
            updated_vcf.append(text)
        return updated_vcf


class VcfWriter:
    def __init__(self, template_vcf_name: Path, vcf_contents: list[str], file_name: Path):
        self.template_vcf_name = template_vcf_name
        self.vcf_contents = vcf_contents
        self.file_name = file_name

    def copy_template_to_new_file(self):
        shutil.copyfile(self.template_vcf_name, self.file_name)

    def write_vcf(self):
        self.copy_template_to_new_file()
        with open(self.file_name, "w") as file:
            file.writelines(self.vcf_contents)
        file.close()


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
    phenopacket_proband_variant_data = PhenopacketUtil(phenopacket_contents).causative_variants(
        phenopacket
    )
    # TODO: check that chosen template is what is expected (VcfPicker)
    chosen_template_vcf = VcfPicker(template_vcf, vcf_dir).pick_file()
    vcf_header = VcfParser(chosen_template_vcf).parse_header()
    incompatible_variants = ProbandVariantChecker(
        phenopacket_proband_variant_data, vcf_header
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
    vcf_contents = VcfSpiker(
        phenopacket,
        chosen_template_vcf,
        output_dir,
        phenopacket_proband_variant_data,
        vcf_header,
    ).construct_header()

    VcfWriter(
        chosen_template_vcf,
        vcf_contents,
        Path(os.path.join(output_dir, phenopacket.name.replace(".json", ".vcf"))),
    ).write_vcf()
    phenopacket_rebuilder = PhenopacketRebuilder(phenopacket_contents)
    phenopacket_rebuilder.add_created_vcf_path(
        Path(os.path.join(output_dir, phenopacket.name.replace(".json", ".vcf"))),
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
