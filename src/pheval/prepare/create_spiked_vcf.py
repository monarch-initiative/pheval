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
    write_phenopacket,
)
from phenopackets import Phenopacket, Family
from .custom_exceptions import InputError, MutuallyExclusiveOptionError
from ..utils.file_utils import all_files, files_with_suffix, is_gzipped

info_log = logging.getLogger("info")

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
        """Selects a file from a directory at random."""
        return secrets.choice(all_files(self.vcf_dir))

    def pick_file(self) -> Path:
        """Selects a VCF file from random when given a directory, if not, template vcf is assigned."""
        return self.pick_file_from_dir() if self.vcf_dir is not None else self.template_vcf


def read_vcf(vcf_file: Path) -> list[str]:
    """Reads the contents of a VCF file into memory - handles both uncompressed and gzipped."""
    open_fn = gzip.open if is_gzipped(vcf_file) else open
    vcf = open_fn(vcf_file)
    vcf_contents = (
        [line.decode() for line in vcf.readlines()] if is_gzipped(vcf_file) else vcf.readlines()
    )
    vcf.close()
    return vcf_contents


class VcfHeaderParser:
    """Parses the header of a VCF file."""

    def __init__(self, vcf_contents: list[str]):
        self.vcf_contents = vcf_contents

    def parse_assembly(self) -> tuple[str, bool]:
        """Parses the genome assembly and format of vcf_records."""
        vcf_assembly = {}
        chr_status = False
        for line in self.vcf_contents:
            if line.startswith("##contig=<ID"):
                line_split = line.split(",")
                chromosome = line_split[0].split("=")[2]
                if "chr" in chromosome:
                    chr_status = True
                    chromosome = chromosome.replace("chr", "")
                contig_length = line_split[1].split("=")[1]
                vcf_assembly[chromosome] = int(contig_length)
                vcf_assembly = {i: vcf_assembly[i] for i in vcf_assembly if i.isdigit()}
        assembly = [k for k, v in genome_assemblies.items() if v == vcf_assembly][0]
        return assembly, chr_status

    def parse_sample_id(self) -> str:
        """Parses the sample ID of the VCF."""
        for line in self.vcf_contents:
            if line.startswith("#CHROM"):
                return line.split("\t")[9].rstrip()

    def parse_vcf_header(self) -> VcfHeader:
        """Parses the header of the VCF."""
        assembly, chr_status = self.parse_assembly()
        sample_id = self.parse_sample_id()
        return VcfHeader(sample_id, assembly, chr_status)


def check_variant_assembly(proband_causative_variants: list[ProbandCausativeVariant], vcf_header: VcfHeader,
                           phenopacket: Path):
    """Checks the assembly of the variant assembly against the VCF."""
    compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
    for variant in proband_causative_variants:
        if variant.assembly not in compatible_genome_assembly:
            raise IncompatibleGenomeAssemblyError(variant.assembly, phenopacket)
        if variant.assembly != vcf_header.assembly:
            raise IncompatibleGenomeAssemblyError(assembly=variant.assembly,
                                                  phenopacket=phenopacket)


class VcfSpiker:
    """Spikes proband variants into template VCF file contents."""

    def __init__(
            self,
            vcf_contents: list[str],
            proband_causative_variants: list[ProbandCausativeVariant],
            vcf_header: VcfHeader,
    ):
        self.vcf_contents = vcf_contents
        self.proband_causative_variants = proband_causative_variants
        self.vcf_header = vcf_header

    def construct_variant_entry(self, proband_variant_data: ProbandCausativeVariant) -> list[str]:
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

    def construct_vcf_records(self):
        """Inserts spiked variant into correct position within VCF."""
        for variant in self.proband_causative_variants:
            variant = self.construct_variant_entry(variant)
            variant_entry_position = [
                                         i
                                         for i, val in enumerate(self.vcf_contents)
                                         if
                                         val.split("\t")[0] == variant[0] and int(val.split("\t")[1]) < int(variant[1])
                                     ][-1] + 1
            self.vcf_contents.insert(variant_entry_position, "\t".join(variant))

    def construct_header(self) -> list[str]:
        """Constructs the header of the VCF."""
        updated_vcf = []
        for line in self.vcf_contents:
            text = line.replace(
                self.vcf_header.sample_id,
                self.proband_causative_variants[0].proband_id,
            )
            updated_vcf.append(text)
        return updated_vcf

    def construct_vcf(self):
        self.construct_vcf_records()
        return self.construct_header()


class VcfWriter:
    def __init__(
            self,
            template_vcf_name: Path,
            vcf_contents: list[str],
            spiked_vcf_file_path: Path,
    ):
        self.template_vcf_name = template_vcf_name
        self.vcf_contents = vcf_contents
        self.spiked_vcf_file_path = spiked_vcf_file_path

    def copy_template_to_new_file(self):
        """Copies template vcf file to a new file."""
        shutil.copyfile(self.template_vcf_name, self.spiked_vcf_file_path)

    def write_gzip(self):
        """Writes gzipped vcf file."""
        encoded_contents = [line.encode() for line in self.vcf_contents]
        with gzip.open(self.spiked_vcf_file_path, "wb") as f:
            for line in encoded_contents:
                f.write(line)
        f.close()

    def write_uncompressed(self):
        """Writes an uncompressed vcf file."""
        with open(self.spiked_vcf_file_path, "w") as file:
            file.writelines(self.vcf_contents)
        file.close()

    def write_vcf_file(self):
        """Writes spiked vcf file."""
        self.copy_template_to_new_file()
        self.write_gzip() if is_gzipped(self.spiked_vcf_file_path) else self.write_uncompressed()


def add_vcf_path_to_phenopacket(
        phenopacket_path: Path, phenopacket: Phenopacket or Family, spiked_vcf_path: Path, vcf_header: VcfHeader
) -> None:
    """Writes the created vcf path to phenopacket"""
    updated_phenopacket = PhenopacketRebuilder(phenopacket).add_spiked_vcf_path(spiked_vcf_path, vcf_header.assembly)
    write_phenopacket(updated_phenopacket, phenopacket_path)


def spike_vcf(phenopacket: Phenopacket or Family, output_dir: Path, template_vcf_path: Path = None,
              vcf_dir: Path = None):
    """Spikes VCF records with variants."""
    # this is a separate function to a click command as it will fail if annotated with click annotations
    # and referenced from another click command
    if template_vcf_path is None and vcf_dir is None:
        raise InputError("Either a template_vcf or vcf_dir must be specified")
    try:
        output_dir.mkdir()
        info_log.info(f" Created a directory {output_dir}")
    except FileExistsError:
        pass
    phenopacket_causative_variants = PhenopacketUtil(phenopacket).causative_variants()
    chosen_template_vcf = VcfPicker(template_vcf_path, vcf_dir).pick_file()
    vcf_contents = read_vcf(chosen_template_vcf)
    vcf_header = VcfHeaderParser(vcf_contents).parse_vcf_header()
    check_variant_assembly(phenopacket_causative_variants, vcf_header, phenopacket)
    return template_vcf_path, vcf_header.assembly, VcfSpiker(vcf_contents, phenopacket_causative_variants,
                                                             vcf_header).construct_vcf()


@click.command()
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="FILE",
    required=True,
    help="Path to phenopacket file",
    type=Path,
)
@click.option(
    "--template-vcf-path",
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
        phenopacket_path: Path, output_dir: Path, template_vcf_path: Path = None, vcf_dir: Path = None
):
    """Spikes variants into a template VCF file for a singular phenopacket."""
    phenopacket = phenopacket_reader(phenopacket_path)
    template_vcf_path, vcf_assembly, spiked_vcf = spike_vcf(phenopacket, output_dir, template_vcf_path, vcf_dir)
    spiked_vcf_path = output_dir.joinpath(phenopacket_path.name.replace(".json", ".vcf.gz")) if is_gzipped(
        template_vcf_path) else output_dir.joinpath(phenopacket_path.name.replace(".json", ".vcf"))
    VcfWriter(
        template_vcf_path, spiked_vcf, spiked_vcf_path).write_vcf_file()
    updated_phenopacket = PhenopacketRebuilder(phenopacket).add_spiked_vcf_path(spiked_vcf_path, vcf_assembly)
    write_phenopacket(updated_phenopacket, phenopacket_path)


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
    "--template-vcf-path",
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
        template_vcf_path: Path = None,
        vcf_dir: Path = None,
):
    """Spikes variants into a template VCF file for a directory of phenopackets."""
    for phenopacket_path in files_with_suffix(phenopacket_dir, ".json"):
        phenopacket = phenopacket_reader(phenopacket_path)
        template_vcf_path, vcf_assembly, spiked_vcf = spike_vcf(phenopacket, output_dir, template_vcf_path, vcf_dir)
        spiked_vcf_path = output_dir.joinpath(phenopacket_path.name.replace(".json", ".vcf.gz")) if is_gzipped(
            template_vcf_path) else output_dir.joinpath(phenopacket_path.name.replace(".json", ".vcf"))
        VcfWriter(
            template_vcf_path, spiked_vcf, spiked_vcf_path).write_vcf_file()
        updated_phenopacket = PhenopacketRebuilder(phenopacket).add_spiked_vcf_path(spiked_vcf_path, vcf_assembly)
        write_phenopacket(updated_phenopacket, phenopacket_path)

    # or made a lambda one-liner for maximum wtf...
    # [spike_vcf(path, output_dir, template_vcf, vcf_dir) for path in phenopacket_dir.iterdir() if path.suffix ==
    # ".json"]
