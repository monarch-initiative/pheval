#!/usr/bin/python
import click
import shutil
import os
import random
import logging
from pathlib import Path
from dataclasses import dataclass
from pheval.utils.file_utils import DirectoryFiles
from pheval.utils.phenopacket_utils import PhenopacketReader, CausativeVariant, IncompatibleGenomeAssemblyError
from pheval.prepare.custom_exceptions import MutuallyExclusiveOptionError, InputError

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
fh = logging.FileHandler(r'pheval_logger.txt')
logger.addHandler(ch)
logger.addHandler(fh)

genome_assemblies = {
    "GRCh38": {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
               "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
               "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
               '19': 58617616, "20": 64444167, '21': 46709983, '22': 50818468},
    "GRCh37": {"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067,
               "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
               "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248,
               '19': 59128983, "20": 63025520, '21': 48129895, '22': 51304566}}


@dataclass
class VcfHeader:
    sample_id: str
    assembly: str
    chr_status: bool


class VcfPicker:
    def __init__(self, template_vcf: Path, vcf_dir: Path):
        if template_vcf is not None:
            self.vcf_file = template_vcf
        if vcf_dir is not None:
            self.vcf_file = random.choice(os.listdir(vcf_dir))


class VcfParser:
    def __init__(self, template_vcf_file: str):
        self.vcf_file = open(template_vcf_file)

    def parse_header(self) -> VcfHeader:
        assembly_dict = {}
        sample_id = ""
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
            if line.startswith("#CHROM"):
                sample_id = line.split("\t")[9].rstrip()
        assembly = [k for k, v in genome_assemblies.items() if v == assembly_dict][0]
        self.vcf_file.close()
        return VcfHeader(sample_id, assembly, chr_status)


class ProbandVariantChecker:
    def __init__(self, list_of_variant_proband_data: list[CausativeVariant], vcf_header: VcfHeader):
        self.list_of_variant_proband_data = list_of_variant_proband_data
        self.vcf_header = vcf_header
        self.compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]

    def check_variant_assembly(self) -> list[CausativeVariant]:
        incorrect_variants = []
        for variant in self.list_of_variant_proband_data:
            if variant.assembly not in self.compatible_genome_assembly:
                raise IncompatibleGenomeAssemblyError(variant.assembly, variant.phenopacket)
            if variant.assembly != self.vcf_header.assembly:
                incorrect_variants.append(variant)
        return incorrect_variants


class VcfSpiker:
    def __init__(self, phenopacket: str, vcf_file: str, output_dir: Path,
                 list_of_variant_proband_data: list[CausativeVariant], vcf_header: VcfHeader):
        self.phenopacket = phenopacket
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.list_of_proband_variant_data = list_of_variant_proband_data
        self.vcf_header = vcf_header
        self.proband_vcf_file_name = os.path.join(self.output_dir, self.phenopacket.replace(".json", ".vcf"))

    def copy_template_vcf(self) -> None:
        shutil.copyfile(self.vcf_file, self.proband_vcf_file_name)

    def construct_variant(self, variant: CausativeVariant) -> list[str]:
        genotype_codes = {"hemizygous": "0/1", "homozygous": "1/1", "heterozygous": "0/1",
                          "compound heterozygous": "0/1"}
        if self.vcf_header.chr_status is True:
            variant.chrom = "chr" + variant.chrom
        variant_data = [variant.chrom, str(variant.pos), ".", variant.ref, variant.alt,
                        "100", "PASS", "SPIKED_VARIANT_" + variant.genotype.upper(), "GT:AD:DP:GQ:PL",
                        genotype_codes[variant.genotype.lower()] + ":0,2:2:12:180,12,0" + "\n"]
        return variant_data

    def construct_vcf_records(self) -> list[str]:
        vcf_file = open(self.vcf_file)
        vcf_contents = vcf_file.readlines()
        vcf_file.close()
        for variant in self.list_of_proband_variant_data:
            variant = self.construct_variant(variant)
            locs = [i for i, val in enumerate(vcf_contents) if
                    val.split("\t")[0] == variant[0] and int(val.split("\t")[1]) < int(variant[1])][-1] + 1
            vcf_contents.insert(locs, "\t".join(variant))
        return vcf_contents

    def construct_header(self) -> list[str]:
        updated_vcf = []
        updated_records = self.construct_vcf_records()
        for line in updated_records:
            text = line.replace(self.vcf_header.sample_id, self.list_of_proband_variant_data[0].proband_id)
            updated_vcf.append(text)
        return updated_vcf

    def write_vcf(self) -> None:
        self.copy_template_vcf()
        new_vcf = self.construct_header()
        with open(self.proband_vcf_file_name, 'w') as file:
            file.writelines(new_vcf)
        file.close()


@click.command()
@click.option("--phenopacket", "-p", metavar='FILE', required=True, help="Path to phenopacket file", type=Path)
@click.option("--template-vcf", "-t", cls=MutuallyExclusiveOptionError, metavar="FILE", required=False,
              help="Template VCF file", mutually_exclusive=["vcf_dir"], type=Path)
@click.option("--vcf-dir", "-v", cls=MutuallyExclusiveOptionError, metavar="PATH",
              help="Directory containing template VCF files", mutually_exclusive=["template_vcf"], type=Path)
@click.option("--output-dir", "-O", metavar="PATH", required=True, help="Path for creation of output directory",
              default="vcf", type=Path)
def create_spiked_vcf(phenopacket: Path, output_dir: Path, template_vcf=None, vcf_dir=None):
    print(phenopacket, output_dir, template_vcf, vcf_dir)
    try:
        # TODO update with path api
        # os.mkdir(os.path.join(output_dir, ''))
        output_dir.mkdir()
        logger.info(f" Created a directory {output_dir}")
    except FileExistsError:
        pass
    phenopacket_proband_variant_data = PhenopacketReader(phenopacket).causative_variants()
    # TODO: check that chosen template is what is expected (VcfPicker)
    chosen_template_vcf = VcfPicker(template_vcf, vcf_dir)
    vcf_header = VcfParser(chosen_template_vcf.vcf_file).parse_header()
    incompatible_variants = ProbandVariantChecker(phenopacket_proband_variant_data,
                                                  vcf_header).check_variant_assembly()
    if len(incompatible_variants) != 0:
        for incompatible_variant in incompatible_variants:
            logger.error(
                f' Skipping... Proband variant does not match Human Genome Build of VCF: {phenopacket.absolute().name}. '
                f' Variant Assembly -> {incompatible_variant.assembly} Expected: {vcf_header.assembly}')
            raise IncompatibleGenomeAssemblyError(assembly=incompatible_variant.assembly,
                                                  phenopacket=phenopacket.absolute().name)
    VcfSpiker(phenopacket.name, chosen_template_vcf.vcf_file, output_dir, phenopacket_proband_variant_data,
              vcf_header).write_vcf()


@click.command()
@click.option("--phenopacket-dir", "-p", metavar='PATH', required=True, help="Path to phenopackets directory", type=Path)
@click.option("--template-vcf", "-t", cls=MutuallyExclusiveOptionError, metavar="PATH", required=False,
              help="Template VCF file", mutually_exclusive=["vcf_dir"], type=Path)
@click.option("--vcf-dir", "-v", cls=MutuallyExclusiveOptionError, metavar="PATH",
              help="Directory containing template VCF files", mutually_exclusive=["template_vcf"], type=Path)
@click.option("--output-dir", "-O", metavar="PATH", required=True, help="Path for creation of output directory",
              default="vcf", type=Path)
def create_spiked_vcfs(phenopacket_dir: Path, output_dir: Path, template_vcf: Path = None, vcf_dir: Path = None):
    """ Spikes variants into a template VCF file. """
    if template_vcf is None and vcf_dir is None:
        raise InputError("VCF")
    try:
        os.mkdir(os.path.join(output_dir, ''))
    except FileExistsError:
        pass
    phenopackets = [path for path in phenopacket_dir.iterdir() if path.suffix == ".json"]
    for phenopacket in phenopackets:
        # phenopacket_full_path = os.path.join(phenopacket_dir, phenopacket)
        create_spiked_vcf(phenopacket, output_dir, template_vcf, vcf_dir)
        # phenopacket_proband_variant_data = PhenopacketReader(phenopacket_full_path).causative_variants()
        # chosen_template_vcf = VcfPicker(template_vcf, vcf_dir)
        # vcf_header = VcfParser(chosen_template_vcf.vcf_file).parse_header()
        # incompatible_variants = ProbandVariantChecker(phenopacket_proband_variant_data,
        #                                               vcf_header).check_variant_assembly()
        # if len(incompatible_variants) != 0:
        #     for incompatible_variant in incompatible_variants:
        #         logger.error(f' Skipping... Proband variant does not match Human Genome Build of VCF: {phenopacket}. '
        #                      f' Variant Assembly -> {incompatible_variant.assembly} Expected: {vcf_header.assembly}')
        #     continue
        # VcfSpiker(phenopacket, chosen_template_vcf.vcf_file, output_dir,
        #           phenopacket_proband_variant_data, vcf_header).write_vcf()
