#!/usr/bin/python
import click
import shutil
import json
import os
import random
import logging
from pathlib import Path
from click import Option, UsageError
from phenopackets import Phenopacket
from dataclasses import dataclass
from google.protobuf.json_format import Parse
from pheval_benchmark.assess_prioritisation import directory_files
import pheval_benchmark.create_batch_commands

logging.basicConfig(filename="log.txt", filemode='w')

genome_assemblies = {
    "GRCh38": {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
               "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
               "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
               '19': 58617616, "20": 64444167, '21': 46709983, '22': 50818468},
    "GRCh37": {"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067,
               "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
               "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248,
               '19': 59128983, "20": 63025520, '21': 48129895, '22': 51304566}}


class InputError(Exception):
    """ Exception raised for missing required inputs."""

    def __init__(self, file, message="Missing required input"):
        self.file: str = file
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} -> {self.file} '


class MutuallyExclusiveOptionError(Option):
    """ Exception raised for when """

    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop('mutually_exclusive', []))
        help_ = kwargs.get('help', '')
        if self.mutually_exclusive:
            ex_str = ', '.join(self.mutually_exclusive)
            kwargs['help'] = help_ + (
                    ' NOTE: This argument is mutually exclusive with '
                    ' arguments: [' + ex_str + '].'
            )
        super(MutuallyExclusiveOptionError, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(
                    self.name,
                    ', '.join(self.mutually_exclusive)
                )
            )

        return super(MutuallyExclusiveOptionError, self).handle_parse_result(
            ctx,
            opts,
            args
        )


@dataclass
class ProbandData:
    proband_id: str
    assembly: str
    chrom: str
    pos: str
    ref: str
    alt: str
    genotype: str


class PhenopacketReader:

    def __init__(self, file: str):
        file_name = Path(file).name
        self.vcf_file_name = file_name.replace(".json", ".vcf")
        with open(file) as f:
            ppacket = json.load(f)
        f.close()
        self.ppacket = Parse(json.dumps(ppacket), Phenopacket())

    def extract_data(self):
        all_variants = []
        for i in self.ppacket.interpretations:
            for g in i.diagnosis.genomic_interpretations:
                record = g.variant_interpretation.variation_descriptor.vcf_record
                genotype = g.variant_interpretation.variation_descriptor.allelic_state
                variant = ProbandData(self.ppacket.subject.id, record.genome_assembly, "chr" + record.chrom, record.pos,
                                      record.ref, record.alt, genotype.label)
                all_variants.append(variant)
            return all_variants


class VcfReader:
    def __init__(self, file: str):
        self.file = open(file)

    def extract_sample_id(self):
        template_vcf_sample_id: str = ""
        for line in self.file:
            if line.startswith("#CHROM"):
                template_vcf_sample_id = line.split("\t")[9].rstrip()
                break
        self.file.seek(0)
        return template_vcf_sample_id

    def close_file(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


class VcfWriter(VcfReader):
    def __init__(self, copied_vcf: str, list_of_proband_data):
        super().__init__(copied_vcf)
        self.sample_id = self.extract_sample_id()
        self.text = self.file.readlines()
        self.variant_info = list_of_proband_data
        self.vcf_name = copied_vcf
        self.updated_vcf = []
        self.chr_status = False

    def check_variants(self):
        incorrect_variants = []
        compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
        for variant in self.variant_info:
            if variant.assembly not in compatible_genome_assembly:
                raise pheval_benchmark.create_batch_commands.IncompatibleGenomeAssemblyError(variant.assembly,
                                                                                             self.vcf_name.split("/")[
                                                                                                 1].replace(
                                                                                                 "vcf", "json"))
            if variant.assembly == "hg19":
                variant.assembly = "GRCh37"
            if variant.assembly == "hg38":
                variant.assembly = "GRCh38"
            for line in self.text:
                if line.startswith("##contig=<ID"):
                    line = line.split(",")
                    if "chr" in line[0]:
                        self.chr_status = True
                    chromosome = line[0].split("=")[2]
                    length = line[1].split("=")[1]
                    try:
                        if self.chr_status:
                            chromosome = chromosome.replace("chr", "")
                        if genome_assemblies[variant.assembly][chromosome] == int(length):
                            pass
                        else:
                            incorrect_variants.append(variant)
                    except KeyError:
                        pass
        return incorrect_variants

    def construct_variant(self, variant):
        genotype_codes = {"hemizygous": "0/1", "homozygous": "1/1", "heterozygous": "0/1",
                          "compound heterozygous": "0/1"}
        if self.chr_status is True:
            variant.chrom = "chr" + variant.chrom
        variant_data = [variant.chrom, str(variant.pos), ".", variant.ref, variant.alt,
                        "100", "PASS", "SPIKED_VARIANT_" + variant.genotype.upper(), "GT:AD:DP:GQ:PL",
                        genotype_codes[variant.genotype.lower()] + ":0,2:2:12:180,12,0" + "\n"]
        return variant_data

    def construct_records(self):
        for variant in self.variant_info:
            variant = self.construct_variant(variant)
            locs = [i for i, val in enumerate(self.text) if
                    val.split("\t")[0] == variant[0] and int(val.split("\t")[1]) < int(variant[1])][-1] + 1
            self.text.insert(locs, "\t".join(variant))
        return self.text

    def construct_header(self):
        updated_records = self.construct_records()
        for line in updated_records:
            text = line.replace(self.sample_id, self.variant_info[0].proband_id)
            self.updated_vcf.append(text)
        return self.updated_vcf

    def write_vcf(self):
        new_vcf = self.construct_header()
        self.close_file()
        with open(self.vcf_name, 'w') as file:
            file.writelines(new_vcf)
        file.close()


def create_basic_vcf(phenopacket_dir, template_vcf, vcf_dir, output_dir) -> list:
    """ Copies a template VCF file to a new output directory with the same file prefix as the phenopacket file."""
    try:
        os.mkdir(os.path.join(output_dir, ''))
    except FileExistsError:
        pass
    phenopackets = directory_files(phenopacket_dir)
    for phenopacket in phenopackets:
        vcf_file_name = phenopacket.replace(".json", ".vcf")
        if template_vcf is not None:
            shutil.copyfile(template_vcf, os.path.join(output_dir, vcf_file_name))
        if vcf_dir is not None:
            file = random.choice(os.listdir(vcf_dir))
            shutil.copyfile(file, os.path.join(output_dir, vcf_file_name))
    return phenopackets


@click.command()
@click.option("--phenopacket-dir", "-p", metavar='PATH', required=True, help="Path to phenopackets directory")
@click.option("--template-vcf", "-t", cls=MutuallyExclusiveOptionError, metavar="FILE", required=False,
              help="Template VCF file", mutually_exclusive=["vcf_dir"])
@click.option("--vcf-dir", "-v", cls=MutuallyExclusiveOptionError, metavar="PATH",
              help="Directory containing template VCF files", mutually_exclusive=["template_vcf"])
@click.option("--output-dir", "-O", metavar="PATH", required=True, help="Path for creation of output directory",
              default="vcf")
def spike_vcf(phenopacket_dir, output_dir, template_vcf=None, vcf_dir=None):
    """ Spikes variants into a template VCF file. """
    if template_vcf is None and vcf_dir is None:
        raise InputError("VCF")
    phenopackets = create_basic_vcf(phenopacket_dir, template_vcf, vcf_dir, output_dir)
    for phenopacket in phenopackets:
        phenopacket_full_path = os.path.join(phenopacket_dir, phenopacket)
        vcf_full_path = os.path.join(output_dir, phenopacket.replace(".json", ".vcf"))
        phenopacket_data = PhenopacketReader(phenopacket_full_path)
        phenopacket_proband_data = phenopacket_data.extract_data()
        incorrect_variants = VcfWriter(vcf_full_path, phenopacket_proband_data).check_variants()
        if incorrect_variants:
            logging.error(f'Proband variant does not match Human Genome Build of VCF.: {phenopacket}')
            continue
        VcfWriter(vcf_full_path, phenopacket_proband_data).write_vcf()
