#!/usr/bin/python
import click
import vcfpy
import shutil
import json
import os
import csv
from pheval_benchmark.assess_prioritisation import directory_files


class VariantWriter:
    def __init__(self, file: str):
        self.file = open(file, 'a')
        self.writer = csv.writer(self.file, delimiter='\t', lineterminator='\n')

    def write_variant(self, chrom, pos, ref, alt):
        try:
            self.writer.writerow([chrom, pos, ".", ref, alt, "100", "PASS", ".", "GT:AD:DP:GQ:PL",
                                  "0", "/1:0,2:2:12:180,12,0"])
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


def create_basic_vcf(phenopacket_dir, template_vcf, output_dir) -> list:
    """ Copies a template VCF file to a new output directory with the same file prefix as the phenopacket file."""
    os.mkdir(output_dir)
    phenopackets = directory_files(phenopacket_dir)
    for phenopacket in phenopackets:
        vcf_file_name = phenopacket.replace(".json", ".vcf")
        shutil.copyfile(template_vcf, output_dir + vcf_file_name)
    return phenopackets


def update_subject_id(phenopacket_full_path, template_vcf, vcf_full_path):
    with open(phenopacket_full_path) as f:
        data = json.load(f)
        subject_id = data["subject"]["id"]
    f.close()
    vcf_reader = vcfpy.Reader(open(template_vcf, 'r'))
    template_vcf_sample_id = vcf_reader.header.samples.names[0]
    with open(vcf_full_path, "r") as vcf_file:
        data = vcf_file.readlines()
    vcf_file.close()
    updated_vcf = []
    for entry in data:
        entry = entry.replace(template_vcf_sample_id, subject_id)
        updated_vcf.append(entry)
    with open(vcf_full_path, 'w') as vcf_file:
        vcf_file.writelines(updated_vcf)
    vcf_file.close()


def add_variants(phenopacket_full_path, vcf_full_path):
    with open(phenopacket_full_path) as f:
        data = json.load(f)
        for i in data["interpretations"]:
            for g in i["diagnosis"]["genomicInterpretations"]:
                chrom = "chr" + g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]["chrom"]
                pos = g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]["pos"]
                ref = g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]["ref"]
                alt = g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]["alt"]
                variant_writer = VariantWriter(vcf_full_path)
                variant_writer.write_variant(chrom, pos, ref, alt)
    variant_writer.close()
    f.close()


@click.command()
@click.option("--phenopacket-dir", "-p", metavar='PATH', required=True, help="Path to phenopackets directory")
@click.option("--template-vcf", "-t", metavar="FILE", required=True, help="Template VCF file")
@click.option("--output-dir", "-O", metavar="PATH", required=True, help="Path for creation of output directory",
              default="vcf/")
def spike_vcf(phenopacket_dir, template_vcf, output_dir):
    """ Spikes variants into a template VCF file. """
    phenopackets = create_basic_vcf(phenopacket_dir, template_vcf, output_dir)
    for phenopacket in phenopackets:
        phenopacket_full_path = phenopacket_dir + phenopacket
        vcf_full_path = output_dir + phenopacket.replace(".json", ".vcf")
        update_subject_id(phenopacket_full_path, template_vcf, vcf_full_path)
        add_variants(phenopacket_full_path, vcf_full_path)
