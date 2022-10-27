#!/usr/bin/python

import os
import yaml
import click
import tempfile


def phenopacket_list(ppacket_dir) -> list:
    """ Returns a list of the full path to phenopackets found in the directory specified. """
    ppackets = []
    ppacket_directory = os.fsencode(ppacket_dir)
    for file in os.listdir(ppacket_directory):
        filename = os.fsdecode(file)
        ppackets.append(ppacket_dir + filename)
    ppackets.sort()
    return ppackets


def vcf_list(ppackets, vcf_dir):
    """ Returns a list of the full path to vcf files and corresponding genome assembly
    retrieved from the phenopacket. """
    vcf = []
    genome_assembly = []
    compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
    for ppacket in ppackets:
        with open(ppacket) as ppacket_file:
            p = yaml.load(ppacket_file, Loader=yaml.FullLoader)
            files = p["files"]
            for file in files:
                if file["fileAttributes"]["fileFormat"].upper() == "VCF":
                    vcf_path = file["uri"]
                    if "/" in vcf_path:
                        vcf_name = vcf_path.rsplit('/')[-1]
                    if "/" not in vcf_path:
                        vcf_name = vcf_path
                    vcf.append(vcf_dir + vcf_name)
                    assembly = file["fileAttributes"]["genomeAssembly"]
                    genome_assembly.append(assembly)
        ppacket_file.close()
    return vcf, genome_assembly


def output_options_list(output_options_dir):
    """ Returns a list of the full path to output options files. """
    output_opt = []
    output_options_directory_ = os.fsencode(output_options_dir)
    for file in os.listdir(output_options_directory_):
        filename = os.fsdecode(file)
        output_opt.append(output_options_dir + filename)
    output_opt.sort()
    return output_opt


def create_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file):
    """ Writes a temporary file of all commands. """
    temp = tempfile.NamedTemporaryFile(delete=False)
    if output_opt is None and output_options_file is None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly))
        if output_opt is None:
            for key, value in arg_dict.items():
                with open(temp.name, "a") as outfile:
                    outfile.write("--analysis " + analysis + " --sample " + key +
                                  " --vcf " + value[0] + " --assembly " + value[1] + "\n")
                outfile.close()
    if output_options_file is not None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly))
        for key, value in arg_dict.items():
            with open(temp.name, "a") as outfile:
                outfile.write("--analysis " + analysis + " --sample " + key +
                              " --vcf " + value[0] + " --assembly " + value[
                                  1] + " --output " + output_options_file + "\n")
            outfile.close()
    if output_opt is not None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly, output_opt))
        for key, value in arg_dict.items():
            with open(temp.name, "a") as outfile:
                outfile.write("--analysis" + " /data/home/preset-exome-analysis.yml" + " --sample " + key +
                              " --vcf " + value[0] + " --assembly " + value[1] + " --output "
                              + value[2] + "\n")
            outfile.close()
    return temp.name


def create_split_batch(max_jobs, prefix, temp_name):
    """ Splits temporary commands file into smaller batch jobs and deletes temp file. """
    lines_per_file, f_name = max_jobs, 0
    splitfile = None
    with open(temp_name) as tmp_file:
        for lineno, line in enumerate(tmp_file):
            if lineno % lines_per_file == 0:
                f_name += 1
                if splitfile:
                    splitfile.close()
                split_filename = prefix + '-exomiser-batch-{}.txt'.format(f_name)
                splitfile = open(split_filename, "w")
            splitfile.write(line)
        if splitfile:
            splitfile.close()
    tmp_file.close()
    os.remove(temp_name)


@click.command()
@click.option("--analysis", "-a",
              required=True, metavar='FILE',
              help="Path to the analysis .yml file.")
@click.option("--phenopacket-dir", "-p",
              required=True, metavar='PATH',
              help="Path to phenopackets.")
@click.option("--vcfs", "-v",
              required=True, metavar='PATH',
              help="Path to VCF files.")
@click.option("--batch-prefix", "-b",
              required=True, metavar='TEXT',
              help="Prefix of generated batch files.")
@click.option("--max-jobs", "-j",
              required=False, metavar='<int>', default=50, show_default=True,
              help="Number of jobs in each file.")
@click.option("--output-options-dir", "-O",
              required=False, metavar='PATH',
              help="Path to the output options directory. ")
@click.option("--output-options-file", "-o",
              required=False, metavar='FILE',
              help="Path to the output options file. ")
def prepare_exomiser_batch(analysis, phenopacket_dir, vcfs, batch_prefix, max_jobs=50, output_options_dir=None,
                           output_options_file=None):
    """ Generate Exomiser batch files. """
    ppackets = phenopacket_list(phenopacket_dir)
    vcf, genome_assembly = vcf_list(ppackets, vcfs)
    if output_options_dir is not None:
        output_opt = output_options_list(output_options_dir)
    else:
        output_opt = None
    temp_name = create_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file)
    create_split_batch(max_jobs, batch_prefix, temp_name)


if __name__ == '__main__':
    prepare_exomiser_batch()
