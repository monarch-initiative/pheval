#!/usr/bin/python

import os
import yaml
import click
import tempfile

class IncompatibleGenomeAssemblyError(Exception):
    """ Exception raised for incompatible genome assembly."""
    def __init__(self, assembly, phenopacket, message="Incompatible Genome Assembly"):
        self.assembly: str = assembly
        self.phenopacket: str = phenopacket
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} -> {self.assembly} in {self.phenopacket}'


class IncorrectFileFormatError(Exception):
    def __init__(self, file, expectation, message="Incorrect File Type"):
        self.file: str = file
        self.expectation: str = expectation
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message} -> {self.file} (expected {self.expectation})'


class CommandsWriter:
    def __init__(self, file: str):
        self.file = open(file, 'a')

    def write_command(self, analysis, sample, vcf, assembly):
        try:
            self.file.write("--analysis " + analysis + " --sample " + sample +
                            " --vcf " + vcf + " --assembly " + assembly + "\n")
        except IOError:
            print("Error writing ", self.file)

    def write_command_output_options(self, analysis, sample, vcf, assembly, output_options):
        try:
            self.file.write("--analysis " + analysis + " --sample " + sample +
                            " --vcf " + vcf + " --assembly " + assembly + " --output " + output_options + "\n")
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)

def phenopacket_list(ppacket_dir) -> list:
    """ Returns a list of the full path to phenopackets found in the directory specified. """
    ppackets = []
    ppacket_directory = os.fsencode(ppacket_dir)
    for file in os.listdir(ppacket_directory):
        filename = os.fsdecode(file)
        ppackets.append(ppacket_dir + filename)
    ppackets.sort()
    return ppackets


def vcf_list(ppackets, vcf_dir) -> tuple[list, list]:
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
                    if not vcf_name.endswith(".vcf") and not vcf_name.endswith(".vcf.gz"):
                        raise IncorrectFileFormatError(vcf_name, ".vcf or .vcf.gz file")                    
                    vcf.append(vcf_dir + vcf_name)
                    assembly = file["fileAttributes"]["genomeAssembly"]
                    if assembly not in compatible_genome_assembly:
                        raise IncompatibleGenomeAssemblyError(assembly, ppacket)                    
                    genome_assembly.append(assembly)
        ppacket_file.close()
    return vcf, genome_assembly


def output_options_list(output_options_dir) -> list:
    """ Returns a list of the full path to output options files. """
    output_opt = []
    output_options_directory_ = os.fsencode(output_options_dir)
    for file in os.listdir(output_options_directory_):
        filename = os.fsdecode(file)
        output_opt.append(output_options_dir + filename)
    output_opt.sort()
    return output_opt


def create_temp_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file):
    """ Writes a temporary file of all commands. """
    temp = tempfile.NamedTemporaryFile(delete=False)
    commands_writer = CommandsWriter(temp.name)
    if output_opt is None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly))
        if output_options_file is None:
            for key, value in arg_dict.items():
                commands_writer.write_command(analysis, key, value[0], value[1])            
        if output_options_file is not None:
            for key, value in arg_dict.items():
                commands_writer.write_command_output_options(analysis, key, value[0], value[1], output_options_file)
    if output_opt is not None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly, output_opt))
        for key, value in arg_dict.items():
            commands_writer.write_command_output_options(analysis, key, value[0], value[1], value[2])
    commands_writer.close()
    return temp.name


def create_all_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file, filename):
    """ Writes a temporary file of all commands. """
    commands_writer = CommandsWriter(filename)
    if output_opt is None:
        arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly))
        if output_options_file is None:
            for key, value in arg_dict.items():
                commands_writer.write_command(analysis, key, value[0], value[1])
        if output_options_file is not None:
            for key, value in arg_dict.items():
                commands_writer.write_command_output_options(analysis, key, value[0], value[1], output_options_file)
        if output_opt is not None:
            arg_dict = dict((z[0], list(z[1:])) for z in zip(ppackets, vcf, genome_assembly, output_opt))
            for key, value in arg_dict.items():
                commands_writer.write_command_output_options(analysis, key, value[0], value[1], value[2])
    commands_writer.close()
    
    
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
              required=False, metavar='<int>', default=False, show_default=True,
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
    if max_jobs:
        temp_name = create_temp_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file)
        create_split_batch(max_jobs, batch_prefix, temp_name)
    else:
        filename = batch_prefix + "-exomiser-batch.txt"
        create_all_commands(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file, filename)


if __name__ == '__main__':
    prepare_exomiser_batch()
