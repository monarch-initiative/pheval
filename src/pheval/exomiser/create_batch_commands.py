#!/usr/bin/python

import os
import click
import tempfile
from pathlib import Path

from pheval.prepare.custom_exceptions import MutuallyExclusiveOptionError
from pheval.utils.file_utils import DirectoryFiles
from pheval.utils.phenopacket_utils import PhenopacketReader


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


class CreateBatchFiles:
    def __init__(self, analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file):
        self.analysis = analysis
        self.ppackets = ppackets
        self.vcf = vcf
        self.genome_assembly = genome_assembly
        self.output_opt = output_opt
        self.output_options_file = output_options_file

    def create_temp_commands(self) -> str:
        temp = tempfile.NamedTemporaryFile(delete=False)
        commands_writer = CommandsWriter(temp.name)
        if self.output_opt is None:
            arg_dict = dict((z[0], list(z[1:])) for z in zip(self.ppackets, self.vcf, self.genome_assembly))
            if self.output_options_file is None:
                for key, value in arg_dict.items():
                    commands_writer.write_command(self.analysis, key, value[0], value[1])
            if self.output_options_file is not None:
                for key, value in arg_dict.items():
                    commands_writer.write_command_output_options(self.analysis, key, value[0], value[1],
                                                                 self.output_options_file)
        if self.output_opt is not None:
            arg_dict = dict(
                (z[0], list(z[1:])) for z in zip(self.ppackets, self.vcf, self.genome_assembly, self.output_opt))
            for key, value in arg_dict.items():
                commands_writer.write_command_output_options(self.analysis, key, value[0], value[1], value[2])
        commands_writer.close()
        return temp.name

    def create_all_commands(self, file_name):
        commands_writer = CommandsWriter(file_name)
        if self.output_opt is None:
            arg_dict = dict((z[0], list(z[1:])) for z in zip(self.ppackets, self.vcf, self.genome_assembly))
            if self.output_options_file is None:
                for key, value in arg_dict.items():
                    commands_writer.write_command(self.analysis, key, value[0], value[1])
            if self.output_options_file is not None:
                for key, value in arg_dict.items():
                    commands_writer.write_command_output_options(self.analysis, key, value[0], value[1],
                                                                 self.output_options_file)
        if self.output_opt is not None:
            arg_dict = dict(
                (z[0], list(z[1:])) for z in zip(self.ppackets, self.vcf, self.genome_assembly, self.output_opt))
            for key, value in arg_dict.items():
                commands_writer.write_command_output_options(self.analysis, key, value[0], value[1], value[2])
        commands_writer.close()

    @staticmethod
    def create_split_batch(max_jobs, prefix, temp_name):
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
              required=True, metavar='FILE', type=Path,
              help="Path to the analysis .yml file.")
@click.option("--phenopacket-dir", "-p",
              required=True, metavar='PATH', type=Path,
              help="Path to phenopackets.")
@click.option("--vcf-dir", "-v",
              required=True, metavar='PATH', type=Path,
              help="Path to VCF files.")
@click.option("--batch-prefix", "-b",
              required=True, metavar='TEXT',
              help="Prefix of generated batch files.")
@click.option("--max-jobs", "-j",
              required=False, metavar='<int>', default=0, show_default=True,
              help="Number of jobs in each file.")
@click.option("--output-options-dir", "-O", cls=MutuallyExclusiveOptionError,
              mutually_exclusive=["output_options_file"],
              required=False, metavar='PATH', type=Path,
              help="Path to the output options directory. ")
@click.option("--output-options-file", "-o", cls=MutuallyExclusiveOptionError,
              mutually_exclusive=["output_options_dir"],
              required=False, metavar='FILE', type=Path,
              help="Path to the output options file. ")
def prepare_exomiser_batch(analysis: Path, phenopacket_dir: Path, vcf_dir: Path, batch_prefix, max_jobs,
                           output_options_dir: Path = None, output_options_file: Path = None):
    """ Generate Exomiser batch files. """
    vcf, genome_assembly = [], []
    ppackets = DirectoryFiles(phenopacket_dir, ".json").obtain_files_full_path()
    for ppacket in ppackets:
        vcf_file, assembly = PhenopacketReader(ppacket).vcf_file_data(vcf_dir)
        vcf.append(vcf_file)
        genome_assembly.append(assembly)
    if output_options_dir is not None:
        output_opt = DirectoryFiles(output_options_dir, "").obtain_all_files()
    else:
        output_opt = None
    if max_jobs != 0:
        temp_name = CreateBatchFiles(analysis, ppackets, vcf, genome_assembly, output_opt,
                                     output_options_file).create_temp_commands()
        CreateBatchFiles.create_split_batch(max_jobs, batch_prefix, temp_name)
    else:
        filename = batch_prefix + "-exomiser-batch.txt"
        CreateBatchFiles(analysis, ppackets, vcf, genome_assembly, output_opt, output_options_file).create_all_commands(
            filename)
