"""
Monarch Initiative
"""

import click


@click.command()
@click.option(
    "--analysis",
    "-a",
    required=True,
    metavar="FILE",
    help="Path to the analysis .yml file.",
)
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    metavar="PATH",
    help="Path to phenopackets.",
)
@click.option("--vcfs", "-v", required=True, metavar="PATH", help="Path to VCF files.")
@click.option(
    "--batch-prefix",
    "-b",
    required=True,
    metavar="TEXT",
    help="Prefix of generated batch files.",
)
@click.option(
    "--max-jobs",
    "-j",
    required=False,
    metavar="<int>",
    default=False,
    show_default=True,
    help="Number of jobs in each file.",
)
@click.option(
    "--output-options-dir",
    "-O",
    required=False,
    metavar="PATH",
    help="Path to the output options directory. ",
)
@click.option(
    "--output-options-file",
    "-o",
    required=False,
    metavar="FILE",
    help="Path to the output options file. ",
)
def run(
    analysis,
    phenopacket_dir,
    vcfs,
    batch_prefix,
    max_jobs,
    output_options_dir,
    output_options_file,
):
    print("running pheval::run command")


def prepare_run(
    prepare_method: str
):
    print("running pheval::run command")
    # see https://stackoverflow.com/questions/67631/how-do-i-import-a-module-given-the-full-path
    module = load_module_from_string(prepare_method)
    runner = module.Runner()
    runner.run()