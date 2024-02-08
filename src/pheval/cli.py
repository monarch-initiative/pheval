"""PhEval CLI Module """

import logging

import click

from .cli_pheval import run
from .cli_pheval_utils import (
    benchmark,
    benchmark_comparison,
    create_spiked_vcfs_command,
    generate_stats_plot,
    scramble_phenopackets_command,
    semsim_scramble_command,
    semsim_to_exomiserdb_command,
    update_phenopackets_command,
)

info_log = logging.getLogger("info")


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet")
def main(verbose=1, quiet=False) -> None:
    """main CLI method for PhEval

    Args:
        verbose (int, optional): Verbose flag.
        quiet (bool, optional): Queit Flag.
    """
    if verbose >= 2:
        info_log.setLevel(level=logging.DEBUG)
    elif verbose == 1:
        info_log.setLevel(level=logging.INFO)
    else:
        info_log.setLevel(level=logging.WARNING)
    if quiet:
        info_log.setLevel(level=logging.ERROR)


@click.group()
def pheval():
    """pheval"""


pheval.add_command(run)


@click.group()
def pheval_utils():
    """pheval_utils"""


pheval_utils.add_command(semsim_scramble_command)
pheval_utils.add_command(scramble_phenopackets_command)
pheval_utils.add_command(update_phenopackets_command)
pheval_utils.add_command(create_spiked_vcfs_command)
pheval_utils.add_command(benchmark)
pheval_utils.add_command(benchmark_comparison)
pheval_utils.add_command(semsim_to_exomiserdb_command)
pheval_utils.add_command(generate_stats_plot)

if __name__ == "__main__":
    main()
