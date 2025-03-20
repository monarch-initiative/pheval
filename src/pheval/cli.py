"""PhEval CLI Module"""

import logging

import click

from pheval.utils.logger import get_logger, initialise_context

from .cli_pheval import run
from .cli_pheval_utils import (
    benchmark,
    create_spiked_vcfs_command,
    generate_plots,
    prepare_corpus_command,
    scramble_phenopackets_command,
    semsim_scramble_command,
    semsim_to_exomiserdb_command,
    update_phenopackets_command,
)

logger = get_logger()


@click.group()
@click.option("-v", "--verbose", count=True)
@click.option("-q", "--quiet", is_flag=True)
@click.pass_context
def main(ctx, verbose=1, quiet=False):
    """Main CLI method for PhEval."""
    initialise_context(ctx)

    if verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    if quiet:
        logger.setLevel(logging.ERROR)


@main.group()
@click.pass_context
def pheval(ctx):
    """pheval"""
    initialise_context(ctx)


@main.group()
@click.pass_context
def pheval_utils(ctx):
    """pheval_utils"""
    initialise_context(ctx)


pheval.add_command(run)

pheval_utils.add_command(semsim_scramble_command)
pheval_utils.add_command(scramble_phenopackets_command)
pheval_utils.add_command(update_phenopackets_command)
pheval_utils.add_command(create_spiked_vcfs_command)
pheval_utils.add_command(benchmark)
pheval_utils.add_command(semsim_to_exomiserdb_command)
pheval_utils.add_command(prepare_corpus_command)
pheval_utils.add_command(generate_plots)

if __name__ == "__main__":
    main()
