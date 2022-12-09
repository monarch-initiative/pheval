"""PhEval utils Command Line Interface"""

import logging
from typing import List

import click

import pheval.utils as utils

info_log = logging.getLogger("info")


@click.command()
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
)
def scramble_semsim():
    """scramble_semsim"""
    print("running pheval_utils::scramble_semsim command")


@click.command()
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the phenopacket to be spiked.",
)
def scramble_phenopacket():
    """scramble_phenopacket"""
    print("running pheval_utils::scramble_phenopacket command")


@click.command()
@click.option(
    "--left",
    "-l",
    required=True,
    metavar="FILE",
    help="Path of first file that will be compared",
)
@click.option(
    "--right",
    "-r",
    required=True,
    metavar="FILE",
    help="Path of second file that will be compared",
)
@click.option(
    "--join-columns",
    "-j",
    required=False,
    multiple=True,
    default=["GENE_SYMBOL"],
    metavar="LIST",
    help="Columns that will be used in comparison",
)
@click.option(
    "--comparison-columns",
    "-c",
    required=True,
    multiple=True,
    metavar="LIST",
    help="Columns that will be used in comparison",
)
@click.option(
    "--sort-columns",
    "-s",
    required=True,
    multiple=True,
    metavar="LIST",
    default=["P-VALUE"],
    help="Columns that will be used for ranking purposes",
)
@click.option(
    "--output",
    "-o",
    required=True,
    metavar="FILE",
    help="Path to the result output",
)
def compare_semsim(
    left: click.Path,
    right: click.Path,
    output: click.Path,
    join_columns: List[str],
    comparison_columns: List[str],
    sort_columns: List[str],
):
    """compare_semsim"""
    try:
        dataframe = utils.calc_semsim(left, right, join_columns, comparison_columns, sort_columns)
        utils.ensure_file_exists(left, right)
        dataframe.to_csv(output, index=False)
        info_log.info("done")
    except (ValueError, FileNotFoundError) as err:
        info_log.error("Error: %s", err)
