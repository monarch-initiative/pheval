"""PhEval utils Command Line Interface"""

from pathlib import Path

import click

from pheval.utils.utils import semsim_randomisation


@click.command()
@click.option(
    "--semsim-left",
    "-L",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
)
@click.option(
    "--semsim-right",
    "-R",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
)
@click.option(
    "--output",
    "-O",
    metavar="FILE",
    required=True,
    help="",
)
@click.option(
    "--scramble-factor",
    "-S",
    metavar=float,
    default=0.5,
    help="",
)
@click.option(
    "--times",
    "-T",
    metavar=int,
    default=1,
    help="",
)
def scramble_semsim(
    semsim_left: Path, semsim_right: Path, output: Path, scramble_factor: float, times: int
):
    """scramble_semsim

    Args:
        semsim_left (Path): [description]
        semsim_right (Path): [description]
        output (Path): [description]
        scramble_factor (float): [description]
        times (int): [description]
    """
    semsim_randomisation(semsim_left, semsim_right, output, scramble_factor, times)


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
