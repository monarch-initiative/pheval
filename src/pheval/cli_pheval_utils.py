"""PhEval utils Command Line Interface"""

from pathlib import Path

import click

from pheval.prepare.create_noisy_phenopackets import scramble_phenopackets
from pheval.prepare.custom_exceptions import InputError, MutuallyExclusiveOptionError
from pheval.utils.semsim_utils import percentage_diff, semsim_heatmap_plot


@click.command()
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
)
def scramble_semsim(input: Path):
    """scramble_semsim"""
    print("running pheval_utils::scramble_semsim command")


@click.command("scramble-phenopackets")
@click.option(
    "--phenopacket-path",
    "-p",
    metavar="PATH",
    help="Path to phenopacket.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_dir"],
)
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    help="Path to phenopackets directory.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--scramble-factor",
    "-s",
    metavar=float,
    required=True,
    default=0.5,
    show_default=True,
    help="Scramble factor for randomising phenopacket phenotypic profiles.",
    type=float,
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="noisy_phenopackets",
    type=Path,
)
def scramble_phenopackets_command(
    phenopacket_path: Path,
    phenopacket_dir: Path,
    scramble_factor: float,
    output_dir: Path,
):
    """Generate noisy phenopackets from existing ones."""
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    else:
        scramble_phenopackets(output_dir, phenopacket_path, phenopacket_dir, scramble_factor)


@click.command("semsim-comparison")
@click.option(
    "--semsim-left",
    "-L",
    required=True,
    metavar="FILE",
    help="Path to the first semantic similarity profile.",
)
@click.option(
    "--semsim-right",
    "-R",
    required=True,
    metavar="FILE",
    help="Path to the second semantic similarity profile.",
)
@click.option(
    "--score-column",
    "-s",
    required=True,
    type=click.Choice(
        ["jaccard_similarity", "dice_similarity", "phenodigm_score"], case_sensitive=False
    ),
    help="Score column that will be used in comparison",
)
@click.option(
    "--analysis",
    "-a",
    required=True,
    type=click.Choice(["heatmap", "percentage_diff"], case_sensitive=False),
    help="""There are two types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles""",
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    default="percentage_diff.semsim.tsv",
    help="Output path for the difference tsv. Defaults to percentage_diff.semsim.tsv",
)
def semsim_comparison(
    semsim_left: Path,
    semsim_right: Path,
    score_column: str,
    analysis: str,
    output: Path = "percentage_diff.semsim.tsv",
):
    """Compares two semantic similarity profiles

    Args:
        semsim-left (Path): File path of the first semantic similarity profile
        semsim-right (Path): File path of the second semantic similarity profile
        output (Path): Output path for the difference tsv. Defaults to "percentage_diff.semsim.tsv".
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        analysis (str): There are two types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles.
    """
    if analysis == "heatmap":
        return semsim_heatmap_plot(semsim_left, semsim_right, score_column)
    if analysis == "percentage_diff":
        percentage_diff(semsim_left, semsim_right, score_column, output)
