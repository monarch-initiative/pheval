"""PhEval utils Command Line Interface"""

from pathlib import Path

import click

from pheval.prepare.create_noisy_phenopackets import create_scrambled_phenopackets
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
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    required=True,
    help="Path to phenopackets directory",
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
    "--output-file-suffix",
    "-o",
    metavar="<str>",
    required=True,
    help="Suffix to append to output file",
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="noisy_phenopackets",
)
def scramble_phenopackets_command(
    phenopacket_dir: Path,
    scramble_factor: float,
    output_file_suffix: str,
    output_dir: Path,
):
    """Generate noisy phenopackets from existing ones."""
    create_scrambled_phenopackets(output_dir, output_file_suffix, phenopacket_dir, scramble_factor)


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
    type=str,
    help="Score column that will be used in comparison (e.g jaccard_similarity)",
)
@click.option(
    "--analysis",
    "-a",
    default="heatmap",
    type=str,
    help="""There are two types of analysis:
        heatmap - Generates a heatmap plot that shows the differences between the semantic similarity profiles using the
        score column for this purpose. Defaults to "heatmap".
        percentage_diff - Calculates the score column percentage difference between the semantic similarity profiles""",
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    default="diff.semsim.tsv",
    help="Output path for the difference tsv. Defaults to percentage_diff.semsim.tsv",
)
def semsim_comparison(
    semsim_left: Path,
    semsim_right: Path,
    score_column: str,
    output: Path = "percentage_diff.semsim.tsv",
    analysis: str = "heatmap",
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
    analysis_type = ["heatmap", "percentage_diff"]
    if analysis not in analysis_type:
        raise ValueError(
            f"""Invalid analysis type: {analysis}.
    These are the valid analysis types:{','.join(analysis_type)}"""
        )

    if analysis == "heatmap":
        return semsim_heatmap_plot(semsim_left, semsim_right, score_column)
    if analysis == "percentage_diff":
        percentage_diff(semsim_left, semsim_right, score_column, output)
