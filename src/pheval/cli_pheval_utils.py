"""PhEval utils Command Line Interface"""

from pathlib import Path

import click

from pheval.prepare.create_noisy_phenopackets import create_scrambled_phenopackets
from pheval.utils.semsim_utils import diff, semsim_heatmap_plot


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
    default='heatmap',
    type=str,
    help="Score column that will be used in comparison (e.g jaccard_similarity)",
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    default='semsim_analysis.tsv',
    help="Score column that will be used in comparison (e.g jaccard_similarity)",
)
def semsim_comparison(semsim_left: click.Path,
                      semsim_right: click.Path,
                      score_column: str,
                      analysis: str,
                      output: click.Path = 'semsim_analysis.tsv'):
    """Compares two semantic similarity profiles

    Args:
        semsim-left (click.Path): File path of the first semantic similarity profile
        semsim-right (click.Path): File path of the second semantic similarity profile
        score_column (str): Score column that will be computed (e.g. jaccard_similarity)
        analysis (str): [description]
        output (click.Path): Output path for the difference tsv
    """

    match analysis:
        case 'heatmap':
            semsim_heatmap_plot(semsim_left, semsim_right, score_column)
        case 'diff':
            diff(semsim_left, semsim_right, score_column, output)
