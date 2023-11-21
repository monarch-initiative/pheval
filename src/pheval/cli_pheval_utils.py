"""PhEval utils Command Line Interface"""

from pathlib import Path
from typing import List

import click

from pheval.analyse.analysis import (
    TrackInputOutputDirectories,
    benchmark_directory,
    benchmark_run_comparisons,
)
from pheval.analyse.run_data_parser import parse_run_data_text_file
from pheval.prepare.create_noisy_phenopackets import scramble_phenopackets
from pheval.prepare.create_spiked_vcf import spike_vcfs
from pheval.prepare.custom_exceptions import InputError, MutuallyExclusiveOptionError
from pheval.prepare.update_phenopacket import update_phenopackets
from pheval.utils.semsim_utils import percentage_diff, semsim_heatmap_plot
from pheval.utils.utils import semsim_convert, semsim_scramble


@click.command("semsim-scramble")
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the semantic similarity profile to be scrambled.",
    type=Path,
)
@click.option(
    "--output",
    "-o",
    metavar="FILE",
    required=True,
    help="Path where the scrambled semsim file will be written.",
    type=Path,
)
@click.option(
    "--score-column",
    "-c",
    required=True,
    multiple=True,
    type=click.Choice(
        ["jaccard_similarity", "dice_similarity", "phenodigm_score"], case_sensitive=False
    ),
    help="Score column that will be scrambled",
)
@click.option(
    "--scramble-factor",
    "-s",
    metavar=float,
    default=0.5,
    show_default=True,
    type=float,
    help="""Scramble Magnitude (noise)
    that will be applied to semantic similarity score column (e.g. jaccard similarity).""",
)
def semsim_scramble_command(
    input: Path, output: Path, score_column: List[str], scramble_factor: float
):
    """Scrambles semsim profile multiplying score value by scramble factor
    Args:
        input (Path): Path file that points out to the semsim profile
        output (Path): Path file that points out to the output file
        score_column (List[str]): Score column(s) that will be scrambled
        scramble_factor (float): Scramble Magnitude
    """
    semsim_scramble(input, output, score_column, scramble_factor)


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
    "-c",
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


@click.command("update-phenopackets")
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
    help="Path to phenopacket directory for updating.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--output-dir",
    "-o",
    metavar="PATH",
    required=True,
    help="Path to write phenopacket.",
    type=Path,
)
@click.option(
    "--gene-identifier",
    "-g",
    required=False,
    default="ensembl_id",
    show_default=True,
    help="Gene identifier to add to phenopacket",
    type=click.Choice(["ensembl_id", "entrez_id", "hgnc_id"]),
)
def update_phenopackets_command(
    phenopacket_path: Path, phenopacket_dir: Path, output_dir: Path, gene_identifier: str
):
    """Update gene symbols and identifiers for phenopackets."""
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    update_phenopackets(gene_identifier, phenopacket_path, phenopacket_dir, output_dir)


@click.command("create-spiked-vcfs")
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
    help="Path to phenopacket directory for updating.",
    type=Path,
    cls=MutuallyExclusiveOptionError,
    mutually_exclusive=["phenopacket_path"],
)
@click.option(
    "--template-vcf-path",
    "-t",
    cls=MutuallyExclusiveOptionError,
    metavar="PATH",
    required=False,
    help="Template VCF file",
    mutually_exclusive=["vcf_dir"],
    type=Path,
)
@click.option(
    "--vcf-dir",
    "-v",
    cls=MutuallyExclusiveOptionError,
    metavar="PATH",
    help="Directory containing template VCF files",
    mutually_exclusive=["template_vcf"],
    type=Path,
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="vcf",
    type=Path,
)
def create_spiked_vcfs_command(
    phenopacket_path: Path,
    phenopacket_dir: Path,
    output_dir: Path,
    template_vcf_path: Path = None,
    vcf_dir: Path = None,
):
    """Spikes variants into a template VCF file for a directory of phenopackets."""
    if phenopacket_path is None and phenopacket_dir is None:
        raise InputError("Either a phenopacket or phenopacket directory must be specified")
    spike_vcfs(output_dir, phenopacket_path, phenopacket_dir, template_vcf_path, vcf_dir)


@click.command("semsim-convert")
@click.option(
    "--input",
    "-i",
    required=True,
    metavar="FILE",
    help="Path to the semsim file.",
    type=Path,
)
@click.option(
    "--output",
    "-o",
    required=True,
    metavar="FILE",
    help="Path where converted semsim will be written.",
    type=Path,
)
@click.option(
    "--subject-prefix",
    "-s",
    required=True,
    metavar="FILE",
    help="Subject Prefix that will be mapped to the database",
    type=str,
)
@click.option(
    "--object-prefix",
    "-b",
    required=True,
    metavar="FILE",
    help="Object Prefix that will be mapped to the database.",
    type=str,
)
@click.option(
    "--output-format",
    "-O",
    required=True,
    metavar=str,
    help="Output file format. Available formats: (exomiserdb)",
    type=click.Choice(["exomiserdb"], case_sensitive=False),
)
def semsim_convert_command(
    input: Path, output: Path, subject_prefix: str, object_prefix: str, output_format: str
):
    """convert semsim profile to an exomiser database file"""
    semsim_convert(input, output, subject_prefix, object_prefix, output_format)


@click.command()
@click.option(
    "--directory",
    "-d",
    required=True,
    metavar="PATH",
    help="General results directory to be benchmarked, assumes contains subdirectories of pheval_gene_results/,"
    "pheval_variant_results/ or pheval_disease_results/. ",
    type=Path,
)
@click.option(
    "--phenopacket-dir",
    "-p",
    required=True,
    metavar="PATH",
    help="Full path to directory containing input phenopackets.",
    type=Path,
)
@click.option(
    "--output-prefix",
    "-o",
    metavar="<str>",
    required=True,
    help=" Output file prefix. ",
)
@click.option(
    "--score-order",
    "-so",
    required=True,
    help="Ordering of results for ranking.",
    type=click.Choice(["ascending", "descending"]),
    default="descending",
    show_default=True,
)
@click.option(
    "--threshold",
    "-t",
    metavar="<float>",
    default=float(0.0),
    required=False,
    help="Score threshold.",
    type=float,
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for variant prioritisation",
)
@click.option(
    "--disease-analysis/--no-disease-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for disease prioritisation",
)
@click.option(
    "--plot-type",
    "-y",
    default="bar_stacked",
    show_default=True,
    type=click.Choice(["bar_stacked", "bar_cumulative", "bar_non_cumulative"]),
    help="Bar chart type to output.",
)
def benchmark(
    directory: Path,
    phenopacket_dir: Path,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
):
    """Benchmark the gene/variant/disease prioritisation performance for a single run."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify at least one of gene/variant/disease analysis.")
    benchmark_directory(
        TrackInputOutputDirectories(results_dir=directory, phenopacket_dir=phenopacket_dir),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        disease_analysis,
        plot_type,
    )


@click.command()
@click.option(
    "--run-data",
    "-r",
    required=True,
    metavar="PATH",
    help="Path to .txt file containing testdata phenopacket directory "
    "and corresponding results directory separated by tab."
    "Each run contained to a new line with the input testdata listed first and on the same line separated by a tab"
    "the results directory.",
    type=Path,
)
@click.option(
    "--output-prefix",
    "-o",
    metavar="<str>",
    required=True,
    help=" Output file prefix. ",
)
@click.option(
    "--score-order",
    "-so",
    required=True,
    help="Ordering of results for ranking.",
    type=click.Choice(["ascending", "descending"]),
    default="descending",
    show_default=True,
)
@click.option(
    "--threshold",
    "-t",
    metavar="<float>",
    default=float(0.0),
    required=False,
    help="Score threshold.",
    type=float,
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for variant prioritisation",
)
@click.option(
    "--disease-analysis/--no-disease-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for disease prioritisation",
)
@click.option(
    "--plot-type",
    "-y",
    default="bar_stacked",
    show_default=True,
    type=click.Choice(["bar_stacked", "bar_cumulative", "bar_non_cumulative"]),
    help="Bar chart type to output.",
)
def benchmark_comparison(
    run_data: Path,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
):
    """Benchmark the gene/variant/disease prioritisation performance for two runs."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify at least one of gene/variant/disease analysis.")
    benchmark_run_comparisons(
        parse_run_data_text_file(run_data),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        disease_analysis,
        plot_type,
    )
