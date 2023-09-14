from collections import defaultdict
from pathlib import Path

import click

from pheval.analyse.disease_prioritisation_analysis import assess_phenopacket_disease_prioritisation
from pheval.analyse.gene_prioritisation_analysis import assess_phenopacket_gene_prioritisation
from pheval.analyse.generate_plots import AnalysisResults, TrackRunPrioritisation
from pheval.analyse.generate_summary_outputs import (
    RankStatsWriter,
    generate_benchmark_comparison_disease_output,
    generate_benchmark_comparison_gene_output,
    generate_benchmark_comparison_variant_output,
    generate_benchmark_disease_output,
    generate_benchmark_gene_output,
    generate_benchmark_variant_output,
)
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.run_data_parser import TrackInputOutputDirectories, _parse_run_data_text_file
from pheval.analyse.variant_prioritisation_analysis import assess_phenopacket_variant_prioritisation
from pheval.prepare.custom_exceptions import InputError
from pheval.utils.file_utils import files_with_suffix


def _assess_prioritisation_for_results_directory(
    results_directory_and_input: TrackInputOutputDirectories,
    score_order: str,
    threshold: float,
    gene_rank_comparison: defaultdict,
    variant_rank_comparison: defaultdict,
    disease_rank_comparison: defaultdict,
    gene_stats_writer: RankStatsWriter,
    variants_stats_writer: RankStatsWriter,
    disease_stats_writer: RankStatsWriter,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
) -> TrackRunPrioritisation:
    """Assess prioritisation for a single results directory."""
    gene_rank_stats, variant_rank_stats, disease_rank_stats = RankStats(), RankStats(), RankStats()
    if gene_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_gene_results/"), ".tsv"
        ):
            assess_phenopacket_gene_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                gene_rank_stats,
                gene_rank_comparison,
            )
        gene_stats_writer.write_row(results_directory_and_input.results_dir, gene_rank_stats)
    if variant_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_variant_results/"),
            ".tsv",
        ):
            assess_phenopacket_variant_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                variant_rank_stats,
                variant_rank_comparison,
            )
        variants_stats_writer.write_row(results_directory_and_input.results_dir, variant_rank_stats)
    if disease_analysis:
        for standardised_result in files_with_suffix(
            results_directory_and_input.results_dir.joinpath("pheval_disease_results/"),
            ".tsv",
        ):
            assess_phenopacket_disease_prioritisation(
                standardised_result,
                score_order,
                results_directory_and_input,
                threshold,
                disease_rank_stats,
                disease_rank_comparison,
            )
        disease_stats_writer.write_row(results_directory_and_input.results_dir, disease_rank_stats)
    return TrackRunPrioritisation(
        gene_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=gene_rank_comparison,
            rank_stats=gene_rank_stats,
        ),
        variant_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=variant_rank_comparison,
            rank_stats=variant_rank_stats,
        ),
        disease_prioritisation=AnalysisResults(
            results_dir=results_directory_and_input.results_dir,
            ranks=disease_rank_comparison,
            rank_stats=disease_rank_stats,
        ),
    )


def benchmark_directory(
    results_dir_and_input: TrackInputOutputDirectories,
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark prioritisation performance for a single directory."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    disease_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-disease_summary.tsv")) if disease_analysis else None
    )
    gene_rank_comparison, variant_rank_comparison, disease_rank_comparison = (
        defaultdict(dict),
        defaultdict(dict),
        defaultdict(dict),
    )
    prioritisation_data = _assess_prioritisation_for_results_directory(
        results_dir_and_input,
        score_order,
        threshold,
        gene_rank_comparison,
        variant_rank_comparison,
        disease_rank_comparison,
        gene_stats_writer,
        variants_stats_writer,
        disease_stats_writer,
        gene_analysis,
        variant_analysis,
        disease_analysis,
    )
    generate_benchmark_gene_output(prioritisation_data, plot_type) if gene_analysis else None
    generate_benchmark_variant_output(prioritisation_data, plot_type) if variant_analysis else None
    generate_benchmark_disease_output(prioritisation_data, plot_type) if disease_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None
    disease_stats_writer.close() if disease_analysis else None


def benchmark_runs(
    results_directories: [TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
) -> None:
    """Benchmark several result directories."""
    gene_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-gene_summary.tsv")) if gene_analysis else None
    )
    variants_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-variant_summary.tsv")) if variant_analysis else None
    )
    disease_stats_writer = (
        RankStatsWriter(Path(output_prefix + "-disease_summary.tsv")) if disease_analysis else None
    )
    prioritisation_stats_for_runs = []
    for results_dir_and_input in results_directories:
        gene_rank_comparison, variant_rank_comparison, disease_rank_comparison = (
            defaultdict(dict),
            defaultdict(dict),
            defaultdict(dict),
        )
        prioritisation_stats = _assess_prioritisation_for_results_directory(
            results_dir_and_input,
            score_order,
            threshold,
            gene_rank_comparison,
            variant_rank_comparison,
            disease_rank_comparison,
            gene_stats_writer,
            variants_stats_writer,
            disease_stats_writer,
            gene_analysis,
            variant_analysis,
            disease_analysis,
        )
        prioritisation_stats_for_runs.append(prioritisation_stats)
    generate_benchmark_comparison_gene_output(
        prioritisation_stats_for_runs, plot_type
    ) if gene_analysis else None
    generate_benchmark_comparison_variant_output(
        prioritisation_stats_for_runs, plot_type
    ) if variant_analysis else None
    generate_benchmark_comparison_disease_output(
        prioritisation_stats_for_runs, plot_type
    ) if disease_analysis else None
    gene_stats_writer.close() if gene_analysis else None
    variants_stats_writer.close() if variant_analysis else None
    disease_stats_writer.close() if disease_analysis else None


@click.command()
@click.option(
    "--directory",
    "-d",
    required=True,
    metavar="PATH",
    help="General results directory to be benchmarked, assumes contains subdirectories of pheval_gene_results/"
    "pheval_variant_results and the tool specific results directory. ",
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
    """Benchmark the gene/variant prioritisation performance for a single run."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify gene analysis and/or variant and/or disease analysis.")
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
    help="Path to .txt file containing testdata directory and corresponding results directory separated by tab."
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
    """Benchmark the gene/variant prioritisation performance for two runs."""
    if not gene_analysis and not variant_analysis and not disease_analysis:
        raise InputError("Need to specify gene analysis and/or variant and/or disease analysis.")
    benchmark_runs(
        _parse_run_data_text_file(run_data),
        score_order,
        output_prefix,
        threshold,
        gene_analysis,
        variant_analysis,
        disease_analysis,
        plot_type,
    )
