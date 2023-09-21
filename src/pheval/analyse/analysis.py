from collections import defaultdict
from pathlib import Path

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
from pheval.analyse.run_data_parser import TrackInputOutputDirectories
from pheval.analyse.variant_prioritisation_analysis import assess_phenopacket_variant_prioritisation
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
