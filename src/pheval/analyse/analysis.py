from collections import defaultdict
from pathlib import Path
from typing import List

from pheval.analyse.benchmark_generator import (
    BenchmarkRunOutputGenerator,
    DiseaseBenchmarkRunOutputGenerator,
    GeneBenchmarkRunOutputGenerator,
    VariantBenchmarkRunOutputGenerator,
)
from pheval.analyse.generate_summary_outputs import (
    generate_benchmark_comparison_output,
    generate_benchmark_output,
)
from pheval.analyse.rank_stats import RankStatsWriter
from pheval.analyse.run_data_parser import TrackInputOutputDirectories


def _run_benchmark(
    results_dir_and_input: TrackInputOutputDirectories,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Run a benchmark on a result directory.

    Args:
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories for tracking results.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
        benchmark_generator (BenchmarkRunOutputGenerator): Generator for benchmark run output.
    """
    stats_writer = RankStatsWriter(
        Path(output_prefix + benchmark_generator.stats_comparison_file_suffix)
    )
    rank_comparison = defaultdict(dict)
    benchmark_result = benchmark_generator.generate_benchmark_run_results(
        results_dir_and_input, score_order, threshold, rank_comparison
    )
    stats_writer.write_row(
        results_dir_and_input.results_dir,
        benchmark_result.rank_stats,
        benchmark_result.binary_classification_stats,
    )
    generate_benchmark_output(benchmark_result, plot_type, benchmark_generator)
    stats_writer.close()


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
    """
    Benchmark prioritisation performance for a single run.

    Args:
        results_dir_and_input (TrackInputOutputDirectories): Input and output directories for tracking results.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        gene_analysis (bool): Boolean flag indicating whether to benchmark gene results.
        variant_analysis (bool): Boolean flag indicating whether to benchmark variant results.
        disease_analysis (bool): Boolean flag indicating whether to benchmark disease results.
        plot_type (str): Type of plot for benchmark visualisation.
    """
    if gene_analysis:
        _run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if variant_analysis:
        _run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if disease_analysis:
        _run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )


def _run_benchmark_comparison(
    results_directories: List[TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Run a benchmark on several result directories.

    Args:
        results_directories (List[TrackInputOutputDirectories]): List of input and output directories
            for tracking results across multiple directories.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
        benchmark_generator (BenchmarkRunOutputGenerator): Generator for benchmark run output.
    """
    stats_writer = RankStatsWriter(
        Path(output_prefix + benchmark_generator.stats_comparison_file_suffix)
    )
    benchmarking_results = []
    for results_dir_and_input in results_directories:
        rank_comparison = defaultdict(dict)
        benchmark_result = benchmark_generator.generate_benchmark_run_results(
            results_dir_and_input, score_order, threshold, rank_comparison
        )
        stats_writer.write_row(
            results_dir_and_input.results_dir,
            benchmark_result.rank_stats,
            benchmark_result.binary_classification_stats,
        )
        benchmarking_results.append(benchmark_result)
    generate_benchmark_comparison_output(benchmarking_results, plot_type, benchmark_generator)
    stats_writer.close()


def benchmark_run_comparisons(
    results_directories: List[TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
) -> None:
    """
    Benchmark prioritisation performance for several runs.

    Args:
        results_directories (List[TrackInputOutputDirectories]): Input and output directories for tracking results.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        gene_analysis (bool): Boolean flag indicating whether to benchmark gene results.
        variant_analysis (bool): Boolean flag indicating whether to benchmark variant results.
        disease_analysis (bool): Boolean flag indicating whether to benchmark disease results.
        plot_type (str): Type of plot for benchmark visualisation.
    """
    if gene_analysis:
        _run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if variant_analysis:
        _run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if disease_analysis:
        _run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )
