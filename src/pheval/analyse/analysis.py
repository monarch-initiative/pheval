from collections import defaultdict
from pathlib import Path

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
from pheval.analyse.rank_stats_writer import RankStatsWriter
from pheval.analyse.run_data_parser import TrackInputOutputDirectories


def run_benchmark(
    results_dir_and_input: TrackInputOutputDirectories,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Run a benchmark on a result directory."""
    stats_writer = RankStatsWriter(
        Path(output_prefix + benchmark_generator.rank_comparison_file_suffix)
    )
    rank_comparison = defaultdict(dict)
    benchmark_result = benchmark_generator.benchmark_function(
        results_dir_and_input, score_order, threshold, rank_comparison, stats_writer
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
    """Benchmark prioritisation performance for a result directory."""
    if gene_analysis:
        run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if variant_analysis:
        run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if disease_analysis:
        run_benchmark(
            results_dir_and_input=results_dir_and_input,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )


def run_benchmark_comparison(
    results_directories: [TrackInputOutputDirectories],
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Run a benchmark on several result directories."""
    stats_writer = RankStatsWriter(
        Path(output_prefix + benchmark_generator.rank_comparison_file_suffix)
    )
    benchmarking_results = []
    for results_dir_and_input in results_directories:
        rank_comparison = defaultdict(dict)
        benchmark_result = benchmark_generator.benchmark_function(
            results_dir_and_input, score_order, threshold, rank_comparison, stats_writer
        )
        benchmarking_results.append(benchmark_result)
    generate_benchmark_comparison_output(benchmarking_results, plot_type, benchmark_generator)
    stats_writer.close()


def benchmark_run_comparisons(
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
    if gene_analysis:
        run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if variant_analysis:
        run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if disease_analysis:
        run_benchmark_comparison(
            results_directories=results_directories,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )
