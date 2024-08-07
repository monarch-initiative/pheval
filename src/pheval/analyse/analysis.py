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
from pheval.analyse.parse_corpus import CorpusParser
from pheval.analyse.rank_stats import RankStatsWriter
from pheval.analyse.run_data_parser import Config, RunConfig


def _run_benchmark(
    run_config: RunConfig,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Run a benchmark on a result directory.

    Args:
        run_config (RunConfig): Run configuration.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
        benchmark_generator (BenchmarkRunOutputGenerator): Generator for benchmark run output.
    """
    CorpusParser(run_config.phenopacket_dir).parse_corpus(benchmark_generator)
    stats_writer = RankStatsWriter(
        str(output_prefix + benchmark_generator.stats_comparison_file_suffix)
    )
    benchmark_result = benchmark_generator.generate_benchmark_run_results(
        run_config, score_order, threshold
    )
    stats_writer.add_statistics_entry(
        run_config.run_identifier,
        benchmark_result.rank_stats,
        benchmark_result.binary_classification_stats,
    )
    generate_benchmark_output(benchmark_result, plot_type, benchmark_generator)


def benchmark_directory(
    run_config: RunConfig,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
) -> None:
    """
    Benchmark prioritisation performance for a single run.

    Args:
        run_config (RunConfig): Run configuration.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
    """
    if run_config.gene_analysis:
        _run_benchmark(
            run_config=run_config,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if run_config.variant_analysis:
        _run_benchmark(
            run_config=run_config,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if run_config.disease_analysis:
        _run_benchmark(
            run_config=run_config,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )


def _run_benchmark_comparison(
    run_config: Config,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Run a benchmark on several result directories.

    Args:
        run_config (List[TrackInputOutputDirectories]): List of input and output directories
            for tracking results across multiple directories.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
        benchmark_generator (BenchmarkRunOutputGenerator): Generator for benchmark run output.
    """
    stats_writer = RankStatsWriter(
        str(output_prefix + benchmark_generator.stats_comparison_file_suffix)
    )
    unique_test_corpora_directories = set([result.phenopacket_dir for result in run_config.runs])
    [
        CorpusParser(test_corpora_directory).parse_corpus(benchmark_generator)
        for test_corpora_directory in unique_test_corpora_directories
    ]
    benchmarking_results = []
    for run in run_config.runs:
        benchmark_result = benchmark_generator.generate_benchmark_run_results(
            run, score_order, threshold
        )
        stats_writer.add_statistics_entry(
            run.run_identifier,
            benchmark_result.rank_stats,
            benchmark_result.binary_classification_stats,
        )
        benchmarking_results.append(benchmark_result)
    run_identifiers = [run.run_identifier for run in run_config.runs]
    [
        generate_benchmark_comparison_output(
            benchmarking_results,
            run_identifiers,
            plot_type,
            benchmark_generator,
            f"{unique_test_corpora_directory.parents[0].name}_"
            f"{benchmark_generator.prioritisation_type_string}",
        )
        for unique_test_corpora_directory in unique_test_corpora_directories
    ]


def benchmark_run_comparisons(
    run_config: Config,
    score_order: str,
    output_prefix: str,
    threshold: float,
    plot_type: str,
) -> None:
    """
    Benchmark prioritisation performance for several runs.

    Args:
        run_config (Config): Run configurations.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        output_prefix (str): Prefix for the benchmark output file names.
        threshold (float): The threshold for benchmark evaluation.
        plot_type (str): Type of plot for benchmark visualisation.
    """
    gene_analysis_runs = Config(runs=[run for run in run_config.runs if run.gene_analysis])
    variant_analysis_runs = Config(runs=[run for run in run_config.runs if run.variant_analysis])
    disease_analysis_runs = Config(runs=[run for run in run_config.runs if run.disease_analysis])
    if gene_analysis_runs:
        _run_benchmark_comparison(
            run_config=gene_analysis_runs,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(),
        )
    if variant_analysis_runs:
        _run_benchmark_comparison(
            run_config=variant_analysis_runs,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(),
        )
    if disease_analysis_runs:
        _run_benchmark_comparison(
            run_config=disease_analysis_runs,
            score_order=score_order,
            output_prefix=output_prefix,
            threshold=threshold,
            plot_type=plot_type,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(),
        )
