from pheval.analyse.benchmark_generator import (
    BenchmarkRunOutputGenerator,
    DiseaseBenchmarkRunOutputGenerator,
    GeneBenchmarkRunOutputGenerator,
    VariantBenchmarkRunOutputGenerator,
)
from pheval.analyse.generate_summary_outputs import (
    generate_benchmark_comparison_output,
)
from pheval.analyse.parse_corpus import CorpusParser
from pheval.analyse.rank_stats import RankStatsWriter
from pheval.analyse.run_data_parser import Config


def _run_benchmark_comparison(
        run_config: Config,
        score_order: str,
        threshold: float,
        benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Run a benchmark on several result directories.

    Args:
        run_config (List[TrackInputOutputDirectories]): List of input and output directories
            for tracking results across multiple directories.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        threshold (float): The threshold for benchmark evaluation.
        benchmark_generator (BenchmarkRunOutputGenerator): Generator for benchmark run output.
    """
    stats_writer = RankStatsWriter(run_config.benchmark_name, benchmark_generator.stats_comparison_file)
    unique_test_corpora_directories = set([result.phenopacket_dir for result in run_config.runs])
    [
        CorpusParser(run_config.benchmark_name, test_corpora_directory).parse_corpus(benchmark_generator)
        for test_corpora_directory in unique_test_corpora_directories
    ]
    benchmarking_results = []
    for run in run_config.runs:
        benchmark_result = benchmark_generator.generate_benchmark_run_results(
            run_config.benchmark_name, run, score_order, threshold
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
            run_config.benchmark_name,
            benchmarking_results,
            run_identifiers,
            benchmark_generator,
            f"{unique_test_corpora_directory.parents[0].name}_"
            f"{benchmark_generator.prioritisation_type_string}",
        )
        for unique_test_corpora_directory in unique_test_corpora_directories
    ]


def benchmark_run_comparisons(
        run_config: Config,
        score_order: str,
        threshold: float,
) -> None:
    """
    Benchmark prioritisation performance for several runs.

    Args:
        run_config (Config): Run configurations.
        score_order (str): The order in which scores are arranged, this can be either ascending or descending.
        threshold (float): The threshold for benchmark evaluation.
    """
    gene_analysis_runs = Config(benchmark_name=run_config.benchmark_name,
                                runs=[run for run in run_config.runs if run.gene_analysis],
                                plot_customisation=run_config.plot_customisation)
    variant_analysis_runs = Config(benchmark_name=run_config.benchmark_name,
                                   runs=[run for run in run_config.runs if run.variant_analysis],
                                   plot_customisation=run_config.plot_customisation)
    disease_analysis_runs = Config(benchmark_name=run_config.benchmark_name,
                                   runs=[run for run in run_config.runs if run.disease_analysis],
                                   plot_customisation=run_config.plot_customisation)
    if gene_analysis_runs.runs:
        _run_benchmark_comparison(
            run_config=gene_analysis_runs,
            score_order=score_order,
            threshold=threshold,
            benchmark_generator=GeneBenchmarkRunOutputGenerator(
                plot_customisation=gene_analysis_runs.plot_customisation.gene_plots),
        )
    if variant_analysis_runs.runs:
        _run_benchmark_comparison(
            run_config=variant_analysis_runs,
            score_order=score_order,
            threshold=threshold,
            benchmark_generator=VariantBenchmarkRunOutputGenerator(
                plot_customisation=variant_analysis_runs.plot_customisation.variant_plots),
        )
    if disease_analysis_runs.runs:
        _run_benchmark_comparison(
            run_config=disease_analysis_runs,
            score_order=score_order,
            threshold=threshold,
            benchmark_generator=DiseaseBenchmarkRunOutputGenerator(
                plot_customisation=disease_analysis_runs.plot_customisation.disease_plots),
        )
