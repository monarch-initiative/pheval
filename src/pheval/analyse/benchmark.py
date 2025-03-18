import time
from pathlib import Path
from typing import List, Tuple

import duckdb
import polars as pl

from pheval.analyse.benchmark_db_manager import write_table
from pheval.analyse.benchmark_output_type import BenchmarkOutputType, BenchmarkOutputTypeEnum
from pheval.analyse.binary_classification_curves import compute_curves
from pheval.analyse.binary_classification_stats import compute_confusion_matrix
from pheval.analyse.generate_plots import generate_plots
from pheval.analyse.generate_rank_comparisons import calculate_rank_changes
from pheval.analyse.rank_stats import compute_rank_stats
from pheval.analyse.run_data_parser import Config, RunConfig, parse_run_config
from pheval.utils.logger import get_logger


def scan_directory(run: RunConfig, benchmark_type: BenchmarkOutputType) -> pl.LazyFrame:
    """
    Scan a results directory containing pheval parquet standardised results and return a LazyFrame object.
    Args:
        run (RunConfig): RunConfig object.
        benchmark_type (BenchmarkOutputTypeEnum): Benchmark output type.
    Returns:
        pl.LazyFrame: LazyFrame object containing all the results in the directory..
    """
    logger = get_logger()
    logger.info(f"Analysing results in {run.results_dir.joinpath(benchmark_type.result_directory)}")
    return (
        pl.scan_parquet(
            run.results_dir.joinpath(benchmark_type.result_directory),
            include_file_paths="file_path",
        ).with_columns(
            pl.col("rank").cast(pl.Int64),
            pl.col("file_path").str.extract(r"([^/\\]+)$").alias("result_file"),
            pl.col("true_positive").fill_null(False),
        )
    ).filter(
        (
            pl.col("score") >= run.threshold
            if run.score_order.lower() == "descending"
            else pl.col("score") <= run.threshold
        )
        if run.threshold is not None
        else True
    )


def process_stats(
    runs: List[RunConfig], benchmark_type: BenchmarkOutputType
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Processes stats outputs for specified runs to compare.
    Args:
        runs (List[RunConfig]): List of runs to benchmark.
        benchmark_type (BenchmarkOutputTypeEnum): Benchmark output type.
    Returns:
        Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]: The stats for all runs.
    """
    stats, curve_results, true_positive_cases = [], [], []
    for run in runs:
        result_scan = scan_directory(run, benchmark_type)
        stats.append(
            compute_rank_stats(run.run_identifier, result_scan).join(
                compute_confusion_matrix(run.run_identifier, result_scan), on="run_identifier"
            )
        )
        curve_results.append(compute_curves(run.run_identifier, result_scan))
        true_positive_cases.append(
            result_scan.filter(pl.col("true_positive")).select(
                ["result_file", *benchmark_type.columns, pl.col("rank").alias(run.run_identifier)]
            )
        )
    return (
        pl.concat(stats, how="vertical").collect(),
        pl.concat(curve_results, how="vertical").collect(),
        pl.concat(true_positive_cases, how="align_inner").collect(),
    )


def benchmark(config: Config, benchmark_type: BenchmarkOutputType) -> None:
    """
    Benchmark results for specified runs for a specified prioritisation type for comparison.
    Args:
        config (Config): Configuration for benchmarking.
        benchmark_type (BenchmarkOutputType): Benchmark output type.
    """
    conn = duckdb.connect(f"{config.benchmark_name}.duckdb")
    stats, curve_results, true_positive_cases = process_stats(config.runs, benchmark_type)
    write_table(
        conn, stats, f"{config.benchmark_name}_{benchmark_type.prioritisation_type_string}_summary"
    )
    write_table(
        conn,
        curve_results,
        f"{config.benchmark_name}_{benchmark_type.prioritisation_type_string}_binary_classification_curves",
    )
    calculate_rank_changes(
        conn, [run.run_identifier for run in config.runs], true_positive_cases, benchmark_type
    )
    generate_plots(
        config.benchmark_name, stats, curve_results, benchmark_type, config.plot_customisation
    )
    conn.close()


def benchmark_runs(benchmark_config_file: Path) -> None:
    """
    Benchmark results for specified runs for comparison.
    Args:
        benchmark_config_file (Path): Path to benchmark config file.
    """
    logger = get_logger()
    start_time = time.perf_counter()
    logger.info("Initiated benchmarking process.")
    config = parse_run_config(benchmark_config_file)
    gene_analysis_runs = [run for run in config.runs if run.gene_analysis]
    variant_analysis_runs = [run for run in config.runs if run.variant_analysis]
    disease_analysis_runs = [run for run in config.runs if run.disease_analysis]
    if gene_analysis_runs:
        logger.info("Initiating benchmarking for gene results.")
        benchmark(
            Config(
                benchmark_name=config.benchmark_name,
                runs=gene_analysis_runs,
                plot_customisation=config.plot_customisation,
            ),
            BenchmarkOutputTypeEnum.GENE.value,
        )
        logger.info("Finished benchmarking for gene results.")
    if variant_analysis_runs:
        logger.info("Initiating benchmarking for variant results")
        benchmark(
            Config(
                benchmark_name=config.benchmark_name,
                runs=variant_analysis_runs,
                plot_customisation=config.plot_customisation,
            ),
            BenchmarkOutputTypeEnum.VARIANT.value,
        )
        logger.info("Finished benchmarking for variant results.")
    if disease_analysis_runs:
        logger.info("Initiating benchmarking for disease results")
        benchmark(
            Config(
                benchmark_name=config.benchmark_name,
                runs=disease_analysis_runs,
                plot_customisation=config.plot_customisation,
            ),
            BenchmarkOutputTypeEnum.DISEASE.value,
        )
        logger.info("Finished benchmarking for disease results.")
    logger.info(
        f"Finished benchmarking! Total time: {time.perf_counter() - start_time:.2f} seconds."
    )
