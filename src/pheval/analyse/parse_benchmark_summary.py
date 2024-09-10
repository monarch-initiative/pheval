from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats


@dataclass
class BenchmarkSummaryResults:
    gene_results: List[BenchmarkRunResults]
    disease_results: List[BenchmarkRunResults]
    variant_results: List[BenchmarkRunResults]


def parse_benchmark_results(benchmark_summary_table: pd.DataFrame) -> List[BenchmarkRunResults]:
    """
    Parse benchmark results from a DataFrame.

    Args:
        benchmark_summary_table (pd.DataFrame): DataFrame containing benchmark results.

    Returns:
        List[BenchmarkRunResults]: A list of BenchmarkRunResults objects parsed from the DataFrame.
    """
    results = []
    for _, row in benchmark_summary_table.iterrows():
        benchmarking_result = BenchmarkRunResults(
            rank_stats=RankStats(
                top=row["top"],
                top3=row["top3"],
                top5=row["top5"],
                top10=row["top10"],
                found=row["found"],
                total=row["total"],
                mrr=row["mean_reciprocal_rank"],
            ),
            benchmark_name=row["results_directory_path"],
            binary_classification_stats=BinaryClassificationStats(),
        )
        results.append(benchmarking_result)
    return results


def parse_benchmark_db(benchmarking_db: Path) -> BenchmarkSummaryResults:
    """
    Read the summary benchmark TSV output generated from the benchmark-comparison command.

    Args:
        benchmarking_db (Path): Path to the benchmark db.

    Returns:
        BenchmarkSummaryResults: A dataclass containing all benchmarking results contained in the db.
    """
    db_connector = BenchmarkDBManager(benchmarking_db)
    gene_benchmarking_results, disease_benchmarking_results, variant_benchmarking_results = (
        None,
        None,
        None,
    )
    if db_connector.check_table_exists("gene_summary"):
        gene_benchmarking_results = parse_benchmark_results(
            db_connector.conn.execute("SELECT * FROM gene_summary").fetchdf()
        )
    if db_connector.check_table_exists("disease_summary"):
        disease_benchmarking_results = parse_benchmark_results(
            db_connector.conn.execute("SELECT * FROM disease_summary").fetchdf()
        )
    if db_connector.check_table_exists("variant_summary"):
        variant_benchmarking_results = parse_benchmark_results(
            db_connector.conn.execute("SELECT * FROM variant_summary").fetchdf()
        )
    return BenchmarkSummaryResults(
        gene_results=gene_benchmarking_results,
        disease_results=disease_benchmarking_results,
        variant_results=variant_benchmarking_results,
    )
