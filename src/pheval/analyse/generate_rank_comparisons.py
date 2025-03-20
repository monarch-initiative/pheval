from itertools import combinations
from typing import List

import polars as pl
from duckdb.duckdb import DuckDBPyConnection

from pheval.analyse.benchmark_db_manager import write_table
from pheval.analyse.benchmark_output_type import BenchmarkOutputType
from pheval.utils.logger import get_logger


def calculate_rank_changes(
    conn: DuckDBPyConnection,
    run_identifiers: List[str],
    true_positive_cases: pl.DataFrame,
    benchmark_type: BenchmarkOutputType,
) -> None:
    """
    Calculate rank changes between runs.
    Args:
        conn (DuckDBPyConnection): DuckDB connection.
        run_identifiers (List[str]): List of run identifiers.
        true_positive_cases (pl.LazyFrame): All true positive cases for a benchmark.
        benchmark_type (BenchmarkOutputType): Type of benchmark output.
    """
    logger = get_logger()
    pairwise_comparisons = list(combinations(run_identifiers, 2))
    for col1, col2 in pairwise_comparisons:
        logger.info(f"Comparing rank changes: {col1} vs. {col2}")
        rank_change_lf = true_positive_cases.with_columns(
            [
                pl.when((pl.col(col1) == 0) & (pl.col(col2) != 0))
                .then(pl.lit("GAINED"))
                .when((pl.col(col1) != 0) & (pl.col(col2) == 0))
                .then(pl.lit("LOST"))
                .otherwise((pl.col(col1) - pl.col(col2)).cast(pl.Int64))
                .alias("rank_change")
            ]
        ).select(["result_file", *benchmark_type.columns, col1, col2, "rank_change"])
        write_table(
            conn,
            rank_change_lf,
            f"{col1}_vs_{col2}_{benchmark_type.prioritisation_type_string}_rank_changes",
        )
