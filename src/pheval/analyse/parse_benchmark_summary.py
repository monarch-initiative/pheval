from pathlib import Path
from typing import List

import pandas as pd

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats


def read_benchmark_tsv_result_summary(benchmarking_tsv: Path) -> pd.DataFrame:
    """
    Read the summary benchmark TSV output generated from the benchmark-comparison command.

    Args:
        benchmarking_tsv (Path): Path to the summary benchmark TSV output file.

    Returns:
        pd.DataFrame: A pandas DataFrame containing specific columns from the TSV file, including:
                      'results_directory_path', 'top', 'top3', 'top5', 'top10', 'found',
                      'total', 'mean_reciprocal_rank'.
    """
    return pd.read_csv(
        benchmarking_tsv,
        delimiter="\t",
        usecols=[
            "results_directory_path",
            "top",
            "top3",
            "top5",
            "top10",
            "found",
            "total",
            "mean_reciprocal_rank",
        ],
    )


def parse_benchmark_result_summary(benchmarking_df: pd.DataFrame) -> List[BenchmarkRunResults]:
    """
    Parse the summary benchmark DataFrame into a list of BenchmarkRunResults.

    Args:
        benchmarking_df (pd.DataFrame): Summary benchmark DataFrame containing columns such as
                                        'results_directory_path', 'top', 'top3', 'top5', 'top10',
                                        'found', 'total', 'mean_reciprocal_rank'.

    Returns:
        List[BenchmarkRunResults]: A list of BenchmarkRunResults instances generated from the DataFrame.
    """
    benchmarking_results = []
    for _, row in benchmarking_df.iterrows():
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
            ranks={},
            benchmark_name=row["results_directory_path"],
            binary_classification_stats=BinaryClassificationStats(),
        )
        benchmarking_results.append(benchmarking_result)
    return benchmarking_results
