from pathlib import Path

import pandas as pd

from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.rank_stats import RankStats


def read_benchmark_tsv_result_summary(benchmarking_tsv: Path) -> pd.DataFrame:
    """Read the summary benchmark tsv output from the benchmark-comparison command."""
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


def parse_benchmark_result_summary(benchmarking_df: pd.DataFrame) -> [BenchmarkRunResults]:
    """Parse the summary benchmark dataframe into BenchmarkRunResults."""
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
        )
        benchmarking_results.append(benchmarking_result)
    return benchmarking_results
