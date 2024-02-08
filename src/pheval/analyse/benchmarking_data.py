from dataclasses import dataclass
from pathlib import Path

from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.rank_stats import RankStats


@dataclass
class BenchmarkRunResults:
    """
    Benchmarking results for a run.

    Attributes:
        ranks (dict): Dictionary containing recorded ranks for samples.
        rank_stats (RankStats): Statistics related to benchmark.
        results_dir (Path, optional): Path to the result directory. Defaults to None.
        benchmark_name (str, optional): Name of the benchmark run. Defaults to None.
    """

    ranks: dict
    rank_stats: RankStats
    binary_classification_stats: BinaryClassificationStats
    results_dir: Path = None
    benchmark_name: str = None
