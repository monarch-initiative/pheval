from dataclasses import dataclass
from pathlib import Path

from pheval.analyse.rank_stats import RankStats


@dataclass
class BenchmarkRunResults:
    """Benchmarking results for a run."""

    ranks: dict
    rank_stats: RankStats
    results_dir: Path = None
    benchmark_name: str = None
