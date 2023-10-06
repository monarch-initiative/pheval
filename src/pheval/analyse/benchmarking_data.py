from dataclasses import dataclass
from pathlib import Path

from pheval.analyse.rank_stats import RankStats


@dataclass
class BenchmarkRunResults:
    """Analysis results for a run."""

    results_dir: Path
    ranks: dict
    rank_stats: RankStats
