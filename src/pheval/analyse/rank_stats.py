import csv
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean


@dataclass
class RankStats:
    """Class for keeping track of the rank stats."""

    top: int = 0
    top3: int = 0
    top5: int = 0
    top10: int = 0
    found: int = 0
    total: int = 0
    reciprocal_ranks: list = field(default_factory=list)
    mrr: float = None

    def add_rank(self, rank: int) -> None:
        """Add rank for phenopacket."""
        self.reciprocal_ranks.append(1 / rank)
        self.found += 1
        if rank == 1:
            self.top += 1
        if rank != "" and rank <= 3:
            self.top3 += 1
        if rank != "" and rank <= 5:
            self.top5 += 1
        if rank != "" and rank <= 10:
            self.top10 += 1

    def percentage_rank(self, value: int) -> float:
        """Return a percentage rank."""
        return 100 * value / self.total

    def percentage_top(self) -> float:
        """Return percentage of top matches."""
        return self.percentage_rank(self.top)

    def percentage_top3(self) -> float:
        """Return percentage of matches in the top3."""
        return self.percentage_rank(self.top3)

    def percentage_top5(self) -> float:
        """Return percentage of matches in the top5."""
        return self.percentage_rank(self.top5)

    def percentage_top10(self) -> float:
        """Return percentage of matches in the top10."""
        return self.percentage_rank(self.top10)

    def percentage_found(self) -> float:
        """Return percentage of matches found."""
        return self.percentage_rank(self.found)

    @staticmethod
    def percentage_difference(percentage_value_1: float, percentage_value_2: float) -> float:
        """Return percentage difference between two percentage values"""
        return percentage_value_1 - percentage_value_2

    def mean_reciprocal_rank(self) -> float:
        """Calculate the mean reciprocal rank."""
        if len(self.reciprocal_ranks) != self.total:
            missing_cases = self.total - self.found
            self.reciprocal_ranks.extend([0] * missing_cases)
            return mean(self.reciprocal_ranks)
        return mean(self.reciprocal_ranks)

    def return_mean_reciprocal_rank(self) -> float:
        """Return the mean reciprocal rank."""
        if self.mrr is not None:
            return self.mrr
        else:
            return self.mean_reciprocal_rank()


class RankStatsWriter:
    """Write the rank stats for each run."""

    def __init__(self, file: Path):
        self.file = open(file, "w")
        self.writer = csv.writer(self.file, delimiter="\t")
        self.writer.writerow(
            [
                "results_directory_path",
                "top",
                "top3",
                "top5",
                "top10",
                "found",
                "total",
                "mean_reciprocal_rank",
                "percentage_top",
                "percentage_top3",
                "percentage_top5",
                "percentage_top10",
                "percentage_found",
            ]
        )

    def write_row(self, directory: Path, rank_stats: RankStats) -> None:
        """Write summary rank stats row for run."""
        try:
            self.writer.writerow(
                [
                    directory,
                    rank_stats.top,
                    rank_stats.top3,
                    rank_stats.top5,
                    rank_stats.top10,
                    rank_stats.found,
                    rank_stats.total,
                    rank_stats.mean_reciprocal_rank(),
                    rank_stats.percentage_top(),
                    rank_stats.percentage_top3(),
                    rank_stats.percentage_top5(),
                    rank_stats.percentage_top10(),
                    rank_stats.percentage_found(),
                ]
            )
        except IOError:
            print("Error writing ", self.file)

    def close(self) -> None:
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)
