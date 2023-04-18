from dataclasses import dataclass, field
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
        return 100 * value / self.found

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
        return 100 * self.found / self.total

    @staticmethod
    def percentage_difference(percentage_value_1: float, percentage_value_2: float) -> float:
        """Return percentage difference between two percentage values"""
        return percentage_value_1 - percentage_value_2

    def mean_reciprocal_rank(self) -> float:
        """Return the mean reciprocal rank."""
        return mean(self.reciprocal_ranks)
