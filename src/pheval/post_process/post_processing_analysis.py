# #!/usr/bin/python
import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean

import pandas as pd

from ..utils.phenopacket_utils import VariantData


@dataclass
class GenePrioritisationResultData:
    phenopacket: Path
    gene: str
    rank: int = 0


@dataclass
class VariantPrioritisationResultData:
    phenopacket: Path
    variant: VariantData
    rank: int = 0


@dataclass
class ComparePrioritisationForRuns:
    """Compares the rank of different runs."""

    index: int
    directory: Path
    prioritisation_run_comparison: VariantPrioritisationResultData or GenePrioritisationResultData
    run_comparison: defaultdict

    def record_rank(self) -> None:
        self.run_comparison[self.index][
            "Phenopacket"
        ] = self.prioritisation_run_comparison.phenopacket.name
        if type(self.prioritisation_run_comparison) is GenePrioritisationResultData:
            self.run_comparison[self.index]["Gene"] = self.prioritisation_run_comparison.gene
        if type(self.prioritisation_run_comparison) is VariantPrioritisationResultData:
            variant_info = self.prioritisation_run_comparison.variant
            self.run_comparison[self.index]["Variant"] = "_".join(
                [variant_info.chrom, str(variant_info.pos), variant_info.ref, variant_info.alt]
            )
        self.run_comparison[self.index][self.directory] = self.prioritisation_run_comparison.rank


class GenerateRankComparisonOutput:
    """Writes the run comparison of rank assignment for prioritisation."""

    def __init__(self, run_comparison: defaultdict):
        self.run_comparison = run_comparison

    def generate_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame.from_dict(self.run_comparison, orient="index")

    def calculate_rank_difference(self) -> pd.DataFrame:
        comparison_df = self.generate_dataframe()
        comparison_df["rank_decrease"] = comparison_df.iloc[:, 3] - comparison_df.iloc[:, 2]
        return comparison_df

    def generate_gene_output(self, prefix: str) -> None:
        # comparison_df["absolute_rank_difference"] = pd.Series.abs(
        #     comparison_df[:, 2] - comparison_df[:, 3]
        # )
        self.calculate_rank_difference().to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def generate_variant_output(self, prefix: str) -> None:
        # comparison_df["absolute_rank_difference"] = pd.Series.abs(
        #     comparison_df[:, 2] - comparison_df[:, 3]
        # )
        self.calculate_rank_difference().to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")


@dataclass
class RankStats:
    """Class for keeping track of the rank stats."""

    top: int = 0
    top3: int = 0
    top5: int = 0
    found: int = 0
    total: int = 0
    reciprocal_ranks: list = field(default_factory=list)

    def add_rank(self, rank: int) -> None:
        self.reciprocal_ranks.append(1 / rank)
        self.found += 1
        if rank == 1:
            self.top += 1
        if rank != "" and rank <= 3:
            self.top3 += 1
        if rank != "" and rank <= 5:
            self.top5 += 1

    def percentage_rank(self, value: int) -> float:
        return 100 * value / self.found

    def percentage_top(self) -> float:
        return self.percentage_rank(self.top)

    def percentage_top3(self) -> float:
        return self.percentage_rank(self.top3)

    def percentage_top5(self) -> float:
        return self.percentage_rank(self.top5)

    def percentage_found(self) -> float:
        return 100 * self.found / self.total

    def mean_reciprocal_rank(self) -> float:
        return mean(self.reciprocal_ranks)


class RankStatsWriter:
    """Writes the rank stats for each run."""

    def __init__(self, file: Path):
        self.file = open(file, "w")
        self.writer = csv.writer(self.file, delimiter="\t")
        self.writer.writerow(
            [
                "results_directory_path",
                "top",
                "top3",
                "top5",
                "found",
                "total",
                "mean_reciprocal_rank",
                "percentage_top",
                "percentage_top3",
                "percentage_top5",
                "percentage_found",
            ]
        )

    def write_row(self, directory: Path, rank_stats: RankStats) -> None:
        try:
            self.writer.writerow(
                [
                    directory,
                    rank_stats.top,
                    rank_stats.top3,
                    rank_stats.top5,
                    rank_stats.found,
                    rank_stats.total,
                    rank_stats.mean_reciprocal_rank(),
                    rank_stats.percentage_top(),
                    rank_stats.percentage_top3(),
                    rank_stats.percentage_top5(),
                    rank_stats.percentage_found(),
                ]
            )
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)
