import csv
import itertools
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import pandas as pd

from pheval.analyse.generate_plots import (
    TrackPrioritisation,
    generate_gene_plots,
    generate_variant_plots,
)
from pheval.analyse.rank_stats import RankStats


class RankComparisonGenerator:
    """Write the run comparison of rank assignment for prioritisation."""

    def __init__(self, run_comparison: defaultdict):
        self.run_comparison = run_comparison

    def _generate_dataframe(self) -> pd.DataFrame:
        """Generate pandas dataframe."""
        return pd.DataFrame.from_dict(self.run_comparison, orient="index")

    def _calculate_rank_difference(self) -> pd.DataFrame:
        """Calculate the rank decrease for runs - taking the first directory as a baseline."""
        comparison_df = self._generate_dataframe()
        print(len(comparison_df.columns))
        comparison_df["rank_decrease"] = comparison_df.iloc[:, 3] - comparison_df.iloc[:, 2]
        return comparison_df

    def generate_gene_output(self, prefix: str) -> None:
        """Generate the output for gene prioritisation ranks."""
        self._generate_dataframe().to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def generate_variant_output(self, prefix: str) -> None:
        """Generate the output for variant prioritisation ranks."""
        self._generate_dataframe().to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")

    def generate_gene_comparison_output(self, prefix: str) -> None:
        """Generate the output for gene prioritisation rank comparison."""
        self._calculate_rank_difference().to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def generate_variant_comparison_output(self, prefix: str) -> None:
        """Generate the output for variant prioritisation rank comparison."""
        self._calculate_rank_difference().to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")


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


def generate_benchmark_gene_output(
    prioritisation_data: TrackPrioritisation, plot_type: str
) -> None:
    """Generate gene prioritisation outputs for benchmarking single run."""
    RankComparisonGenerator(prioritisation_data.gene_prioritisation.ranks).generate_gene_output(
        f"{prioritisation_data.gene_prioritisation.results_dir.name}"
    )
    generate_gene_plots([prioritisation_data], plot_type)


def generate_benchmark_variant_output(
    prioritisation_data: TrackPrioritisation, plot_type: str
) -> None:
    """Generate variant prioritisation outputs for benchmarking single run."""
    RankComparisonGenerator(
        prioritisation_data.variant_prioritisation.ranks
    ).generate_variant_output(f"{prioritisation_data.gene_prioritisation.results_dir.name}")
    generate_variant_plots([prioritisation_data], plot_type)


def merge_results(result1: dict, result2: dict) -> dict:
    """Merge two nested dictionaries containing results on commonalities."""
    for key, val in result1.items():
        if type(val) == dict:
            if key in result2 and type(result2[key] == dict):
                merge_results(result1[key], result2[key])
        else:
            if key in result2:
                result1[key] = result2[key]

    for key, val in result2.items():
        if key not in result1:
            result1[key] = val
    return result1


def generate_gene_rank_comparisons(comparison_ranks: [tuple]) -> None:
    """Generate the gene rank comparison of two result directories."""
    for pair in comparison_ranks:
        merged_results = merge_results(
            deepcopy(pair[0].gene_prioritisation.ranks), deepcopy(pair[1].gene_prioritisation.ranks)
        )
        RankComparisonGenerator(merged_results).generate_gene_comparison_output(
            f"{pair[0].gene_prioritisation.results_dir.parents[0].name}_"
            f"{pair[0].gene_prioritisation.results_dir.name}"
            f"__v__{pair[1].gene_prioritisation.results_dir.parents[0].name}_"
            f"{pair[1].gene_prioritisation.results_dir.name}"
        )


def generate_variant_rank_comparisons(comparison_ranks: [tuple]) -> None:
    """Generate the variant rank comparison of two result directories."""
    for pair in comparison_ranks:
        merged_results = merge_results(
            deepcopy(pair[0].variant_prioritisation.ranks),
            deepcopy(pair[1].variant_prioritisation.ranks),
        )
        RankComparisonGenerator(merged_results).generate_variant_comparison_output(
            f"{pair[0].gene_prioritisation.results_dir.parents[0].name}_"
            f"{pair[0].variant_prioritisation.results_dir.name}"
            f"__v__{pair[0].gene_prioritisation.results_dir.parents[0].name}_"
            f"{pair[1].variant_prioritisation.results_dir.name}"
        )


def generate_benchmark_comparison_gene_output(
    prioritisation_stats_for_runs: [TrackPrioritisation], plot_type: str
) -> None:
    """Generate gene prioritisation outputs for benchmarking multiple runs."""
    generate_gene_rank_comparisons(list(itertools.combinations(prioritisation_stats_for_runs, 2)))
    generate_gene_plots(prioritisation_stats_for_runs, plot_type)


def generate_benchmark_comparison_variant_output(
    prioritisation_stats_for_runs: [TrackPrioritisation], plot_type: str
) -> None:
    """Generate variant prioritisation outputs for benchmarking multiple runs."""
    generate_variant_rank_comparisons(
        list(itertools.combinations(prioritisation_stats_for_runs, 2))
    )
    generate_variant_plots(prioritisation_stats_for_runs, plot_type)
