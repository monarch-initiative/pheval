import csv
import itertools
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import pandas as pd

from pheval.analyse.benchmark_generator import BenchmarkPrioritisationOutputGenerator
from pheval.analyse.benchmarking_data import TrackRunPrioritisation
from pheval.analyse.generate_plots import generate_plots
from pheval.analyse.rank_stats import RankStats
from pheval.constants import RANK_COMPARISON_FILE_SUFFIX


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
        comparison_df["rank_decrease"] = comparison_df.iloc[:, 3] - comparison_df.iloc[:, 2]
        return comparison_df

    def generate_output(self, prefix: str, suffix: str) -> None:
        self._generate_dataframe().to_csv(prefix + suffix, sep="\t")

    def generate_comparison_output(self, prefix: str, suffix: str) -> None:
        self._calculate_rank_difference().to_csv(prefix + suffix, sep="\t")


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


def generate_benchmark_output(
    prioritisation_data: TrackRunPrioritisation,
    plot_type: str,
    benchmark_generator: BenchmarkPrioritisationOutputGenerator,
) -> None:
    """Generate prioritisation outputs for benchmarking single run."""
    rank_comparison_data = benchmark_generator.return_function(prioritisation_data).ranks
    results_dir_name = benchmark_generator.return_function(prioritisation_data).results_dir.name
    RankComparisonGenerator(rank_comparison_data).generate_output(
        f"{results_dir_name}",
        f"-{benchmark_generator.prioritisation_type_file_prefix}{RANK_COMPARISON_FILE_SUFFIX}",
    )
    generate_plots(
        [prioritisation_data],
        benchmark_generator,
        plot_type,
    )


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


def generate_benchmark_comparison_output(
    prioritisation_stats_for_runs: [TrackRunPrioritisation],
    plot_type: str,
    benchmark_generator: BenchmarkPrioritisationOutputGenerator,
) -> None:
    """Generate prioritisation outputs for benchmarking multiple runs."""
    output_prefix = benchmark_generator.prioritisation_type_file_prefix
    for pair in itertools.combinations(prioritisation_stats_for_runs, 2):
        result1 = benchmark_generator.return_function(pair[0])
        result2 = benchmark_generator.return_function(pair[1])
        merged_results = merge_results(
            deepcopy(result1.ranks),
            deepcopy(result2.ranks),
        )
        RankComparisonGenerator(merged_results).generate_comparison_output(
            f"{result1.results_dir.parents[0].name}_"
            f"{result1.results_dir.name}"
            f"_vs_{result2.results_dir.parents[0].name}_"
            f"{result2.results_dir.name}",
            f"-{output_prefix}{RANK_COMPARISON_FILE_SUFFIX}",
        )

    generate_plots(
        prioritisation_stats_for_runs,
        benchmark_generator,
        plot_type,
    )
