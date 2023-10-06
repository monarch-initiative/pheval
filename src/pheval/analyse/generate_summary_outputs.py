import itertools
from collections import defaultdict
from copy import deepcopy

import pandas as pd

from pheval.analyse.benchmark_generator import BenchmarkRunOutputGenerator
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.generate_plots import generate_plots
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


def generate_benchmark_output(
    benchmarking_results: BenchmarkRunResults,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Generate prioritisation outputs for benchmarking single run."""
    rank_comparison_data = benchmarking_results.ranks
    results_dir_name = benchmarking_results.results_dir.name
    RankComparisonGenerator(rank_comparison_data).generate_output(
        f"{results_dir_name}",
        f"-{benchmark_generator.prioritisation_type_file_prefix}{RANK_COMPARISON_FILE_SUFFIX}",
    )
    generate_plots(
        [benchmarking_results],
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
    benchmarking_results: [BenchmarkRunResults],
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """Generate prioritisation outputs for benchmarking multiple runs."""
    output_prefix = benchmark_generator.prioritisation_type_file_prefix
    for pair in itertools.combinations(benchmarking_results, 2):
        result1 = pair[0]
        result2 = pair[1]
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
        benchmarking_results,
        benchmark_generator,
        plot_type,
    )
