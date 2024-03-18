import itertools
from collections import defaultdict
from copy import deepcopy
from typing import List

import numpy as np
import pandas as pd

from pheval.analyse.benchmark_generator import BenchmarkRunOutputGenerator
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.generate_plots import generate_plots
from pheval.constants import RANK_COMPARISON_FILE_SUFFIX


class RankComparisonGenerator:
    """Class for writing the run comparison of rank assignment for prioritisation."""

    def __init__(self, run_comparison: defaultdict):
        """
        Initialise the RankComparisonGenerator class.

        Args:
            run_comparison (defaultdict): A nested dictionary containing the run comparison data.
        """
        self.run_comparison = run_comparison

    def _generate_dataframe(self) -> pd.DataFrame:
        """
        Generate a Pandas DataFrame based on the run comparison data.

        Returns:
            pd.DataFrame: DataFrame containing the run comparison data.
        """
        return pd.DataFrame.from_dict(self.run_comparison, orient="index")

    def _calculate_rank_difference(self) -> pd.DataFrame:
        """
        Calculate the rank decrease for runs, taking the first directory as a baseline.

        Returns:
            pd.DataFrame: DataFrame containing the calculated rank differences.
        """
        comparison_df = self._generate_dataframe()
        comparison_df["rank_change"] = comparison_df.iloc[:, 2] - comparison_df.iloc[:, 3]
        comparison_df["rank_change"] = np.where(
            (comparison_df.iloc[:, 2] == 0) & (comparison_df.iloc[:, 3] != 0),
            "GAINED",
            np.where(
                (comparison_df.iloc[:, 3] == 0) & (comparison_df.iloc[:, 2] != 0),
                "LOST",
                comparison_df["rank_change"],
            ),
        )
        comparison_df["rank_change"] = comparison_df["rank_change"].apply(
            lambda x: int(x) if str(x).lstrip("-").isdigit() else x
        )
        return comparison_df

    def generate_output(self, prefix: str, suffix: str) -> None:
        """
        Generate output file from the run comparison data.

        Args:
            prefix (str): Prefix for the output file name.
            suffix (str): Suffix for the output file name.
        """
        self._generate_dataframe().to_csv(prefix + suffix, sep="\t")

    def generate_comparison_output(self, prefix: str, suffix: str) -> None:
        """
        Generate output file with calculated rank differences.

        Args:
            prefix (str): Prefix for the output file name.
            suffix (str): Suffix for the output file name.
        """
        self._calculate_rank_difference().to_csv(prefix + suffix, sep="\t")


def generate_benchmark_output(
    benchmarking_results: BenchmarkRunResults,
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Generate prioritisation outputs for a single benchmarking run.

    Args:
        benchmarking_results (BenchmarkRunResults): Results of a benchmarking run.
        plot_type (str): Type of plot to generate.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
    """
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


def merge_results(result1: dict, result2: dict) -> defaultdict:
    """
    Merge two nested dictionaries containing results on commonalities.

    This function merges two dictionaries, `result1` and `result2`, containing nested structures.
    It traverses the dictionaries recursively and merges their contents based on common keys.
    If a key is present in both dictionaries and points to another dictionary, the function
    will further merge their nested contents. If a key exists in `result2` but not in `result1`,
    it will be added to `result1`.

    Args:
        result1 (dict): The first dictionary to be merged.
        result2 (dict): The second dictionary to be merged.

    Returns:
        defaultdict: The merged dictionary containing the combined contents of `result1` and `result2`.
    """
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
    benchmarking_results: List[BenchmarkRunResults],
    plot_type: str,
    benchmark_generator: BenchmarkRunOutputGenerator,
) -> None:
    """
    Generate prioritisation outputs for benchmarking multiple runs.

    This function generates comparison outputs for benchmarking multiple runs. It compares the results
    between pairs of `BenchmarkRunResults` instances in `benchmarking_results` and generates rank
    comparison outputs using `RankComparisonGenerator` for each pair.

    Args:
        benchmarking_results (List[BenchmarkRunResults]): A list containing BenchmarkRunResults instances
            representing the benchmarking results of multiple runs.
        plot_type (str): The type of plot to be generated.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
    """
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
