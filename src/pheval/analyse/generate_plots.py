from pathlib import Path
from typing import List

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.metrics import auc, precision_recall_curve, roc_curve

from pheval.analyse.benchmark_generator import (
    BenchmarkRunOutputGenerator,
    DiseaseBenchmarkRunOutputGenerator,
    GeneBenchmarkRunOutputGenerator,
    VariantBenchmarkRunOutputGenerator,
)
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.parse_benchmark_summary import (
    parse_benchmark_result_summary,
    read_benchmark_tsv_result_summary,
)
from pheval.constants import PHEVAL_RESULTS_DIRECTORY_SUFFIX


def trim_corpus_results_directory_suffix(corpus_results_directory: Path) -> Path:
    """
    Trim the suffix from the corpus results directory name.

    Args:
        corpus_results_directory (Path): The directory path containing corpus results.

    Returns:
        Path: The Path object with the suffix removed from the directory name.
    """
    return Path(str(corpus_results_directory).replace(PHEVAL_RESULTS_DIRECTORY_SUFFIX, ""))


class PlotGenerator:
    """Class to generate plots."""

    palette_hex_codes = [
        "#f4ae3d",
        "#ee5825",
        "#2b7288",
        "#9a84b2",
        "#0c604c",
        "#c94c4c",
        "#3d8e83",
        "#725ac1",
        "#e7ba52",
        "#1b9e77",
    ]

    def __init__(
        self,
    ):
        """
        Initialise the PlotGenerator class.
        Note:
            `self.stats` will be used to store statistics data.
            `self.mrr` will store Mean Reciprocal Rank (MRR) values.
            Matplotlib settings are configured to remove the right and top axes spines
            for generated plots.
        """
        self.stats, self.mrr = [], []
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False

    @staticmethod
    def _create_run_identifier(results_dir: Path) -> str:
        """
        Create a run identifier from a path.

        Args:
            results_dir (Path): The directory path for results.

        Returns:
            str: A string representing the run identifier created from the given path.
        """
        return f"{Path(results_dir).parents[0].name}_{trim_corpus_results_directory_suffix(Path(results_dir).name)}"

    def return_benchmark_name(self, benchmark_result: BenchmarkRunResults) -> str:
        """
        Return the benchmark name for a run.

        Args:
            benchmark_result (BenchmarkRunResults): The benchmarking results for a run.

        Returns:
            str: The benchmark name obtained from the given BenchmarkRunResults instance.
        """
        return (
            benchmark_result.benchmark_name
            if benchmark_result.results_dir is None
            else self._create_run_identifier(benchmark_result.results_dir)
        )

    def _generate_stacked_bar_plot_data(self, benchmark_result: BenchmarkRunResults) -> None:
        """
        Generate data in the correct format for dataframe creation for a stacked bar plot,
        appending to the self.stats attribute of the class.

        Args:
            benchmark_result (BenchmarkRunResults): The benchmarking results for a run.
        """
        rank_stats = benchmark_result.rank_stats
        self.stats.append(
            {
                "Run": self.return_benchmark_name(benchmark_result),
                "Top": benchmark_result.rank_stats.percentage_top(),
                "2-3": rank_stats.percentage_difference(
                    rank_stats.percentage_top3(), rank_stats.percentage_top()
                ),
                "4-5": rank_stats.percentage_difference(
                    rank_stats.percentage_top5(), rank_stats.percentage_top3()
                ),
                "6-10": rank_stats.percentage_difference(
                    rank_stats.percentage_top10(), rank_stats.percentage_top5()
                ),
                ">10": rank_stats.percentage_difference(
                    rank_stats.percentage_found(), rank_stats.percentage_top10()
                ),
                "Missed": rank_stats.percentage_difference(100, rank_stats.percentage_found()),
            }
        )

    def _generate_stats_mrr_bar_plot_data(self, benchmark_result: BenchmarkRunResults) -> None:
        """
        Generate data in the correct format for dataframe creation for MRR (Mean Reciprocal Rank) bar plot,
        appending to the self.mrr attribute of the class.

        Args:
            benchmark_result (BenchmarkRunResults): The benchmarking results for a run.
        """
        self.mrr.extend(
            [
                {
                    "Rank": "MRR",
                    "Percentage": benchmark_result.rank_stats.return_mean_reciprocal_rank(),
                    "Run": self.return_benchmark_name(benchmark_result),
                }
            ]
        )

    def generate_stacked_bar_plot(
        self,
        benchmarking_results: List[BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
        title: str = None,
    ) -> None:
        """
        Generate a stacked bar plot and Mean Reciprocal Rank (MRR) bar plot.

        Args:
            benchmarking_results (List[BenchmarkRunResults]): List of benchmarking results for multiple runs.
            benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
            title (str, optional): Title for the generated plot. Defaults to None.
        """
        for benchmark_result in benchmarking_results:
            self._generate_stacked_bar_plot_data(benchmark_result)
            self._generate_stats_mrr_bar_plot_data(benchmark_result)
        stats_df = pd.DataFrame(self.stats)
        stats_df.set_index("Run").plot(
            kind="bar",
            stacked=True,
            color=self.palette_hex_codes,
            ylabel=benchmark_generator.y_label,
            edgecolor="white",
        ).legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        if title is None:
            plt.title(
                f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} Rank Stats"
            )
        else:
            plt.title(title, loc="center", fontsize=15)
        plt.ylim(0, 100)
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )

        mrr_df = pd.DataFrame(self.mrr)
        mrr_df.set_index("Run").plot(
            kind="bar",
            color=self.palette_hex_codes,
            ylabel=f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} mean reciprocal rank",
            legend=False,
            edgecolor="white",
        )
        plt.title(
            f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} results - mean reciprocal rank"
        )
        plt.ylim(0, 1)
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_mrr.svg",
            format="svg",
            bbox_inches="tight",
        )

    def _generate_cumulative_bar_plot_data(self, benchmark_result: BenchmarkRunResults):
        """
        Generate data in the correct format for dataframe creation for a cumulative bar plot,
        appending to the self.stats attribute of the class.

        Args:
            benchmark_result (BenchmarkRunResults): The benchmarking results for a run.
        """
        rank_stats = benchmark_result.rank_stats
        run_identifier = self.return_benchmark_name(benchmark_result)
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Top3",
                    "Percentage": rank_stats.percentage_top3() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Top5",
                    "Percentage": rank_stats.percentage_top5() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Top10",
                    "Percentage": rank_stats.percentage_top10() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Found",
                    "Percentage": rank_stats.percentage_found() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Missed",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.return_mean_reciprocal_rank(),
                    "Run": run_identifier,
                },
            ]
        )

    def generate_cumulative_bar(
        self,
        benchmarking_results: List[BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
        title: str = None,
    ) -> None:
        """
        Generate a cumulative bar plot.

        Args:
            benchmarking_results (List[BenchmarkRunResults]): List of benchmarking results for multiple runs.
            benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
            title (str, optional): Title for the generated plot. Defaults to None.
        """
        for benchmark_result in benchmarking_results:
            self._generate_cumulative_bar_plot_data(benchmark_result)
        stats_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=stats_df,
            kind="bar",
            x="Rank",
            y="Percentage",
            hue="Run",
            palette=self.palette_hex_codes,
            edgecolor="white",
            legend=False,
        ).set(xlabel="Rank", ylabel=benchmark_generator.y_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, title="Run")
        if title is None:
            plt.title(
                f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} Cumulative Rank Stats"
            )
        else:
            plt.title(title, loc="center", fontsize=15)
        plt.ylim(0, 1)
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )

    def _generate_non_cumulative_bar_plot_data(
        self, benchmark_result: BenchmarkRunResults
    ) -> [dict]:
        """
        Generate data in the correct format for dataframe creation for a non-cumulative bar plot,
        appending to the self.stats attribute of the class.

        Args:
            benchmark_result (BenchmarkRunResults): The benchmarking results for a run.
        """
        rank_stats = benchmark_result.rank_stats
        run_identifier = self.return_benchmark_name(benchmark_result)
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "2-3",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top3(), rank_stats.percentage_top()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "4-5",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top5(), rank_stats.percentage_top3()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "6-10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top10(), rank_stats.percentage_top5()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": ">10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_found(), rank_stats.percentage_top10()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "Missed",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": run_identifier,
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.return_mean_reciprocal_rank(),
                    "Run": run_identifier,
                },
            ]
        )

    def generate_roc_curve(
        self,
        benchmarking_results: List[BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
    ):
        """
        Generate and plot Receiver Operating Characteristic (ROC) curves for binary classification benchmark results.

        Args:
            benchmarking_results (List[BenchmarkRunResults]): List of benchmarking results for multiple runs.
            benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
        """
        for i, benchmark_result in enumerate(benchmarking_results):
            fpr, tpr, thresh = roc_curve(
                benchmark_result.binary_classification_stats.labels,
                benchmark_result.binary_classification_stats.scores,
                pos_label=1,
            )
            roc_auc = auc(fpr, tpr)

            plt.plot(
                fpr,
                tpr,
                label=f"{self.return_benchmark_name(benchmark_result)} ROC Curve (AUC = {roc_auc:.2f})",
                color=self.palette_hex_codes[i],
            )

        plt.plot(linestyle="--", color="gray")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("Receiver Operating Characteristic (ROC) Curve")
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_roc_curve.svg",
            format="svg",
            bbox_inches="tight",
        )

    def generate_precision_recall(
        self,
        benchmarking_results: List[BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
    ):
        """
        Generate and plot Precision-Recall curves for binary classification benchmark results.

        Args:
            benchmarking_results (List[BenchmarkRunResults]): List of benchmarking results for multiple runs.
            benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
        """
        plt.figure()
        for i, benchmark_result in enumerate(benchmarking_results):
            precision, recall, thresh = precision_recall_curve(
                benchmark_result.binary_classification_stats.labels,
                benchmark_result.binary_classification_stats.scores,
            )
            precision_recall_auc = auc(recall, precision)
            plt.plot(
                recall,
                precision,
                label=f"{self.return_benchmark_name(benchmark_result)} Precision-Recall Curve "
                f"(AUC = {precision_recall_auc:.2f})",
                color=self.palette_hex_codes[i],
            )

        plt.plot(linestyle="--", color="gray")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title("Precision-Recall Curve")
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_precision_recall_curve.svg",
            format="svg",
            bbox_inches="tight",
        )

    def generate_non_cumulative_bar(
        self,
        benchmarking_results: List[BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
        title: str = None,
    ) -> None:
        """
        Generate a non-cumulative bar plot.

        Args:
            benchmarking_results (List[BenchmarkRunResults]): List of benchmarking results for multiple runs.
            benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
            title (str, optional): Title for the generated plot. Defaults to None.
        """
        for benchmark_result in benchmarking_results:
            self._generate_non_cumulative_bar_plot_data(benchmark_result)

        stats_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=stats_df,
            kind="bar",
            x="Rank",
            y="Percentage",
            hue="Run",
            palette=self.palette_hex_codes,
            edgecolor="white",
            legend=False,
        ).set(xlabel="Rank", ylabel=benchmark_generator.y_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, title="Run")
        if title is None:
            plt.title(
                f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} Non-Cumulative Rank Stats"
            )
        else:
            plt.title(title, loc="center", fontsize=15)
        plt.ylim(0, 1)
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )


def generate_plots(
    benchmarking_results: List[BenchmarkRunResults],
    benchmark_generator: BenchmarkRunOutputGenerator,
    plot_type: str,
    title: str = None,
) -> None:
    """
    Generate summary statistics bar plots for prioritisation.

    This method generates summary statistics bar plots based on the provided benchmarking results and plot type.

    Args:
        benchmarking_results (list[BenchmarkRunResults]): List of benchmarking results for multiple runs.
        benchmark_generator (BenchmarkRunOutputGenerator): Object containing benchmarking output generation details.
        plot_type (str): Type of plot to be generated ("bar_stacked", "bar_cumulative", "bar_non_cumulative").
        title (str, optional): Title for the generated plot. Defaults to None.
    """
    plot_generator = PlotGenerator()
    plot_generator.generate_roc_curve(benchmarking_results, benchmark_generator)
    plot_generator.generate_precision_recall(benchmarking_results, benchmark_generator)
    if plot_type == "bar_stacked":
        plot_generator.generate_stacked_bar_plot(benchmarking_results, benchmark_generator, title)
    elif plot_type == "bar_cumulative":
        plot_generator.generate_cumulative_bar(benchmarking_results, benchmark_generator, title)
    elif plot_type == "bar_non_cumulative":
        plot_generator.generate_non_cumulative_bar(benchmarking_results, benchmark_generator, title)


def generate_plots_from_benchmark_summary_tsv(
    benchmark_summary_tsv: Path,
    gene_analysis: bool,
    variant_analysis: bool,
    disease_analysis: bool,
    plot_type: str,
    title: str,
):
    """
    Generate bar plot from summary benchmark results.

    Reads a summary of benchmark results from a TSV file and generates a bar plot
    based on the analysis type and plot type.

    Args:
        benchmark_summary_tsv (Path): Path to the summary TSV file containing benchmark results.
        gene_analysis (bool): Flag indicating whether to analyse gene prioritisation.
        variant_analysis (bool): Flag indicating whether to analyse variant prioritisation.
        disease_analysis (bool): Flag indicating whether to analyse disease prioritisation.
        plot_type (str): Type of plot to be generated ("bar_stacked", "bar_cumulative", "bar_non_cumulative").
        title (str): Title for the generated plot.
    Raises:
         ValueError: If an unsupported plot type is specified.
    """
    benchmark_stats_summary = read_benchmark_tsv_result_summary(benchmark_summary_tsv)
    benchmarking_results = parse_benchmark_result_summary(benchmark_stats_summary)
    if gene_analysis:
        benchmark_generator = GeneBenchmarkRunOutputGenerator()
    elif variant_analysis:
        benchmark_generator = VariantBenchmarkRunOutputGenerator()
    elif disease_analysis:
        benchmark_generator = DiseaseBenchmarkRunOutputGenerator()
    else:
        raise ValueError(
            "Specify one analysis type (gene_analysis, variant_analysis, or disease_analysis)"
        )
    generate_plots(benchmarking_results, benchmark_generator, plot_type, title)
