from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pheval.analyse.benchmark_generator import BenchmarkPrioritisationOutputGenerator
from pheval.analyse.benchmarking_data import AnalysisResults, TrackRunPrioritisation
from pheval.constants import PHEVAL_RESULTS_DIRECTORY_SUFFIX


def trim_corpus_results_directory_suffix(corpus_results_directory: Path) -> Path:
    """Trim the end of the corpus results directory name."""
    return Path(str(corpus_results_directory).replace(PHEVAL_RESULTS_DIRECTORY_SUFFIX, ""))


class PlotGenerator:
    def __init__(
        self,
    ):
        self.stats, self.mrr = [], []
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False

    def _generate_stacked_bar_plot_data(self, prioritisation_result: AnalysisResults) -> None:
        """Generate data in correct format for dataframe creation for stacked bar plot."""
        rank_stats = prioritisation_result.rank_stats
        self.stats.append(
            {
                "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                f"{trim_corpus_results_directory_suffix(prioritisation_result.results_dir.name)}",
                "Top": prioritisation_result.rank_stats.percentage_top(),
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

    def _generate_stats_mrr_bar_plot_data(self, prioritisation_result: AnalysisResults) -> None:
        """Generate data in correct format for dataframe creation for MRR bar plot."""
        self.mrr.extend(
            [
                {
                    "Rank": "MRR",
                    "Percentage": prioritisation_result.rank_stats.mean_reciprocal_rank(),
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trim_corpus_results_directory_suffix(prioritisation_result.results_dir.name)}",
                }
            ]
        )

    def generate_stacked_bar_plot(
        self,
        prioritisation_data: [TrackRunPrioritisation],
        benchmark_generator: BenchmarkPrioritisationOutputGenerator,
    ) -> None:
        """Generate stacked bar plot."""
        for prioritisation_result in prioritisation_data:
            self._generate_stacked_bar_plot_data(
                benchmark_generator.return_function(prioritisation_result)
            )
            self._generate_stats_mrr_bar_plot_data(
                benchmark_generator.return_function(prioritisation_result)
            )
        stats_df = pd.DataFrame(self.stats)
        stats_df.set_index("Run").plot(
            kind="bar", stacked=True, colormap="tab10", ylabel=benchmark_generator.y_label
        ).legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )

        mrr_df = pd.DataFrame(self.mrr)
        mrr_df.set_index("Run").plot(
            kind="bar",
            colormap="tab10",
            ylabel=f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} mean reciprocal rank",
            legend=False,
        )
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_mrr.svg",
            format="svg",
            bbox_inches="tight",
        )

    def _generate_cumulative_bar_plot_data(self, prioritisation_result: AnalysisResults):
        """Generate data in correct format for dataframe creation for cumulative bar plot."""
        rank_stats = prioritisation_result.rank_stats
        trimmed_corpus_results_dir = trim_corpus_results_directory_suffix(
            prioritisation_result.results_dir.name
        )
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top3",
                    "Percentage": rank_stats.percentage_top3() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top5",
                    "Percentage": rank_stats.percentage_top5() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top10",
                    "Percentage": rank_stats.percentage_top10() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Found",
                    "Percentage": rank_stats.percentage_found() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Missed",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.mean_reciprocal_rank(),
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
            ]
        )

    def generate_cumulative_bar(
        self,
        prioritisation_data: [TrackRunPrioritisation],
        benchmark_generator: BenchmarkPrioritisationOutputGenerator,
    ) -> None:
        """Generate cumulative bar plot."""
        for prioritisation_result in prioritisation_data:
            self._generate_cumulative_bar_plot_data(
                benchmark_generator.return_function(prioritisation_result)
            )
        stats_df = pd.DataFrame(self.stats)
        sns.catplot(data=stats_df, kind="bar", x="Rank", y="Percentage", hue="Run").set(
            xlabel="Rank", ylabel=benchmark_generator.y_label
        )
        plt.title(
            f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} Cumulative Rank Stats"
        )
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )

    def _generate_non_cumulative_bar_plot_data(
        self, prioritisation_result: AnalysisResults
    ) -> [dict]:
        """Generate data in correct format for dataframe creation for non-cumulative bar plot."""
        rank_stats = prioritisation_result.rank_stats
        trimmed_corpus_results_dir = trim_corpus_results_directory_suffix(
            prioritisation_result.results_dir.name
        )
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "2-3",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top3(), rank_stats.percentage_top()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "4-5",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top5(), rank_stats.percentage_top3()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "6-10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top10(), rank_stats.percentage_top5()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": ">10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_found(), rank_stats.percentage_top10()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Missed",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.mean_reciprocal_rank(),
                    "Run": f"{prioritisation_result.results_dir.parents[0].name}_"
                    f"{trimmed_corpus_results_dir}",
                },
            ]
        )

    def generate_non_cumulative_bar(
        self,
        prioritisation_data: [TrackRunPrioritisation],
        benchmark_generator: BenchmarkPrioritisationOutputGenerator,
    ) -> None:
        """Generate non-cumulative bar plot."""
        for prioritisation_result in prioritisation_data:
            self._generate_non_cumulative_bar_plot_data(
                benchmark_generator.return_function(prioritisation_result)
            )

        stats_df = pd.DataFrame(self.stats)
        sns.catplot(data=stats_df, kind="bar", x="Rank", y="Percentage", hue="Run").set(
            xlabel="Rank", ylabel=benchmark_generator.y_label
        )
        plt.title(
            f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} Non-Cumulative Rank Stats"
        )
        plt.savefig(
            f"{benchmark_generator.prioritisation_type_file_prefix}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )


def generate_plots(
    prioritisation_data: [TrackRunPrioritisation],
    benchmark_generator: BenchmarkPrioritisationOutputGenerator,
    plot_type: str,
) -> None:
    """Generate summary stats bar plots for prioritisation."""
    plot_generator = PlotGenerator()
    if plot_type == "bar_stacked":
        plot_generator.generate_stacked_bar_plot(
            prioritisation_data,
            benchmark_generator,
        )
    elif plot_type == "bar_cumulative":
        plot_generator.generate_cumulative_bar(
            prioritisation_data,
            benchmark_generator,
        )
    elif plot_type == "bar_non_cumulative":
        plot_generator.generate_non_cumulative_bar(
            prioritisation_data,
            benchmark_generator,
        )
