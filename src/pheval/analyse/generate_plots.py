from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pheval.analyse.benchmark_generator import (
    BenchmarkRunOutputGenerator, GeneBenchmarkRunOutputGenerator,
    VariantBenchmarkRunOutputGenerator, DiseaseBenchmarkRunOutputGenerator
)
from pheval.analyse.benchmarking_data import BenchmarkRunResults
from pheval.analyse.parse_benchmark_summary import read_benchmark_tsv_result_summary, parse_benchmark_result_summary
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

    @staticmethod
    def create_run_identifier(results_dir: Path) -> str:
        """Create a run identifier from a path."""
        if Path(results_dir).exists():
            return f"{Path(results_dir).parents[0].name}_{trim_corpus_results_directory_suffix(Path(results_dir).name)}"
        return results_dir

    def _generate_stacked_bar_plot_data(self, benchmark_result: BenchmarkRunResults) -> None:
        """Generate data in correct format for dataframe creation for stacked bar plot."""
        rank_stats = benchmark_result.rank_stats
        self.stats.append(
            {
                "Run": self.create_run_identifier(benchmark_result.results_dir),
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
        """Generate data in correct format for dataframe creation for MRR bar plot."""
        self.mrr.extend(
            [
                {
                    "Rank": "MRR",
                    "Percentage": benchmark_result.rank_stats.return_mean_reciprocal_rank(),
                    "Run": self.create_run_identifier(benchmark_result.results_dir),

                }
            ]
        )

    def generate_stacked_bar_plot(
            self,
            benchmarking_results: [BenchmarkRunResults],
            benchmark_generator: BenchmarkRunOutputGenerator,
            title: str = None
    ) -> None:
        """Generate stacked bar plot."""
        for benchmark_result in benchmarking_results:
            self._generate_stacked_bar_plot_data(benchmark_result)
            self._generate_stats_mrr_bar_plot_data(benchmark_result)
        stats_df = pd.DataFrame(self.stats)
        stats_df.set_index("Run").plot(
            kind="bar", stacked=True, colormap="tab10", ylabel=benchmark_generator.y_label
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
            colormap="tab10",
            ylabel=f"{benchmark_generator.prioritisation_type_file_prefix.capitalize()} mean reciprocal rank",
            legend=False,
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
        """Generate data in correct format for dataframe creation for cumulative bar plot."""
        rank_stats = benchmark_result.rank_stats
        run_identifier = self.create_run_identifier(benchmark_result.results_dir)
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
            benchmarking_results: [BenchmarkRunResults],
            benchmark_generator: BenchmarkRunOutputGenerator,
            title: str = None
    ) -> None:
        """Generate cumulative bar plot."""
        for benchmark_result in benchmarking_results:
            self._generate_cumulative_bar_plot_data(benchmark_result)
        stats_df = pd.DataFrame(self.stats)
        sns.catplot(data=stats_df, kind="bar", x="Rank", y="Percentage", hue="Run").set(
            xlabel="Rank", ylabel=benchmark_generator.y_label
        )
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
        """Generate data in correct format for dataframe creation for non-cumulative bar plot."""
        rank_stats = benchmark_result.rank_stats
        run_identifier = self.create_run_identifier(benchmark_result.results_dir)
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

    def generate_non_cumulative_bar(
            self,
            benchmarking_results: [BenchmarkRunResults],
            benchmark_generator: BenchmarkRunOutputGenerator,
            title: str = None
    ) -> None:
        """Generate non-cumulative bar plot."""
        for benchmark_result in benchmarking_results:
            self._generate_non_cumulative_bar_plot_data(benchmark_result)

        stats_df = pd.DataFrame(self.stats)
        sns.catplot(data=stats_df, kind="bar", x="Rank", y="Percentage", hue="Run").set(
            xlabel="Rank", ylabel=benchmark_generator.y_label
        )
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
        benchmarking_results: [BenchmarkRunResults],
        benchmark_generator: BenchmarkRunOutputGenerator,
        plot_type: str,
        title: str = None
) -> None:
    """Generate summary stats bar plots for prioritisation."""
    plot_generator = PlotGenerator()
    if plot_type == "bar_stacked":
        plot_generator.generate_stacked_bar_plot(
            benchmarking_results,
            benchmark_generator,
            title
        )
    elif plot_type == "bar_cumulative":
        plot_generator.generate_cumulative_bar(
            benchmarking_results,
            benchmark_generator,
            title
        )
    elif plot_type == "bar_non_cumulative":
        plot_generator.generate_non_cumulative_bar(
            benchmarking_results,
            benchmark_generator,
            title
        )


def generate_plots_from_benchmark_summary_tsv(
        benchmark_summary_tsv: Path,
        gene_analysis: bool,
        variant_analysis: bool,
        disease_analysis: bool,
        plot_type: str,
        title: str,
):
    """Generate bar plot from summary benchmark results."""
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
