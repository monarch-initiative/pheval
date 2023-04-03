from dataclasses import dataclass
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from pheval.analyse.rank_stats import RankStats
from pheval.constants import PHEVAL_RESULTS_DIRECTORY_SUFFIX


def trim_corpus_results_directory_suffix(corpus_results_directory: Path) -> Path:
    """Trim the end of the corpus results directory name."""
    return Path(str(corpus_results_directory).replace(PHEVAL_RESULTS_DIRECTORY_SUFFIX, ""))


@dataclass
class TrackGenePrioritisation:
    """Track gene prioritisation for a run."""

    results_dir: Path
    ranks: dict
    rank_stats: RankStats


@dataclass
class TrackVariantPrioritisation:
    """Track variant prioritisation for a run."""

    results_dir: Path
    ranks: dict
    rank_stats: RankStats


@dataclass
class TrackPrioritisation:
    """Track prioritisation for a run."""

    gene_prioritisation: TrackGenePrioritisation
    variant_prioritisation: TrackVariantPrioritisation


class PlotGenerator:
    def __init__(self, gene_analysis: bool):
        self.gene_analysis = gene_analysis
        self.stats, self.mrr = [], []
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False

    def _retrieve_prioritisation_data(self, prioritisation_result: TrackPrioritisation):
        """Return either gene prioritisation or variant prioritisation stats."""
        return (
            prioritisation_result.gene_prioritisation
            if self.gene_analysis
            else prioritisation_result.variant_prioritisation
        )

    def _generate_stacked_bar_plot_data(self, prioritisation_result: TrackPrioritisation) -> None:
        """Generate data in correct format for dataframe creation for stacked bar plot."""
        result = self._retrieve_prioritisation_data(prioritisation_result)
        rank_stats = result.rank_stats
        self.stats.append(
            {
                "Run": f"{result.results_dir.parents[0].name}_"
                f"{trim_corpus_results_directory_suffix(result.results_dir.name)}",
                "Top": result.rank_stats.percentage_top(),
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
                "FO/NP": rank_stats.percentage_difference(100, rank_stats.percentage_found()),
            }
        )

    def _generate_stats_mrr_bar_plot_data(self, prioritisation_result: TrackPrioritisation) -> None:
        """Generate data in correct format for dataframe creation for MRR bar plot."""
        result = self._retrieve_prioritisation_data(prioritisation_result)
        self.mrr.extend(
            [
                {
                    "Rank": "MRR",
                    "Percentage": result.rank_stats.mean_reciprocal_rank(),
                    "Run": f"{result.results_dir.parents[0].name}_"
                    f"{trim_corpus_results_directory_suffix(result.results_dir.name)}",
                }
            ]
        )

    def generate_stacked_bar_gene(self, prioritisation_data: [TrackPrioritisation]) -> None:
        """Generate stacked bar plot and MRR bar plot for gene prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_stacked_bar_plot_data(prioritisation_result)
            self._generate_stats_mrr_bar_plot_data(prioritisation_result)
        gene_prioritisation_stats_df = pd.DataFrame(self.stats)
        gene_prioritisation_stats_df.set_index("Run").plot(
            kind="bar",
            stacked=True,
            colormap="tab10",
            ylabel="Disease-causing genes (%)",
            figsize=(10, 8),
            # rot=45,
        ).legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.savefig("gene_rank_stats.svg", format="svg", bbox_inches="tight")
        gene_mrr_df = pd.DataFrame(self.mrr)
        gene_mrr_df.set_index("Run").plot(
            kind="bar",
            colormap="tab10",
            ylabel="Gene prioritisation mean reciprocal rank",
            legend=False,
        )
        plt.savefig("gene_mrr.svg", format="svg", bbox_inches="tight")

    def generate_stacked_bar_variant(self, prioritisation_data: [TrackPrioritisation]):
        """Generate stacked bar plot and MRR bar plot for variant prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_stacked_bar_plot_data(prioritisation_result)
            self._generate_stats_mrr_bar_plot_data(prioritisation_result)
        variant_prioritisation_stats_df = pd.DataFrame(self.stats)

        variant_prioritisation_stats_df.set_index("Run").plot(
            kind="bar", stacked=True, colormap="tab10", ylabel="Disease-causing variants (%)"
        ).legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.savefig("variant_rank_stats.svg", format="svg", bbox_inches="tight")
        gene_mrr_df = pd.DataFrame(self.mrr)
        gene_mrr_df.set_index("Run").plot(
            kind="bar",
            colormap="tab10",
            ylabel="Variant prioritisation mean reciprocal rank",
            legend=False,
        )
        plt.savefig("variant_mrr.svg", format="svg", bbox_inches="tight")

    def _generate_cumulative_bar_plot_data(self, prioritisation_result: TrackPrioritisation):
        """Generate data in correct format for dataframe creation for cumulative bar plot."""
        result = self._retrieve_prioritisation_data(prioritisation_result)
        rank_stats = result.rank_stats
        trimmed_corpus_results_dir = trim_corpus_results_directory_suffix(result.results_dir.name)
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top3",
                    "Percentage": rank_stats.percentage_top3() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top5",
                    "Percentage": rank_stats.percentage_top5() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Top10",
                    "Percentage": rank_stats.percentage_top10() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "Found",
                    "Percentage": rank_stats.percentage_found() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "FO/NP",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.mean_reciprocal_rank(),
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
            ]
        )

    def generate_cumulative_bar_gene(self, prioritisation_data: [TrackPrioritisation]):
        """Generate cumulative bar plot for gene prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_cumulative_bar_plot_data(prioritisation_result)
        gene_prioritisation_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=gene_prioritisation_df, kind="bar", x="Rank", y="Percentage", hue="Run"
        ).set(xlabel="Rank", ylabel="Disease-causing genes (%)")
        plt.savefig("gene_rank_stats.svg", format="svg", bbox_inches="tight")

    def generate_cumulative_bar_variant(self, prioritisation_data: [TrackPrioritisation]):
        """Generate cumulative bar plot for variant prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_cumulative_bar_plot_data(prioritisation_result)
        variant_prioritisation_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=variant_prioritisation_df, kind="bar", x="Rank", y="Percentage", hue="Run"
        ).set(xlabel="Rank", ylabel="Disease-causing variants (%)")
        plt.savefig("variant_rank_stats.svg", format="svg", bbox_inches="tight")

    def _generate_non_cumulative_bar_plot_data(
        self, prioritisation_result: TrackPrioritisation
    ) -> [dict]:
        """Generate data in correct format for dataframe creation for non-cumulative bar plot."""
        result = self._retrieve_prioritisation_data(prioritisation_result)
        rank_stats = result.rank_stats
        trimmed_corpus_results_dir = trim_corpus_results_directory_suffix(result.results_dir.name)
        self.stats.extend(
            [
                {
                    "Rank": "Top",
                    "Percentage": rank_stats.percentage_top() / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "2-3",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top3(), rank_stats.percentage_top()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "4-5",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top5(), rank_stats.percentage_top3()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "6-10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_top10(), rank_stats.percentage_top5()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": ">10",
                    "Percentage": rank_stats.percentage_difference(
                        rank_stats.percentage_found(), rank_stats.percentage_top10()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "FO/NP",
                    "Percentage": rank_stats.percentage_difference(
                        100, rank_stats.percentage_found()
                    )
                    / 100,
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
                {
                    "Rank": "MRR",
                    "Percentage": rank_stats.mean_reciprocal_rank(),
                    "Run": f"{result.results_dir.parents[0].name}_" f"{trimmed_corpus_results_dir}",
                },
            ]
        )

    def generate_non_cumulative_bar_gene(self, prioritisation_data: [TrackPrioritisation]):
        """Generate non-cumulative bar plot for gene prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_non_cumulative_bar_plot_data(prioritisation_result)
        gene_prioritisation_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=gene_prioritisation_df, kind="bar", x="Rank", y="Percentage", hue="Run"
        ).set(xlabel="Rank", ylabel="Disease-causing genes (%)")
        plt.savefig("gene_rank_stats.svg", format="svg", bbox_inches="tight")

    def generate_non_cumulative_bar_variant(self, prioritisation_data: [TrackPrioritisation]):
        """Generate non-cumulative bar plot for variant prioritisation stats."""
        for prioritisation_result in prioritisation_data:
            self._generate_non_cumulative_bar_plot_data(prioritisation_result)
        variant_prioritisation_df = pd.DataFrame(self.stats)
        sns.catplot(
            data=variant_prioritisation_df, kind="bar", x="Rank", y="Percentage", hue="Run"
        ).set(xlabel="Rank", ylabel="Disease-causing variants (%)")
        plt.savefig("variant_rank_stats.svg", format="svg", bbox_inches="tight")


def generate_gene_plots(prioritisation_data: [TrackPrioritisation], plot_type: str) -> None:
    """Generate summary stats bar plot for gene prioritisation."""
    plot_generator = PlotGenerator(gene_analysis=True)
    if plot_type == "bar_stacked":
        plot_generator.generate_stacked_bar_gene(prioritisation_data)
    elif plot_type == "bar_cumulative":
        plot_generator.generate_cumulative_bar_gene(prioritisation_data)
    elif plot_type == "bar_non_cumulative":
        plot_generator.generate_non_cumulative_bar_gene(prioritisation_data)


def generate_variant_plots(prioritisation_data: [TrackPrioritisation], plot_type: str) -> None:
    """Generate summary stats bar plot for variant prioritisation."""
    plot_generator = PlotGenerator(gene_analysis=False)
    if plot_type == "bar_stacked":
        plot_generator.generate_stacked_bar_variant(prioritisation_data)
    elif plot_type == "bar_cumulative":
        plot_generator.generate_cumulative_bar_variant(prioritisation_data)
    elif plot_type == "bar_non_cumulative":
        plot_generator.generate_non_cumulative_bar_variant(prioritisation_data)
