import unittest
from pathlib import Path, PosixPath

from pheval.analyse.generate_plots import (
    PlotGenerator,
    TrackGenePrioritisation,
    TrackPrioritisation,
    TrackVariantPrioritisation,
)
from pheval.analyse.rank_stats import RankStats


class TestPlotGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.gene_plot_generator = PlotGenerator(gene_analysis=True)
        self.variant_plot_generator = PlotGenerator(gene_analysis=False)
        self.track_prioritisation = TrackPrioritisation(
            gene_prioritisation=TrackGenePrioritisation(
                results_dir=Path("/path/to/tool/corpus_results"),
                ranks={},
                rank_stats=RankStats(
                    top=1,
                    top3=2,
                    top5=3,
                    top10=9,
                    found=20,
                    total=30,
                    reciprocal_ranks=[1, 1 / 3, 1 / 5, 1 / 10, 1 / 50],
                ),
            ),
            variant_prioritisation=TrackVariantPrioritisation(
                results_dir=Path("/path/to/tool/corpus_results"),
                ranks={},
                rank_stats=RankStats(
                    top=1,
                    top3=2,
                    top5=3,
                    top10=4,
                    found=5,
                    total=2,
                    reciprocal_ranks=[1, 1 / 3, 1 / 5, 1 / 10, 1 / 12],
                ),
            ),
        )

    def test__retrieve_prioritisation_data_gene(self):
        self.assertEqual(
            self.gene_plot_generator._retrieve_prioritisation_data(self.track_prioritisation),
            TrackGenePrioritisation(
                results_dir=Path("/path/to/tool/corpus_results"),
                ranks={},
                rank_stats=RankStats(
                    top=1,
                    top3=2,
                    top5=3,
                    top10=9,
                    found=20,
                    total=30,
                    reciprocal_ranks=[1, 0.3333333333333333, 0.2, 0.1, 0.02],
                ),
            ),
        )

    def test__retrieve_prioritisation_data_variant(self):
        self.assertEqual(
            self.variant_plot_generator._retrieve_prioritisation_data(self.track_prioritisation),
            TrackVariantPrioritisation(
                results_dir=PosixPath("/path/to/tool/corpus_results"),
                ranks={},
                rank_stats=RankStats(
                    top=1,
                    top3=2,
                    top5=3,
                    top10=4,
                    found=5,
                    total=2,
                    reciprocal_ranks=[1, 0.3333333333333333, 0.2, 0.1, 0.08333333333333333],
                ),
            ),
        )

    def test__generate_stacked_bar_plot_data(self):
        self.gene_plot_generator._generate_stacked_bar_plot_data(self.track_prioritisation)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {
                    "2-3": 5.0,
                    "4-5": 5.0,
                    "6-10": 30.0,
                    ">10": 21.66666666666667,
                    "FO/NP": 33.33333333333333,
                    "Run": "tool_corpus_results",
                    "Top": 5.0,
                }
            ],
        )

    def test__generate_stats_mrr_bar_plot_data(self):
        self.gene_plot_generator._generate_stats_mrr_bar_plot_data(self.track_prioritisation)
        self.assertEqual(
            self.gene_plot_generator.mrr,
            [{"Rank": "MRR", "Percentage": 0.33066666666666666, "Run": "tool_corpus_results"}],
        )

    def test__generate_cumulative_bar_plot_data(self):
        self.gene_plot_generator._generate_cumulative_bar_plot_data(self.track_prioritisation)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {"Rank": "Top", "Percentage": 0.05, "Run": "tool_corpus_results"},
                {"Rank": "Top3", "Percentage": 0.1, "Run": "tool_corpus_results"},
                {"Rank": "Top5", "Percentage": 0.15, "Run": "tool_corpus_results"},
                {"Rank": "Top10", "Percentage": 0.45, "Run": "tool_corpus_results"},
                {"Rank": "Found", "Percentage": 0.6666666666666667, "Run": "tool_corpus_results"},
                {"Rank": "FO/NP", "Percentage": 0.33333333333333326, "Run": "tool_corpus_results"},
                {"Rank": "MRR", "Percentage": 0.33066666666666666, "Run": "tool_corpus_results"},
            ],
        )

    def test__generate_non_cumulative_bar_plot_data(self):
        self.gene_plot_generator._generate_non_cumulative_bar_plot_data(self.track_prioritisation)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {"Rank": "Top", "Percentage": 0.05, "Run": "tool_corpus_results"},
                {"Rank": "2-3", "Percentage": 0.05, "Run": "tool_corpus_results"},
                {"Rank": "4-5", "Percentage": 0.05, "Run": "tool_corpus_results"},
                {"Rank": "6-10", "Percentage": 0.3, "Run": "tool_corpus_results"},
                {"Rank": ">10", "Percentage": 0.2166666666666667, "Run": "tool_corpus_results"},
                {"Rank": "FO/NP", "Percentage": 0.33333333333333326, "Run": "tool_corpus_results"},
                {"Rank": "MRR", "Percentage": 0.33066666666666666, "Run": "tool_corpus_results"},
            ],
        )
