import shutil
import tempfile
import unittest
from pathlib import Path

from pheval.analyse.generate_plots import BenchmarkRunResults, PlotGenerator
from pheval.analyse.rank_stats import RankStats


class TestPlotGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.gene_plot_generator = PlotGenerator()
        self.variant_plot_generator = PlotGenerator()
        self.disease_plot_generator = PlotGenerator()
        self.benchmarking_result = BenchmarkRunResults(
            results_dir="tool_corpus",
            ranks={},
            rank_stats=RankStats(
                top=1,
                top3=2,
                top5=4,
                top10=8,
                found=10,
                total=20,
                reciprocal_ranks=[1, 1 / 3, 1 / 5, 1 / 10, 1 / 50],
            ),
        )
        self.test_dir = tempfile.mkdtemp()
        self.results_dir = Path(self.test_dir).joinpath("tool1-default/corpus1-default")
        self.results_dir.mkdir(exist_ok=True, parents=True)

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_create_run_identifier(self):
        self.assertEqual(
            self.gene_plot_generator.create_run_identifier(self.results_dir),
            "tool1-default_corpus1-default",
        )

    def test__generate_stacked_bar_plot_data(self):
        self.gene_plot_generator._generate_stacked_bar_plot_data(self.benchmarking_result)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {
                    "2-3": 5.0,
                    "4-5": 10.0,
                    "6-10": 20.0,
                    ">10": 10.0,
                    "Missed": 50.0,
                    "Run": "tool_corpus",
                    "Top": 5.0,
                }
            ],
        )

    def test__generate_stats_mrr_bar_plot_data(self):
        self.gene_plot_generator._generate_stats_mrr_bar_plot_data(self.benchmarking_result)
        self.assertEqual(
            self.gene_plot_generator.mrr,
            [{"Rank": "MRR", "Percentage": 0.11022222222222222, "Run": "tool_corpus"}],
        )

    def test__generate_cumulative_bar_plot_data(self):
        self.gene_plot_generator._generate_cumulative_bar_plot_data(self.benchmarking_result)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {"Percentage": 0.05, "Rank": "Top", "Run": "tool_corpus"},
                {"Percentage": 0.1, "Rank": "Top3", "Run": "tool_corpus"},
                {"Percentage": 0.2, "Rank": "Top5", "Run": "tool_corpus"},
                {"Percentage": 0.4, "Rank": "Top10", "Run": "tool_corpus"},
                {"Percentage": 0.5, "Rank": "Found", "Run": "tool_corpus"},
                {"Percentage": 0.5, "Rank": "Missed", "Run": "tool_corpus"},
                {"Percentage": 0.11022222222222222, "Rank": "MRR", "Run": "tool_corpus"},
            ],
        )

    def test__generate_non_cumulative_bar_plot_data(self):
        self.gene_plot_generator._generate_non_cumulative_bar_plot_data(self.benchmarking_result)
        self.assertEqual(
            self.gene_plot_generator.stats,
            [
                {"Percentage": 0.05, "Rank": "Top", "Run": "tool_corpus"},
                {"Percentage": 0.05, "Rank": "2-3", "Run": "tool_corpus"},
                {"Percentage": 0.1, "Rank": "4-5", "Run": "tool_corpus"},
                {"Percentage": 0.2, "Rank": "6-10", "Run": "tool_corpus"},
                {"Percentage": 0.1, "Rank": ">10", "Run": "tool_corpus"},
                {"Percentage": 0.5, "Rank": "Missed", "Run": "tool_corpus"},
                {"Percentage": 0.11022222222222222, "Rank": "MRR", "Run": "tool_corpus"},
            ],
        )
