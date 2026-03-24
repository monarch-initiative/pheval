import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import matplotlib.pyplot as plt
import polars as pl
from duckdb import DuckDBPyConnection

from pheval.analyse.benchmark_output_type import BenchmarkOutputTypeEnum
from pheval.analyse.generate_plots import PlotGenerator, generate_plots
from pheval.analyse.run_data_parser import PlotCustomisation, SinglePlotCustomisation


class TestPlotGenerator(unittest.TestCase):
    def setUp(self):
        self.plot_generator = PlotGenerator("test_benchmark", Path("out"))

        self.benchmarking_results_df = pl.DataFrame(
            {
                "run_identifier": ["Run1", "Run2"],
                "percentage@1": [50, 60],
                "percentage@3": [70, 75],
                "percentage@5": [80, 85],
                "percentage@10": [90, 92],
                "percentage_found": [95, 98],
                "mrr": [0.8, 0.85],
            }
        )

        self.curves_df = pl.DataFrame(
            {
                "run_identifier": ["Run1"],
                "fpr": [[0.0, 0.1, 0.2, 1.0]],
                "tpr": [[0.0, 0.4, 0.8, 1.0]],
                "precision": [[0.8, 0.75, 0.6, 0.5]],
                "recall": [[1.0, 0.9, 0.7, 0.0]],
            }
        )

        self.mock_single_plot_customisation = MagicMock(
            SinglePlotCustomisation(
                plot_type="bar_stacked",
                rank_plot_title="Test Title",
                roc_curve_title="Test ROC Title",
                precision_recall_title="Test PR Title",
            )
        )
        self.multi_plot_customisation = PlotCustomisation(
            gene_plots=SinglePlotCustomisation(
                plot_type="bar_stacked",
                rank_plot_title="Test Title",
                roc_curve_title="Test ROC Title",
                precision_recall_title="Test PR Title",
            ),
            variant_plots=SinglePlotCustomisation(
                plot_type="bar_stacked",
                rank_plot_title="Test Title",
                roc_curve_title="Test ROC Title",
                precision_recall_title="Test PR Title",
            ),
            disease_plots=SinglePlotCustomisation(
                plot_type="bar_stacked",
                rank_plot_title="Test Title",
                roc_curve_title="Test ROC Title",
                precision_recall_title="Test PR Title",
            ),
        )
        self.mock_benchmark_output_type = MagicMock(BenchmarkOutputTypeEnum.GENE.value)
        self.mock_conn = MagicMock(spec=DuckDBPyConnection)

    def test__generate_stacked_data(self):
        self.assertTrue(
            self.plot_generator._generate_stacked_data(self.benchmarking_results_df).equals(
                pl.DataFrame(
                    [
                        {
                            "Run": "Run1",
                            "Top": 50,
                            "2-3": 20,
                            "4-5": 10,
                            "6-10": 10,
                            ">10": 5,
                            "Missed": 5,
                        },
                        {
                            "Run": "Run2",
                            "Top": 60,
                            "2-3": 15,
                            "4-5": 10,
                            "6-10": 7,
                            ">10": 6,
                            "Missed": 2,
                        },
                    ]
                )
            )
        )

    def test__generate_cumulative_bar_plot_data(self):
        self.assertTrue(
            self.plot_generator._generate_cumulative_bar_plot_data(self.benchmarking_results_df).equals(
                pl.DataFrame(
                    [
                        {
                            "Run": "Run1",
                            "Top": 0.5,
                            "Top3": 0.7000000000000001,
                            "Top5": 0.8,
                            "Top10": 0.9,
                            "Found": 0.9500000000000001,
                            "MRR": 0.8,
                        },
                        {
                            "Run": "Run2",
                            "Top": 0.6,
                            "Top3": 0.75,
                            "Top5": 0.85,
                            "Top10": 0.92,
                            "Found": 0.98,
                            "MRR": 0.85,
                        },
                    ]
                )
            )
        )

    def test__generate_non_cumulative_bar_plot_data(self):
        self.assertTrue(
            self.plot_generator._generate_non_cumulative_bar_plot_data(self.benchmarking_results_df).equals(
                pl.DataFrame(
                    [
                        {
                            "Run": "Run1",
                            "Top": 50,
                            "2-3": 20,
                            "4-5": 10,
                            "6-10": 10,
                            ">10": 5,
                            "Missed": 5,
                            "MRR": 0.8,
                        },
                        {
                            "Run": "Run2",
                            "Top": 60,
                            "2-3": 15,
                            "4-5": 10,
                            "6-10": 7,
                            ">10": 6,
                            "Missed": 2,
                            "MRR": 0.85,
                        },
                    ]
                )
            )
        )

    def test__extract_mrr_data(self):
        self.assertTrue(
            self.plot_generator._extract_mrr_data(self.benchmarking_results_df).equals(
                pl.DataFrame([{"Run": "Run1", "Percentage": 0.8}, {"Run": "Run2", "Percentage": 0.85}])
            )
        )

    @patch.object(plt, "savefig")
    def test_generate_stacked_bar_plot(self, mock_savefig):
        self.plot_generator.generate_stacked_bar_plot(
            self.benchmarking_results_df,
            self.mock_benchmark_output_type,
            self.mock_single_plot_customisation,
        )
        mock_savefig.asser_called_with()

    @patch.object(plt, "savefig")
    def test_generate_roc_curve(self, mock_savefig):
        self.plot_generator.generate_roc_curve(
            self.curves_df, self.mock_benchmark_output_type, self.mock_single_plot_customisation
        )
        mock_savefig.assert_called_once()

    @patch.object(plt, "savefig")
    def test_generate_precision_recall(self, mock_savefig):
        self.plot_generator.generate_precision_recall(
            self.curves_df, self.mock_benchmark_output_type, self.mock_single_plot_customisation
        )
        mock_savefig.assert_called_once()

    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_roc_curve")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_precision_recall")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_stacked_bar_plot")
    def test_generate_plots_single_run_no_rank_changes(self, mock_stacked_plot, mock_pr_curve, mock_roc_curve):
        """With a single run, ROC/PR/stacked plots are generated but no rank change plots."""
        generate_plots(
            "test_benchmark",
            self.benchmarking_results_df,
            self.curves_df,
            BenchmarkOutputTypeEnum.GENE.value,
            self.multi_plot_customisation,
            Path("out"),
            False,
            self.mock_conn,
            ["Run1"],
        )
        mock_roc_curve.assert_called_once()
        mock_pr_curve.assert_called_once()
        mock_stacked_plot.assert_called_once()

    @patch("pheval.analyse.generate_plots.load_table_lazy")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_rank_change_plot")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_roc_curve")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_precision_recall")
    @patch("pheval.analyse.generate_plots.PlotGenerator.generate_stacked_bar_plot")
    def test_generate_plots_two_runs_triggers_rank_change_plot(
        self, mock_stacked_plot, mock_pr_curve, mock_roc_curve, mock_rank_change_plot, mock_load_table
    ):
        """With two runs, rank change plot is generated in addition to the standard plots."""
        mock_load_table.return_value.collect.return_value = pl.DataFrame(
            {"run1": [1], "run2": [2], "rank_change": ["1"]}
        )
        generate_plots(
            "test_benchmark",
            self.benchmarking_results_df,
            self.curves_df,
            BenchmarkOutputTypeEnum.GENE.value,
            self.multi_plot_customisation,
            Path("out"),
            False,
            self.mock_conn,
            ["Run1", "Run2"],
        )
        mock_roc_curve.assert_called_once()
        mock_pr_curve.assert_called_once()
        mock_stacked_plot.assert_called_once()
        mock_rank_change_plot.assert_called_once()


class TestClassifyRankChanges(unittest.TestCase):
    def setUp(self):
        self.plot_generator = PlotGenerator("test_benchmark", Path("out"))
        self.rank_changes_df = pl.DataFrame(
            {
                "result_file": ["file1", "file2", "file3", "file4", "file5"],
                "run1": [0, 2, 3, 1, 5],
                "run2": [1, 0, 1, 2, 5],
                "rank_change": ["GAINED", "LOST", "2", "-1", "0"],
            }
        )

    def test_gained_is_improved_with_run2_as_reference(self):
        result = self.plot_generator._classify_rank_changes(self.rank_changes_df, "run1", "run2")
        row = result.filter(pl.col("result_file") == "file1")
        self.assertEqual(row["outcome"][0], "Improved")
        self.assertEqual(row["reference_rank"][0], 1)  # run2 rank

    def test_lost_is_dropped_with_run1_as_reference(self):
        result = self.plot_generator._classify_rank_changes(self.rank_changes_df, "run1", "run2")
        row = result.filter(pl.col("result_file") == "file2")
        self.assertEqual(row["outcome"][0], "Dropped")
        self.assertEqual(row["reference_rank"][0], 2)  # run1 rank

    def test_positive_delta_is_improved(self):
        result = self.plot_generator._classify_rank_changes(self.rank_changes_df, "run1", "run2")
        row = result.filter(pl.col("result_file") == "file3")
        self.assertEqual(row["outcome"][0], "Improved")
        self.assertEqual(row["reference_rank"][0], 3)  # run1 rank

    def test_negative_delta_is_dropped(self):
        result = self.plot_generator._classify_rank_changes(self.rank_changes_df, "run1", "run2")
        row = result.filter(pl.col("result_file") == "file4")
        self.assertEqual(row["outcome"][0], "Dropped")
        self.assertEqual(row["reference_rank"][0], 1)  # run1 rank

    def test_zero_delta_is_unchanged(self):
        result = self.plot_generator._classify_rank_changes(self.rank_changes_df, "run1", "run2")
        row = result.filter(pl.col("result_file") == "file5")
        self.assertEqual(row["outcome"][0], "Unchanged")
        self.assertEqual(row["reference_rank"][0], 5)  # run1 rank


class TestComputeRankChangePlotData(unittest.TestCase):
    def setUp(self):
        self.plot_generator = PlotGenerator("test_benchmark", Path("out"))
        self.classified_df = pl.DataFrame(
            {
                "reference_rank": [1, 1, 1, 2, 15, 15, 250],
                "outcome": ["Improved", "Improved", "Dropped", "Unchanged", "Improved", "Dropped", "Improved"],
            }
        )

    def test_cases_beyond_200_are_excluded(self):
        result = self.plot_generator._compute_rank_change_plot_data(self.classified_df)
        self.assertNotIn("250", result["rank_bin"].cast(pl.Utf8).to_list())

    def test_bin_counts_are_correct(self):
        result = self.plot_generator._compute_rank_change_plot_data(self.classified_df)
        n_by_bin = {row["rank_bin"]: row["n"] for row in result.iter_rows(named=True)}
        self.assertEqual(n_by_bin["1"], 3)
        self.assertEqual(n_by_bin["2"], 1)
        self.assertEqual(n_by_bin["11-20"], 2)

    def test_proportions_sum_to_one_per_bin(self):
        result = self.plot_generator._compute_rank_change_plot_data(self.classified_df)
        outcome_cols = [c for c in result.columns if c in ("Improved", "Unchanged", "Dropped")]
        for row in result.iter_rows(named=True):
            total = sum(row[col] for col in outcome_cols)
            self.assertAlmostEqual(total, 1.0)

    def test_bins_are_ordered(self):
        result = self.plot_generator._compute_rank_change_plot_data(self.classified_df)
        bins = result["rank_bin"].cast(pl.Utf8).to_list()
        expected_order = [b for b in self.plot_generator.rank_bins if b in bins]
        self.assertEqual(bins, expected_order)

    def test_rank_1_proportions(self):
        result = self.plot_generator._compute_rank_change_plot_data(self.classified_df)
        row = result.filter(pl.col("rank_bin").cast(pl.Utf8) == "1").to_dicts()[0]
        self.assertAlmostEqual(row["Improved"], 2 / 3)
        self.assertAlmostEqual(row["Dropped"], 1 / 3)


class TestGenerateRankChangePlot(unittest.TestCase):
    def setUp(self):
        self.plot_generator = PlotGenerator("test_benchmark", Path("out"))
        self.mock_benchmark_output_type = MagicMock(BenchmarkOutputTypeEnum.GENE.value)
        self.mock_benchmark_output_type.prioritisation_type_string = "gene"
        self.rank_changes_df = pl.DataFrame(
            {
                "result_file": ["file1", "file2", "file3"],
                "run1": [0, 2, 3],
                "run2": [1, 0, 1],
                "rank_change": ["GAINED", "LOST", "2"],
            }
        )

    @patch("pheval.analyse.generate_plots.plt.close")
    @patch("matplotlib.figure.Figure.savefig")
    def test_generate_rank_change_plot_saves_file(self, mock_savefig, mock_close):
        self.plot_generator.generate_rank_change_plot(
            self.rank_changes_df, "run1", "run2", self.mock_benchmark_output_type
        )
        mock_savefig.assert_called_once()
        mock_close.assert_called_once()

    @patch("pheval.analyse.generate_plots.logger")
    def test_generate_rank_change_plot_warns_on_empty_data(self, mock_logger):
        empty_df = pl.DataFrame(
            {"result_file": [], "run1": [], "run2": [], "rank_change": []},
            schema={"result_file": pl.Utf8, "run1": pl.Int64, "run2": pl.Int64, "rank_change": pl.Utf8},
        )
        self.plot_generator.generate_rank_change_plot(
            empty_df, "run1", "run2", self.mock_benchmark_output_type
        )
        mock_logger.warning.assert_called_once()
