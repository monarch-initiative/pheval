import unittest
from unittest.mock import MagicMock, patch

import matplotlib.pyplot as plt
import polars as pl

from pheval.analyse.benchmark_output_type import BenchmarkOutputTypeEnum
from pheval.analyse.generate_plots import PlotGenerator, generate_plots
from pheval.analyse.run_data_parser import PlotCustomisation, SinglePlotCustomisation


class TestPlotGenerator(unittest.TestCase):
    def setUp(self):
        self.plot_generator = PlotGenerator("test_benchmark")

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
    def test_generate_plots(self, mock_stacked_plot, mock_pr_curve, mock_roc_curve):
        """Test the full plot generation pipeline."""
        generate_plots(
            "test_benchmark",
            self.benchmarking_results_df,
            self.curves_df,
            BenchmarkOutputTypeEnum.GENE.value,
            self.multi_plot_customisation,
        )

        mock_roc_curve.assert_called_once()
        mock_pr_curve.assert_called_once()
        mock_stacked_plot.assert_called_once()
