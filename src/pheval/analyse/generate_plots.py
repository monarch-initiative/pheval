from enum import Enum
from pathlib import Path

import duckdb
import matplotlib
import polars as pl
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.metrics import auc

from pheval.analyse.benchmark_db_manager import load_table_lazy
from pheval.analyse.benchmark_output_type import (
    BenchmarkOutputType,
    BenchmarkOutputTypeEnum,
)
from pheval.analyse.run_data_parser import (
    PlotCustomisation,
    SinglePlotCustomisation,
    parse_run_config,
)
from pheval.utils.logger import get_logger

logger = get_logger()


class PlotTypes(Enum):
    BAR_STACKED = "bar_stacked"
    BAR_CUMULATIVE = "bar_cumulative"
    BAR_NON_CUMULATIVE = "bar_non_cumulative"


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

    def __init__(self, benchmark_name: str):
        """
        Initialise the PlotGenerator class.
        Note:
            Matplotlib settings are configured to remove the right and top axes spines
            for generated plots.
        """
        self.benchmark_name = benchmark_name
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False

    @staticmethod
    def _generate_stacked_data(benchmarking_stats_df: pl.DataFrame) -> pl.DataFrame:
        """
        Generate stacked data.
        Args:
            benchmarking_stats_df (pl.DataFrame): benchmarking stats dataframe.
        Returns:
            pl.DataFrame: Data formatted for plotting stacked data.
        """
        return benchmarking_stats_df.with_columns(
            [
                pl.col("run_identifier").alias("Run"),
                pl.col("percentage@1").alias("Top"),
                (pl.col("percentage@3") - pl.col("percentage@1")).alias("2-3"),
                (pl.col("percentage@5") - pl.col("percentage@3")).alias("4-5"),
                (pl.col("percentage@10") - pl.col("percentage@5")).alias("6-10"),
                (pl.col("percentage_found") - pl.col("percentage@10")).alias(">10"),
                (100 - pl.col("percentage_found")).alias("Missed"),
            ]
        ).select(["Run", "Top", "2-3", "4-5", "6-10", ">10", "Missed"])

    @staticmethod
    def _extract_mrr_data(benchmarking_results_df: pl.DataFrame) -> pl.DataFrame:
        """
        Generate data in the correct format for dataframe creation for MRR (Mean Reciprocal Rank) bar plot.

        Args:
            benchmarking_results_df (pl.DataFrame): benchmarking stats dataframe.
        Returns:
            pl.DataFrame: Data formatted for plotting MRR bar plot.
        """
        return benchmarking_results_df.select(["run_identifier", "mrr"]).rename(
            {"run_identifier": "Run", "mrr": "Percentage"}
        )

    def _save_fig(
        self, benchmark_output_type: BenchmarkOutputType, y_lower_limit: int, y_upper_limit: int
    ) -> None:
        """
        Save the generated figure.
        Args:
            benchmark_output_type (BenchmarkOutputType): Benchmark output type.
            y_lower_limit (int): Lower limit for the y-axis.
            y_upper_limit (int): Upper limit for the y-axis.
        """
        plt.ylim(y_lower_limit, y_upper_limit)
        plt.savefig(
            f"{self.benchmark_name}_{benchmark_output_type.prioritisation_type_string}_rank_stats.svg",
            format="svg",
            bbox_inches="tight",
        )

    def generate_stacked_bar_plot(
        self,
        benchmarking_results_df: pl.DataFrame,
        benchmark_output_type: BenchmarkOutputType,
        plot_customisation: SinglePlotCustomisation,
    ) -> None:
        """
        Generate a stacked bar plot and Mean Reciprocal Rank (MRR) bar plot.
        Args:
            benchmarking_results_df (pl.DataFrame): benchmarking stats dataframe.
            benchmark_output_type (BenchmarkOutputType): Benchmark output type.
            plot_customisation (SinglePlotCustomisation): Plotting customisation.
        """
        plt.clf()
        stats_df = self._generate_stacked_data(benchmarking_results_df)
        stats_df.to_pandas().set_index("Run").plot(
            kind="bar",
            stacked=True,
            color=self.palette_hex_codes,
            ylabel=benchmark_output_type.y_label,
            edgecolor="white",
        ).legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        plt.title(plot_customisation.rank_plot_title, loc="center", fontsize=15)
        self._save_fig(benchmark_output_type, 0, 100)
        mrr_df = self._extract_mrr_data(benchmarking_results_df)
        mrr_df.to_pandas().set_index("Run").plot(
            kind="bar",
            color=self.palette_hex_codes,
            ylabel=f"{benchmark_output_type.prioritisation_type_string.capitalize()} mean reciprocal rank",
            legend=False,
            edgecolor="white",
        )
        plt.title(
            f"{benchmark_output_type.prioritisation_type_string.capitalize()} results - mean reciprocal rank"
        )
        self._save_fig(benchmark_output_type, 0, 1)

    @staticmethod
    def _generate_cumulative_bar_plot_data(benchmarking_results_df: pl.DataFrame) -> pl.DataFrame:
        """
        Generate data in the correct format for dataframe creation for a cumulative bar plot,
        appending to the self.stats attribute of the class.
        """
        return benchmarking_results_df.select(
            [
                pl.col("run_identifier").alias("Run"),
                pl.col("percentage@1").alias("Top") / 100,
                pl.col("percentage@3").alias("Top3") / 100,
                pl.col("percentage@5").alias("Top5") / 100,
                pl.col("percentage@10").alias("Top10") / 100,
                pl.col("percentage_found").alias("Found") / 100,
                pl.col("mrr").alias("MRR"),
            ]
        )

    def _plot_bar_plot(
        self,
        benchmark_output_type: BenchmarkOutputType,
        stats_df: pl.DataFrame,
        plot_customisation: SinglePlotCustomisation,
    ) -> None:
        stats_df = stats_df.to_pandas().melt(
            id_vars=["Run"],
            value_vars=["Top", "Top3", "Top5", "Top10", "Found", "MRR"],
            var_name="Rank",
            value_name="Percentage",
        )
        sns.catplot(
            data=stats_df,
            kind="bar",
            x="Rank",
            y="Percentage",
            hue="Run",
            palette=self.palette_hex_codes,
            edgecolor="white",
            legend=False,
        ).set(xlabel="Rank", ylabel=benchmark_output_type.y_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, title="Run")
        plt.title(plot_customisation.rank_plot_title, loc="center", fontsize=15)
        self._save_fig(benchmark_output_type, 0, 1)

    def _generate_non_cumulative_bar_plot_data(
        self, benchmarking_results_df: pl.DataFrame
    ) -> pl.DataFrame:
        """
        Generate data in the correct format for dataframe creation for a non-cumulative bar plot,
        appending to the self.stats attribute of the class.
        """
        return self._generate_stacked_data(benchmarking_results_df).hstack(
            self._extract_mrr_data(benchmarking_results_df).select(
                pl.col("Percentage").alias("MRR")
            )
        )

    def generate_cumulative_bar(
        self,
        benchmarking_results_df: pl.DataFrame,
        benchmark_generator: BenchmarkOutputType,
        plot_customisation: SinglePlotCustomisation,
    ) -> None:
        """
        Generate a cumulative bar plot.
        """
        plt.clf()
        stats_df = self._generate_cumulative_bar_plot_data(benchmarking_results_df)
        self._plot_bar_plot(benchmark_generator, stats_df, plot_customisation)

    def generate_non_cumulative_bar(
        self,
        benchmarking_results_df: pl.DataFrame,
        benchmark_generator: BenchmarkOutputType,
        plot_customisation: SinglePlotCustomisation,
    ) -> None:
        """
        Generate a non-cumulative bar plot.
        """
        plt.clf()
        stats_df = self._generate_non_cumulative_bar_plot_data(benchmarking_results_df)
        self._plot_bar_plot(benchmark_generator, stats_df, plot_customisation)

    def generate_roc_curve(
        self,
        curves: pl.DataFrame,
        benchmark_generator: BenchmarkOutputType,
        plot_customisation: SinglePlotCustomisation,
    ):
        """
        Generate and plot Receiver Operating Characteristic (ROC) curves for binary classification benchmark results.

        Args:
        """
        plt.clf()
        for i, row in enumerate(curves.iter_rows(named=True)):
            run_identifier = row["run_identifier"]
            fpr = row["fpr"]
            tpr = row["tpr"]
            roc_auc = auc(fpr, tpr)
            plt.plot(
                fpr,
                tpr,
                label=f"{run_identifier} ROC Curve (AUC = {roc_auc:.2f})",
                color=self.palette_hex_codes[i],
            )
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(plot_customisation.roc_curve_title)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            f"{self.benchmark_name}_{benchmark_generator.prioritisation_type_string}_roc_curve.svg",
            format="svg",
            bbox_inches="tight",
        )

    def generate_precision_recall(
        self,
        curves: pl.DataFrame,
        benchmark_generator: BenchmarkOutputType,
        plot_customisation: SinglePlotCustomisation,
    ):
        """
        Generate and plot Precision-Recall curves for binary classification benchmark results.
        """
        plt.clf()
        plt.figure()
        for i, row in enumerate(curves.iter_rows(named=True)):
            run_identifier = row["run_identifier"]
            precision = row["precision"]
            recall = row["recall"]
            pr_auc = auc(recall[::-1], precision[::-1])
            plt.plot(
                recall,
                precision,
                label=f"{run_identifier} Precision-Recall Curve (AUC = {pr_auc:.2f})",
                color=self.palette_hex_codes[i],
            )
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(plot_customisation.precision_recall_title)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            f"{self.benchmark_name}_{benchmark_generator.prioritisation_type_string}_pr_curve.svg",
            format="svg",
            bbox_inches="tight",
        )


def generate_plots(
    benchmark_name: str,
    benchmarking_results_df: pl.DataFrame,
    curves: pl.DataFrame,
    benchmark_output_type: BenchmarkOutputType,
    plot_customisation: PlotCustomisation,
) -> None:
    """
    Generate summary statistics bar plots for prioritisation.

    This method generates summary statistics bar plots based on the provided benchmarking results and plot type.
    """
    plot_generator = PlotGenerator(benchmark_name)
    plot_customisation_type = getattr(
        plot_customisation, f"{benchmark_output_type.prioritisation_type_string}_plots"
    )
    logger.info("Generating ROC curve visualisations.")
    plot_generator.generate_roc_curve(curves, benchmark_output_type, plot_customisation_type)
    logger.info("Generating Precision-Recall curves visualisations.")
    plot_generator.generate_precision_recall(curves, benchmark_output_type, plot_customisation_type)
    plot_type = PlotTypes(plot_customisation_type.plot_type)
    match plot_type:
        case PlotTypes.BAR_STACKED:
            logger.info("Generating stacked bar plot.")
            plot_generator.generate_stacked_bar_plot(
                benchmarking_results_df, benchmark_output_type, plot_customisation_type
            )
        case PlotTypes.BAR_CUMULATIVE:
            logger.info("Generating cumulative bar plot.")
            plot_generator.generate_cumulative_bar(
                benchmarking_results_df, benchmark_output_type, plot_customisation_type
            )
        case PlotTypes.BAR_NON_CUMULATIVE:
            logger.info("Generating non cumulative bar plot.")
            plot_generator.generate_non_cumulative_bar(
                benchmarking_results_df, benchmark_output_type, plot_customisation_type
            )


def generate_plots_from_db(db_path: Path, config: Path) -> None:
    """
    Generate plots from database file.
    Args:
        db_path (Path): Path to the database file.
        config (Path): Path to the benchmarking config file.
    """
    logger.info(f"Generating plots from {db_path}")
    conn = duckdb.connect(db_path)
    logger.info(f"Parsing configurations from {config}")
    benchmark_config_file = parse_run_config(config)
    tables = {
        row[0]
        for row in conn.execute(
            """SELECT table_name FROM duckdb_tables WHERE table_name """
            """LIKE '%_summary%' OR table_name LIKE '%_binary_classification_curves'"""
        ).fetchall()
    }
    for benchmark_output_type in BenchmarkOutputTypeEnum:
        summary_table = (
            f"{benchmark_config_file.benchmark_name}_"
            f"{benchmark_output_type.value.prioritisation_type_string}_summary"
        )
        curve_table = (
            f"{benchmark_config_file.benchmark_name}_"
            f"{benchmark_output_type.value.prioritisation_type_string}_binary_classification_curves"
        )
        if summary_table in tables and curve_table in tables:
            logger.info(
                f"Generating plots for {benchmark_output_type.value.prioritisation_type_string} prioritisation."
            )
            benchmarking_results_df = load_table_lazy(summary_table, conn).collect()
            curves_df = load_table_lazy(curve_table, conn).collect()
            generate_plots(
                benchmark_name=benchmark_config_file.benchmark_name,
                benchmarking_results_df=benchmarking_results_df,
                curves=curves_df,
                benchmark_output_type=benchmark_output_type.value,
                plot_customisation=benchmark_config_file.plot_customisation,
            )
    logger.info("Finished generating plots.")
    conn.close()
