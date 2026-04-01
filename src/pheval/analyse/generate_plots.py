# ruff: noqa: PLR2004
from enum import Enum
from itertools import combinations
from pathlib import Path
from typing import ClassVar

import duckdb
import matplotlib
import matplotlib.colors as mcolors
import polars as pl
import seaborn as sns
from duckdb import DuckDBPyConnection
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

    palette_hex_codes: ClassVar[list[str]] = [
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

    rank_bins: ClassVar[list[str]] = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11-20",
        "21-50",
        "51-100",
        "101-200",
    ]

    def __init__(self, benchmark_name: str, output_dir: Path):
        """
        Initialise the PlotGenerator class.
        Note:
            Matplotlib settings are configured to remove the right and top axes spines
            for generated plots.
        """
        self.benchmark_name = benchmark_name
        self.output_dir = output_dir
        matplotlib.rcParams["axes.spines.right"] = False
        matplotlib.rcParams["axes.spines.top"] = False

    def get_palette(self, n_colors: int) -> list[str]:
        """
        Generates a palette of colour hex codes with the specified number of colors.

        Args:
            n_colors (int): The number of colour hex codes to generate.

        Returns:
            list[str]: A list containing the generated colour hex codes.
        """
        if n_colors <= len(self.palette_hex_codes):
            return self.palette_hex_codes
        cmap = mcolors.LinearSegmentedColormap.from_list("pheval", self.palette_hex_codes, N=n_colors)
        return [mcolors.rgb2hex(cmap(i)) for i in range(n_colors)]

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

    def _save_fig(self, benchmark_output_type: BenchmarkOutputType, y_lower_limit: int, y_upper_limit: int) -> None:
        """
        Save the generated figure.
        Args:
            benchmark_output_type (BenchmarkOutputType): Benchmark output type.
            y_lower_limit (int): Lower limit for the y-axis.
            y_upper_limit (int): Upper limit for the y-axis.
        """
        plt.ylim(y_lower_limit, y_upper_limit)
        plt.savefig(
            self.output_dir.joinpath(
                f"{self.benchmark_name}_{benchmark_output_type.prioritisation_type_string}_rank_stats.svg"
            ),
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
            color=self.get_palette(len(mrr_df)),
            ylabel=f"{benchmark_output_type.prioritisation_type_string.capitalize()} mean reciprocal rank",
            legend=False,
            edgecolor="white",
        )
        plt.title(f"{benchmark_output_type.prioritisation_type_string.capitalize()} results - mean reciprocal rank")
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
            palette=self.get_palette(stats_df["Run"].nunique()),
            edgecolor="white",
            legend=False,
        ).set(xlabel="Rank", ylabel=benchmark_output_type.y_label)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, title="Run")
        plt.title(plot_customisation.rank_plot_title, loc="center", fontsize=15)
        self._save_fig(benchmark_output_type, 0, 1)

    def _generate_non_cumulative_bar_plot_data(self, benchmarking_results_df: pl.DataFrame) -> pl.DataFrame:
        """
        Generate data in the correct format for dataframe creation for a non-cumulative bar plot,
        appending to the self.stats attribute of the class.
        """
        return self._generate_stacked_data(benchmarking_results_df).hstack(
            self._extract_mrr_data(benchmarking_results_df).select(pl.col("Percentage").alias("MRR"))
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

    @staticmethod
    def _classify_rank_changes(rank_changes_df: pl.DataFrame, run1: str, run2: str) -> pl.DataFrame:
        """
        Classify each case as Improved, Unchanged, or Dropped and assign a reference rank.

        GAINED (rank 0 → ranked): Improved, reference rank taken from run2.
        LOST (ranked → rank 0): Dropped, reference rank taken from run1.
        Numeric delta > 0: Improved (run2 rank is better).
        Numeric delta == 0: Unchanged.
        Numeric delta < 0: Dropped.

        Args:
            rank_changes_df: DataFrame with run1/run2 rank columns and a rank_change column.
            run1: Baseline run identifier.
            run2: Comparison run identifier.

        Returns:
            pl.DataFrame: Input rows with added reference_rank and outcome columns.
        """
        rank_change_val = pl.col("rank_change").cast(pl.Int64, strict=False)
        return rank_changes_df.with_columns(
            [
                pl.when(pl.col("rank_change") == "GAINED")
                .then(pl.col(run2))
                .otherwise(pl.col(run1))
                .alias("reference_rank"),
                pl.when(pl.col("rank_change") == "GAINED")
                .then(pl.lit("Improved"))
                .when(pl.col("rank_change") == "LOST")
                .then(pl.lit("Dropped"))
                .when(rank_change_val > 0)
                .then(pl.lit("Improved"))
                .when(rank_change_val == 0)
                .then(pl.lit("Unchanged"))
                .otherwise(pl.lit("Dropped"))
                .alias("outcome"),
            ]
        )

    def _compute_rank_change_plot_data(self, classified_rank_changes_df: pl.DataFrame) -> pl.DataFrame:
        """
        Bin cases by reference rank and compute per-bin outcome proportions.

        Ranks 1-10 are kept individually; higher ranks are grouped into
        11-20, 21-50, 51-100, and 101-200. Cases outside 1-200 are excluded.

        Args:
            classified_rank_changes_df: DataFrame with reference_rank and outcome columns.

        Returns:
            pl.DataFrame: Pivoted DataFrame with one row per rank bin containing
                outcome proportions and total case counts (n), ordered by rank.
        """
        ref = pl.col("reference_rank")
        bin_expr = (
            pl.when(ref == 1)
            .then(pl.lit("1"))
            .when(ref == 2)
            .then(pl.lit("2"))
            .when(ref == 3)
            .then(pl.lit("3"))
            .when(ref == 4)
            .then(pl.lit("4"))
            .when(ref == 5)
            .then(pl.lit("5"))
            .when(ref == 6)
            .then(pl.lit("6"))
            .when(ref == 7)
            .then(pl.lit("7"))
            .when(ref == 8)
            .then(pl.lit("8"))
            .when(ref == 9)
            .then(pl.lit("9"))
            .when(ref == 10)
            .then(pl.lit("10"))
            .when((ref >= 11) & (ref <= 20))
            .then(pl.lit("11-20"))
            .when((ref >= 21) & (ref <= 50))
            .then(pl.lit("21-50"))
            .when((ref >= 51) & (ref <= 100))
            .then(pl.lit("51-100"))
            .when((ref >= 101) & (ref <= 200))
            .then(pl.lit("101-200"))
            .otherwise(None)
        )
        df = classified_rank_changes_df.with_columns(bin_expr.alias("rank_bin")).filter(
            pl.col("rank_bin").is_not_null()
        )
        totals = df.group_by("rank_bin").len().rename({"len": "n"})
        props = (
            df.group_by(["rank_bin", "outcome"])
            .len()
            .join(totals, on="rank_bin")
            .with_columns((pl.col("len") / pl.col("n")).alias("proportion"))
        )
        pivot_df = (
            props.pivot(index="rank_bin", on="outcome", values="proportion", aggregate_function="sum")
            .fill_null(0)
            .join(totals, on="rank_bin")
        )
        ordered_bins = [b for b in self.rank_bins if b in pivot_df["rank_bin"].to_list()]
        return pivot_df.with_columns(pl.col("rank_bin").cast(pl.Enum(ordered_bins))).sort("rank_bin")

    def _plot_rank_change_bar(
        self,
        plot_data: pl.DataFrame,
        run1: str,
        run2: str,
        benchmark_output_type: BenchmarkOutputType,
    ) -> None:
        """
        Draw and save the horizontal stacked bar chart for rank changes.

        Args:
            plot_data: Pivoted DataFrame from _compute_rank_change_plot_data.
            run1: Baseline run identifier (used in title and filename).
            run2: Comparison run identifier (used in title and filename).
            benchmark_output_type: Determines the prioritisation type label and filename.
        """
        pdf = plot_data.to_pandas().set_index("rank_bin")
        outcomes = ["Improved", "Unchanged", "Dropped"]
        colors = {"Improved": "#174857", "Unchanged": "#a6c4d4", "Dropped": "#2b7288"}

        total_n = int(pdf["n"].sum())
        total_counts = {
            outcome: int(round((pdf[outcome] * pdf["n"]).sum()))
            for outcome in outcomes
            if outcome in pdf.columns
        }
        total_row = {outcome: count / total_n for outcome, count in total_counts.items()}

        # Total bar sits one row below the rank bins with a gap
        total_y = len(pdf) + 0.5
        fig, ax = plt.subplots(figsize=(8, max(4, (total_y + 1) * 0.55)))
        lefts = [0.0] * len(pdf)
        total_left = 0.0
        for outcome in outcomes:
            if outcome not in pdf.columns:
                continue
            values = pdf[outcome].values
            ax.barh(range(len(pdf)), values, left=lefts, color=colors[outcome], label=outcome, edgecolor="white")
            for i, (left, proportion) in enumerate(zip(lefts, values, strict=True)):
                count = round(proportion * pdf["n"].iloc[i])
                if proportion >= 0.05 and count > 0:
                    ax.text(
                        left + proportion / 2, i, str(count),
                        ha="center", va="center", fontsize=8, color="white"
                    )
            lefts = [left + val for left, val in zip(lefts, values, strict=True)]

            proportion = total_row.get(outcome, 0.0)
            count = total_counts.get(outcome, 0)
            ax.barh(total_y, proportion, left=total_left, color=colors[outcome], edgecolor="white")
            if proportion >= 0.05 and count > 0:
                ax.text(
                    total_left + proportion / 2, total_y, str(count),
                    ha="center", va="center", fontsize=8, color="white"
                )
            total_left += proportion

        ax.axhline(y=total_y - 0.5, color="grey", linewidth=0.8, linestyle="--")
        all_yticks = [*range(len(pdf)), total_y]
        all_labels = [*pdf.index.tolist(), "Total"]
        ax.set_yticks(all_yticks)
        ax.set_yticklabels(all_labels)
        ax.invert_yaxis()
        for i, rank_bin in enumerate(pdf.index):
            ax.text(1.02, i, f"n={pdf.loc[rank_bin, 'n']}", va="center", fontsize=8)
        ax.text(1.02, total_y, f"n={total_n}", va="center", fontsize=8)

        ax.set_xlabel("Proportion")
        ax.set_ylabel("Rank")
        ax.set_xlim(0, 1)
        ax.set_title(f"Rank changes: {run1} \u2192 {run2}\n" f"({benchmark_output_type.prioritisation_type_string})")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
        fig.savefig(
            self.output_dir.joinpath(
                f"{self.benchmark_name}_{run1}_vs_{run2}"
                f"_{benchmark_output_type.prioritisation_type_string}_rank_changes.svg"
            ),
            format="svg",
            bbox_inches="tight",
        )
        plt.close(fig)

    def generate_rank_change_plot(
        self,
        rank_changes_df: pl.DataFrame,
        run1: str,
        run2: str,
        benchmark_output_type: BenchmarkOutputType,
    ) -> None:
        """
        Generate a horizontal stacked bar plot showing how ranks changed between two runs.

        Args:
            rank_changes_df: DataFrame with run1/run2 rank columns and a rank_change column.
            run1: Baseline run identifier.
            run2: Comparison run identifier.
            benchmark_output_type: Type of benchmark output.
        """
        classified_rank_changes = self._classify_rank_changes(rank_changes_df, run1, run2)
        plot_data = self._compute_rank_change_plot_data(classified_rank_changes)
        if plot_data.is_empty():
            logger.warning(f"No rank data to plot for {run1} vs {run2}.")
            return
        self._plot_rank_change_bar(plot_data, run1, run2, benchmark_output_type)

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
        palette = self.get_palette(len(curves))
        for i, row in enumerate(curves.iter_rows(named=True)):
            run_identifier = row["run_identifier"]
            fpr = row["fpr"]
            tpr = row["tpr"]
            roc_auc = auc(fpr, tpr)
            plt.plot(
                fpr,
                tpr,
                label=f"{run_identifier} ROC Curve (AUC = {roc_auc:.2f})",
                color=palette[i],
            )
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(plot_customisation.roc_curve_title)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            self.output_dir.joinpath(
                f"{self.benchmark_name}_{benchmark_generator.prioritisation_type_string}_roc_curve.svg"
            ),
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
        palette = self.get_palette(len(curves))
        for i, row in enumerate(curves.iter_rows(named=True)):
            run_identifier = row["run_identifier"]
            precision = row["precision"]
            recall = row["recall"]
            pr_auc = auc(recall[::-1], precision[::-1])
            plt.plot(
                recall,
                precision,
                label=f"{run_identifier} Precision-Recall Curve (AUC = {pr_auc:.2f})",
                color=palette[i],
            )
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(plot_customisation.precision_recall_title)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15))
        plt.savefig(
            self.output_dir.joinpath(
                f"{self.benchmark_name}_{benchmark_generator.prioritisation_type_string}_pr_curve.svg"
            ),
            format="svg",
            bbox_inches="tight",
        )


def generate_plots(
    benchmark_name: str,
    benchmarking_results_df: pl.DataFrame,
    curves: pl.DataFrame,
    benchmark_output_type: BenchmarkOutputType,
    plot_customisation: PlotCustomisation,
    output_dir: Path,
    no_curves: bool,
    conn: DuckDBPyConnection,
    run_identifiers: list[str],
) -> None:
    """
    Generate all plots for a benchmarking run.

    Generates rank summary bar plots, and optionally ROC/PR curves and pairwise
    rank change plots when a database connection and run identifiers are supplied.
    """
    plot_generator = PlotGenerator(benchmark_name, output_dir)
    plot_customisation_type = getattr(plot_customisation, f"{benchmark_output_type.prioritisation_type_string}_plots")
    if not no_curves:
        logger.info("Generating ROC curve visualisations.")
        plot_generator.generate_roc_curve(curves, benchmark_output_type, plot_customisation_type)
        logger.info("Generating Precision-Recall curves visualisations.")
        plot_generator.generate_precision_recall(curves, benchmark_output_type, plot_customisation_type)
    if no_curves:
        logger.info("No ROC curve visualisations generated.")
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
    if len(run_identifiers) >= 2:
        for run1, run2 in combinations(run_identifiers, 2):
            table_name = f"{run1}_vs_{run2}_{benchmark_output_type.prioritisation_type_string}_rank_changes"
            logger.info(f"Generating rank change plot for {run1} vs {run2}.")
            rank_changes_df = load_table_lazy(table_name, conn).collect()
            plot_generator.generate_rank_change_plot(rank_changes_df, run1, run2, benchmark_output_type)


def generate_plots_from_db(db_path: Path, config: Path, output_dir: Path) -> None:
    """
    Generate plots from database file.
    Args:
        db_path (Path): Path to the database file.
        config (Path): Path to the benchmarking config file.
        output_dir (Path): Path to the output directory.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
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
            f"{benchmark_config_file.benchmark_name}_{benchmark_output_type.value.prioritisation_type_string}_summary"
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
                output_dir=output_dir,
                conn=conn,
                run_identifiers=[run.run_identifier for run in benchmark_config_file.runs],
            )
    logger.info("Finished generating plots.")
    conn.close()
