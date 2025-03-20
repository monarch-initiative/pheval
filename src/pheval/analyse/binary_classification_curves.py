from typing import Tuple

import numpy as np
import polars as pl
from sklearn.metrics import precision_recall_curve, roc_curve

from pheval.utils.logger import get_logger


class BinaryClassificationCurves:
    """Class for computing and storing ROC & Precision-Recall curves in Polars."""

    @staticmethod
    def _compute_finite_bounds(result_scan: pl.LazyFrame) -> Tuple[float, float]:
        """
        Compute min and max finite values in the 'score' column to handle NaN and Inf values.
        Args:
            result_scan (pl.LazyFrame): The LazyFrame containing the results for the directory.

        Returns:
            Tuple[float, float]: The (max_finite, min_finite) values for normalising scores.
        """
        return (
            result_scan.select(
                [
                    pl.col("score").filter(pl.col("score").is_finite()).max().alias("max_finite"),
                    pl.col("score").filter(pl.col("score").is_finite()).min().alias("min_finite"),
                ]
            )
            .collect()
            .row(0)
        )

    @staticmethod
    def _clean_and_extract_data(
        result_scan: pl.LazyFrame, max_finite: float, min_finite: float
    ) -> pl.LazyFrame:
        """
        Normalise the 'score' column (handling NaNs and Inf values) and extract 'true_positive' labels.

        Args:
            result_scan (pl.LazyFrame): The LazyFrame containing the results for the directory.
            max_finite (float): The maximum finite score value.
            min_finite (float): The minimum finite score value.

        Returns:
            pl.LazyFrame: A LazyFrame with cleaned 'score' and binary 'true_positive' columns.
        """
        return result_scan.with_columns(
            [
                pl.when(pl.col("score").is_nan())
                .then(0.0)
                .when(pl.col("score").is_infinite() & (pl.col("score") > 0))
                .then(max_finite)
                .when(pl.col("score").is_infinite() & (pl.col("score") < 0))
                .then(min_finite)
                .otherwise(pl.col("score"))
                .alias("score"),
                pl.when(pl.col("true_positive").is_null())
                .then(0)
                .otherwise(pl.col("true_positive").cast(pl.Int8))
                .alias("true_positive"),
            ]
        )

    @staticmethod
    def _compute_roc_pr_curves(
        run_identifier: str, labels: np.ndarray, scores: np.ndarray
    ) -> pl.LazyFrame:
        """
        Compute ROC and Precision-Recall curves.

        Args:
            labels (np.ndarray): Binary ground truth labels (0 or 1).
            scores (np.ndarray): Prediction scores.

        Returns:
            pl.LazyFrame: A LazyFrame containing the computed FPR, TPR, Precision, Recall, and Thresholds.
        """
        fpr, tpr, roc_thresholds = roc_curve(labels, scores, pos_label=1)
        precision, recall, pr_thresholds = precision_recall_curve(labels, scores, pos_label=1)

        return pl.LazyFrame(
            {
                "run_identifier": [run_identifier],
                "fpr": [fpr.tolist()],
                "tpr": [tpr.tolist()],
                "threshold_roc": [roc_thresholds.tolist()],
                "precision": [precision.tolist()],
                "recall": [recall.tolist()],
                "threshold_pr": [pr_thresholds.tolist()],
            }
        )

    @classmethod
    def process(cls, result_scan: pl.LazyFrame, run_identifier: str) -> pl.LazyFrame:
        """
        Process scores, extract true labels, compute ROC and Precision-Recall curves,
        and store results in a Polars LazyFrame with NumPy arrays.

        Args:
            result_scan (pl.LazyFrame): The LazyFrame containing the results for the directory.
            run_identifier (str): Identifier for this run.

        Returns:
            pl.LazyFrame: A LazyFrame containing ROC & PR curve data with NumPy arrays.
        """
        max_finite, min_finite = cls._compute_finite_bounds(result_scan)
        cleaned_data = (
            cls._clean_and_extract_data(result_scan, max_finite, min_finite)
            .select(["true_positive", "score"])
            .collect()
        )
        return cls._compute_roc_pr_curves(
            run_identifier,
            cleaned_data["true_positive"].to_numpy().flatten(),
            cleaned_data["score"].to_numpy().flatten(),
        )


def compute_curves(run_identifier: str, result_scan: pl.LazyFrame) -> pl.LazyFrame:
    """
    Compute ROC and Precision-Recall curves.
    Args:
        result_scan (pl.LazyFrame): The LazyFrame containing the results for the directory.
        run_identifier (str): Identifier for this run.
    Returns:
        pl.LazyFrame: LazyFrame containing the ROC & Precision-Recall curve data with NumPy arrays.
    """
    logger = get_logger()
    logger.info("Calculating ROC and Precision-Recall metrics")
    return BinaryClassificationCurves.process(result_scan, run_identifier)
