from dataclasses import dataclass
from multiprocessing.util import get_logger

import polars as pl


@dataclass(frozen=True)
class ConfusionMatrix:
    """
    Define logical conditions for computing a confusion matrix using Polars expressions.

    Attributes:
        TRUE_POSITIVES (pl.Expr): Condition identifying true positive cases,
            where `rank == 1` and `true_positive` is `True`.
        FALSE_POSITIVES (pl.Expr): Condition identifying false positive cases,
            where `rank == 1` and `true_positive` is `False`.
        TRUE_NEGATIVES (pl.Expr): Condition identifying true negative cases,
            where `rank != 1` and `true_positive` is `False`.
        FALSE_NEGATIVES (pl.Expr): Condition identifying false negative cases,
            where `rank != 1` and `true_positive` is `True`.
    """

    TRUE_POSITIVES = (pl.col("rank") == 1) & (pl.col("true_positive"))
    FALSE_POSITIVES = (pl.col("rank") == 1) & (~pl.col("true_positive"))
    TRUE_NEGATIVES = (pl.col("rank") != 1) & (~pl.col("true_positive"))
    FALSE_NEGATIVES = (pl.col("rank") != 1) & (pl.col("true_positive"))


@dataclass(frozen=True)
class BinaryClassificationStats:
    """Binary classification statistic expressions."""

    SENSITIVITY = (
        pl.when((pl.col("true_positives") + pl.col("false_negatives")) != 0)
        .then(pl.col("true_positives") / (pl.col("true_positives") + pl.col("false_negatives")))
        .otherwise(0.0)
        .alias("sensitivity")
    )

    SPECIFICITY = (
        pl.when((pl.col("true_negatives") + pl.col("false_positives")) != 0)
        .then(pl.col("true_negatives") / (pl.col("true_negatives") + pl.col("false_positives")))
        .otherwise(0.0)
        .alias("specificity")
    )

    PRECISION = (
        pl.when((pl.col("true_positives") + pl.col("false_positives")) != 0)
        .then(pl.col("true_positives") / (pl.col("true_positives") + pl.col("false_positives")))
        .otherwise(0.0)
        .alias("precision")
    )

    NEGATIVE_PREDICTIVE_VALUE = (
        pl.when((pl.col("true_negatives") + pl.col("false_negatives")) != 0)
        .then(pl.col("true_negatives") / (pl.col("true_negatives") + pl.col("false_negatives")))
        .otherwise(0.0)
        .alias("negative_predictive_value")
    )

    FALSE_POSITIVE_RATE = (
        pl.when((pl.col("false_positives") + pl.col("true_negatives")) != 0)
        .then(pl.col("false_positives") / (pl.col("false_positives") + pl.col("true_negatives")))
        .otherwise(0.0)
        .alias("false_positive_rate")
    )

    FALSE_DISCOVERY_RATE = (
        pl.when((pl.col("false_positives") + pl.col("true_positives")) != 0)
        .then(pl.col("false_positives") / (pl.col("false_positives") + pl.col("true_positives")))
        .otherwise(0.0)
        .alias("false_discovery_rate")
    )

    FALSE_NEGATIVE_RATE = (
        pl.when((pl.col("false_negatives") + pl.col("true_positives")) != 0)
        .then(pl.col("false_negatives") / (pl.col("false_negatives") + pl.col("true_positives")))
        .otherwise(0.0)
        .alias("false_negative_rate")
    )

    ACCURACY = (
        pl.when(
            (
                pl.col("true_positives")
                + pl.col("false_positives")
                + pl.col("true_negatives")
                + pl.col("false_negatives")
            )
            != 0
        )
        .then(
            (pl.col("true_positives") + pl.col("true_negatives"))
            / (
                pl.col("true_positives")
                + pl.col("false_positives")
                + pl.col("true_negatives")
                + pl.col("false_negatives")
            )
        )
        .otherwise(0.0)
        .alias("accuracy")
    )

    F1_SCORE = (
        pl.when(
            2 * (pl.col("true_positives") + pl.col("false_positives") + pl.col("false_negatives"))
            != 0
        )
        .then(
            2
            * pl.col("true_positives")
            / (2 * pl.col("true_positives") + pl.col("false_positives") + pl.col("false_negatives"))
        )
        .otherwise(0.0)
        .alias("f1_score")
    )

    MATTHEWS_CORRELATION_COEFFICIENT = (
        pl.when(
            (
                (pl.col("true_positives") + pl.col("false_positives"))
                * (pl.col("true_positives") + pl.col("false_negatives"))
                * (pl.col("true_negatives") + pl.col("false_positives"))
                * (pl.col("true_negatives") + pl.col("false_negatives"))
            )
            > 0
        )
        .then(
            (
                (pl.col("true_positives") * pl.col("true_negatives"))
                - (pl.col("false_positives") * pl.col("false_negatives"))
            )
            / (
                (pl.col("true_positives") + pl.col("false_positives"))
                * (pl.col("true_positives") + pl.col("false_negatives"))
                * (pl.col("true_negatives") + pl.col("false_positives"))
                * (pl.col("true_negatives") + pl.col("false_negatives"))
            ).sqrt()
        )
        .otherwise(0.0)
        .alias("matthews_correlation_coefficient")
    )


def compute_confusion_matrix(run_identifier: str, result_scan: pl.LazyFrame) -> pl.LazyFrame:
    """
    Computes binary classification statistics.

    Args:
        run_identifier (str): The identifier for the run.
        result_scan (pl.LazyFrame): The LazyFrame containing the results for the directory.

    Returns:
        pl.LazyFrame: The LazyFrame containing the binary classification statistics.
    """
    logger = get_logger()
    logger.info(f"Computing binary classification statistics for {run_identifier}")
    confusion_matrix = result_scan.select(
        [
            pl.lit(run_identifier).alias("run_identifier"),
            ConfusionMatrix.TRUE_POSITIVES.sum().alias("true_positives").cast(pl.Int64),
            ConfusionMatrix.FALSE_POSITIVES.sum().alias("false_positives").cast(pl.Int64),
            ConfusionMatrix.TRUE_NEGATIVES.sum().alias("true_negatives").cast(pl.Int64),
            ConfusionMatrix.FALSE_NEGATIVES.sum().alias("false_negatives").cast(pl.Int64),
        ]
    )
    return confusion_matrix.select(
        [
            pl.col("run_identifier"),
            pl.col("true_positives"),
            pl.col("false_positives"),
            pl.col("true_negatives"),
            pl.col("false_negatives"),
            BinaryClassificationStats.SENSITIVITY,
            BinaryClassificationStats.SPECIFICITY,
            BinaryClassificationStats.PRECISION,
            BinaryClassificationStats.NEGATIVE_PREDICTIVE_VALUE,
            BinaryClassificationStats.FALSE_POSITIVE_RATE,
            BinaryClassificationStats.FALSE_DISCOVERY_RATE,
            BinaryClassificationStats.FALSE_NEGATIVE_RATE,
            BinaryClassificationStats.ACCURACY,
            BinaryClassificationStats.F1_SCORE,
            BinaryClassificationStats.MATTHEWS_CORRELATION_COEFFICIENT,
        ]
    )
