from dataclasses import dataclass
from typing import List

import numpy as np
import polars as pl
from sklearn.metrics import ndcg_score

from pheval.utils.logger import get_logger


@dataclass(frozen=True)
class Ranks:
    """
    Class for calculating ranking statistics.
    """

    TOP_1 = pl.col("rank").eq(1).sum().alias("top1")
    TOP_3 = pl.col("rank").is_between(1, 3, closed="both").sum().alias("top3")
    TOP_5 = pl.col("rank").is_between(1, 5, closed="both").sum().alias("top5")
    TOP_10 = pl.col("rank").is_between(1, 10, closed="both").sum().alias("top10")
    FOUND = pl.col("rank").gt(0).sum().alias("found")
    TOTAL = pl.len().alias("total")
    NUMBER_OF_SAMPLES = pl.col("file_path").n_unique().alias("number_of_samples")
    MRR = ((1 / pl.col("rank").filter(pl.col("rank") > 0)).sum() / pl.len()).alias("mrr")

    @classmethod
    def _filter_results(cls, df: pl.LazyFrame, k: int) -> pl.LazyFrame:
        """
        Filter for ranks within k.
        Args:
            df (pl.LazyFrame): The dataframe to filter.
            k (int): The number upper rank limit.

        Returns:
            pl.LazyFrame: The filtered dataframe.
        """
        df = df.filter(pl.col("rank").is_between(1, k, closed="both"))
        return df.group_by("file_path").agg(
            pl.col("rank").sort().alias("ranks"),
        )

    @classmethod
    def percentage_at_k(cls, k: int) -> pl.Expr:
        """
        Compute percentage at k dynamically.
        Args:
            k (int): The upper rank limit.
        Returns:
            pl.Expr: The expression for calculating percentage at k.
        """
        return (100 * pl.col(f"top{k}") / pl.col("total")).alias(f"percentage@{k}")

    @classmethod
    def percentage_found(cls) -> pl.Expr:
        """
        Compute the percentage of found items.
        Returns:
            pl.Expr: The expression for calculating percentage of found items.
        """
        return (100 * pl.col("found") / pl.col("total")).alias("percentage_found")

    @classmethod
    def precision_at_k(cls, k: int) -> pl.Expr:
        """
        Compute precision at k dynamically.
        Args:
            k (int): The upper rank limit.
        Returns:
            pl.Expr: The expression for calculating precision at k.
        """
        return (pl.col(f"top{k}") / (pl.col("number_of_samples") * k)).alias(f"precision@{k}")

    @classmethod
    def f_beta_score_at_k(cls, k: int) -> pl.Expr:
        """
        Compute f_beta_score at k.
        Args:
            k (int): The upper rank limit.
        Returns:
            pl.Expr: The expression for calculating f_beta_score at k.
        """
        precision_expr = pl.col(f"top{k}") / (pl.col("number_of_samples") * k)
        recall_expr = pl.col(f"top{k}") / pl.col("total")
        return (
            ((2 * precision_expr * recall_expr) / (precision_expr + recall_expr))
            .fill_nan(0)
            .alias(f"f_beta@{k}")
        )

    @classmethod
    def _average_precision_at_k(cls, df: pl.LazyFrame, k: int) -> pl.LazyFrame:
        """
        Compute Average Precision at K (AP@K) for each query.

        AP@K = (1 / min(k, R)) * sum(P(i) * rel(i)) for i ≤ k

        Args:
            df (pl.LazyFrame): The dataframe calculate AP@K for each query.
            k (int): The upper rank limit.
        Returns:
            pl.LazyFrame: The dataframe with AP@K for each query.
        """
        filtered_df = cls._filter_results(df, k)
        df_grouped = filtered_df.with_columns(
            pl.struct("ranks")
            .map_elements(
                lambda row: cls._compute_ap_k(np.array(row["ranks"])), return_dtype=pl.Float64
            )
            .alias(f"ap@{k}")
        )
        return df_grouped.select(["file_path", f"ap@{k}"])

    @staticmethod
    def _compute_ap_k(ranks: np.array) -> np.floating:
        """
        Helper function to compute AP@K for a single query.
        Args:
            ranks (np.array): The ranks to compute AP@K.
        Returns:
            float: The AP@K.
        """
        num_relevant = np.arange(1, len(ranks) + 1)
        precision_at_k = num_relevant / ranks
        return np.mean(precision_at_k)

    @classmethod
    def mean_average_precision_at_k(cls, df: pl.LazyFrame, k: int) -> pl.LazyFrame:
        """
        Compute Mean Average Precision at K (MAP@K) by averaging AP@K scores.
        Args:
            df (pl.LazyFrame): The dataframe calculate MAP@K for each query.
            k (int): The upper rank limit.
        Returns:
            pl.LazyFrame: The dataframe with MAP@K for each query.
        """
        ap_at_k_df = cls._average_precision_at_k(df, k)
        return (
            ap_at_k_df.select(
                pl.col(f"ap@{k}").sum() / df.select(Ranks.NUMBER_OF_SAMPLES).collect()
            )
            .fill_null(0.0)
            .collect()
            .item()
        )

    @classmethod
    def _calculate_ndcg_at_k(cls, ranks: List[int], k: int) -> float:
        """
        Compute NDCG@K for a single query.
        Args:
            ranks (List[int]): The ranks to compute NDCG@K.
            k (int): The upper rank limit.
        Returns:
            float: The NDCG@K.
        """
        result_ranks = np.zeros(k, dtype=int)
        indices = np.array(ranks) - 1
        valid_indices = indices[(indices >= 0) & (indices < k)]
        result_ranks[valid_indices] = 3
        ideal_ranking = np.sort(result_ranks)[::-1]
        return (
            ndcg_score(result_ranks.reshape(1, -1), ideal_ranking.reshape(1, -1))
            if np.sum(result_ranks) > 0
            else 0.0
        )

    @classmethod
    def mean_normalised_discounted_cumulative_gain(cls, df: pl.LazyFrame, k: int) -> pl.Float64:
        """
        Compute mean normalised discounted cumulative gain.
        Args:
            df (pl.LazyFrame): The dataframe to calculate mean normalised cumulative gain.
            k (int): The upper rank limit.
        Returns:
            pl.LazyFrame: The dataframe with mean normalised cumulative gain.
        """
        filtered_df = cls._filter_results(df, k)
        return (
            filtered_df.with_columns(
                pl.struct("ranks")
                .map_elements(
                    lambda row: cls._calculate_ndcg_at_k(row["ranks"], k), return_dtype=pl.Float64
                )
                .alias(f"NDCG@{k}")
            )
            .select(pl.col(f"NDCG@{k}").sum() / df.select(Ranks.NUMBER_OF_SAMPLES).collect())
            .fill_null(0.0)
            .collect()
            .item()
        )


def compute_rank_stats(run_identifier: str, result_scan: pl.LazyFrame) -> pl.LazyFrame:
    """
    Computes ranking statistics for a given benchmarking run.
    Args:
        run_identifier (str): The identifier of the benchmarking run.
        result_scan (pl.LazyFrame): The scan of the directory to compute ranking statistics for.
    """
    logger = get_logger()
    logger.info(f"Generating ranking statistics for {run_identifier}...")
    true_positive_scan = result_scan.filter(pl.col("true_positive"))
    rankings = true_positive_scan.select(
        [
            pl.lit(run_identifier).alias("run_identifier"),
            Ranks.TOP_1.alias("top1"),
            Ranks.TOP_3.alias("top3"),
            Ranks.TOP_5.alias("top5"),
            Ranks.TOP_10.alias("top10"),
            Ranks.FOUND.alias("found"),
            Ranks.TOTAL.alias("total"),
            Ranks.NUMBER_OF_SAMPLES.alias("number_of_samples"),
            Ranks.MRR.alias("mrr"),
        ]
    )

    return rankings.select(
        [
            pl.col("run_identifier"),
            pl.col("top1"),
            pl.col("top3"),
            pl.col("top5"),
            pl.col("top10"),
            pl.col("found"),
            pl.col("total"),
            pl.col("number_of_samples"),
            pl.col("mrr"),
            Ranks.percentage_at_k(1),
            Ranks.percentage_at_k(3),
            Ranks.percentage_at_k(5),
            Ranks.percentage_at_k(10),
            Ranks.percentage_found(),
            Ranks.precision_at_k(1),
            Ranks.precision_at_k(3),
            Ranks.precision_at_k(5),
            Ranks.precision_at_k(10),
            Ranks.f_beta_score_at_k(1),
            Ranks.f_beta_score_at_k(3),
            Ranks.f_beta_score_at_k(5),
            Ranks.f_beta_score_at_k(10),
            pl.lit(Ranks.mean_average_precision_at_k(true_positive_scan, 1)).alias("MAP@1"),
            pl.lit(Ranks.mean_average_precision_at_k(true_positive_scan, 3)).alias("MAP@3"),
            pl.lit(Ranks.mean_average_precision_at_k(true_positive_scan, 5)).alias("MAP@5"),
            pl.lit(Ranks.mean_average_precision_at_k(true_positive_scan, 10)).alias("MAP@10"),
            pl.lit(Ranks.mean_normalised_discounted_cumulative_gain(true_positive_scan, 3)).alias(
                "NDCG@3"
            ),
            pl.lit(Ranks.mean_normalised_discounted_cumulative_gain(true_positive_scan, 5)).alias(
                "NDCG@5"
            ),
            pl.lit(Ranks.mean_normalised_discounted_cumulative_gain(true_positive_scan, 10)).alias(
                "NDCG@10"
            ),
        ]
    )
