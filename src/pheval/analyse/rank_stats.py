from dataclasses import dataclass, field
from statistics import mean
from typing import List

import numpy as np
from duckdb import DuckDBPyConnection
from sklearn.metrics import ndcg_score

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.binary_classification_stats import BinaryClassificationStats


@dataclass
class RankStats:
    """Store statistics related to ranking.

    Attributes:
        top (int): Count of top-ranked matches.
        top3 (int): Count of matches within the top 3 ranks.
        top5 (int): Count of matches within the top 5 ranks.
        top10 (int): Count of matches within the top 10 ranks.
        found (int): Count of found matches.
        total (int): Total count of matches.
        reciprocal_ranks (List[float]): List of reciprocal ranks.
        relevant_ranks List[List[int]]: Nested list of ranks for the known entities for all cases in a run.
        mrr (float): Mean Reciprocal Rank (MRR). Defaults to None.
    """

    top: int = 0
    top3: int = 0
    top5: int = 0
    top10: int = 0
    found: int = 0
    total: int = 0
    reciprocal_ranks: List = field(default_factory=list)
    relevant_result_ranks: List[List[int]] = field(default_factory=list)
    mrr: float = None

    def add_ranks(self, benchmark_name: str, table_name: str, column_name: str) -> None:
        """
        Add ranks to RankStats instance from table.
        Args:
            table_name (str): Name of the table to add ranks from.
            column_name (str): Name of the column to add ranks from.:
        """
        conn = BenchmarkDBManager(benchmark_name).conn
        self.top = self._execute_count_query(conn, table_name, column_name, " = 1")
        self.top3 = self._execute_count_query(conn, table_name, column_name, " BETWEEN 1 AND 3")
        self.top5 = self._execute_count_query(conn, table_name, column_name, " BETWEEN 1 AND 5")
        self.top10 = self._execute_count_query(conn, table_name, column_name, " BETWEEN 1 AND 10")
        self.found = self._execute_count_query(conn, table_name, column_name, " > 0")
        self.total = self._execute_count_query(conn, table_name, column_name, " >= 0")
        self.reciprocal_ranks = self._fetch_reciprocal_ranks(conn, table_name, column_name)
        self.relevant_result_ranks = self._fetch_relevant_ranks(conn, table_name, column_name)
        conn.close()

    @staticmethod
    def _execute_count_query(
        conn: DuckDBPyConnection, table_name: str, column_name: str, condition: str
    ) -> int:
        """
        Execute count query on table.
        Args:
            conn (DuckDBPyConnection): Connection to the database.
            table_name (str): Name of the table to execute count query on.
            column_name (str): Name of the column to execute count query on.
            condition (str): Condition to execute count query.
        Returns:
            int: Count query result.
        """
        query = f'SELECT COUNT(*) FROM {table_name} WHERE "{column_name}" {condition}'
        return conn.execute(query).fetchone()[0]

    @staticmethod
    def _fetch_reciprocal_ranks(
        conn: DuckDBPyConnection, table_name: str, column_name: str
    ) -> List[float]:
        """
        Fetch reciprocal ranks from table.
        Args:
            conn (DuckDBPyConnection): Connection to the database.
            table_name (str): Name of the table to fetch reciprocal ranks from.
            column_name (str): Name of the column to fetch reciprocal ranks from.

        Returns:
            List[float]: List of reciprocal ranks.
        """
        query = f'SELECT "{column_name}" FROM {table_name}'
        return [1 / rank[0] if rank[0] > 0 else 0 for rank in conn.execute(query).fetchall()]

    @staticmethod
    def _fetch_relevant_ranks(
        conn: DuckDBPyConnection, table_name: str, column_name: str
    ) -> List[List[int]]:
        """
        Fetch relevant ranks from table.
        Args:
            conn (DuckDBPyConnection): Connection to the database.
            table_name (str): Name of the table to fetch relevant ranks from.
            column_name (str): Name of the column to fetch relevant ranks from.

        Returns:
            List[List[int]]: List of relevant ranks.
        """
        query = (
            f'SELECT LIST("{column_name}") as values_list FROM {table_name} GROUP BY phenopacket'
        )
        return [rank[0] for rank in conn.execute(query).fetchall()]

    def percentage_rank(self, value: int) -> float:
        """
        Calculate the percentage rank.

        Args:
            value (int): The value for which the percentage rank needs to be calculated.

        Returns:
            float: The calculated percentage rank based on the provided value and the total count.
        """
        return 100 * value / self.total

    def percentage_top(self) -> float:
        """
        Calculate the percentage of top matches.

        Returns:
            float: The percentage of top matches compared to the total count.
        """
        return self.percentage_rank(self.top)

    def percentage_top3(self) -> float:
        """
        Calculate the percentage of matches within the top 3.

        Returns:
            float: The percentage of matches within the top 3 compared to the total count.
        """
        return self.percentage_rank(self.top3)

    def percentage_top5(self) -> float:
        """
        Calculate the percentage of matches within the top 5.

        Returns:
            float: The percentage of matches within the top 5 compared to the total count.
        """
        return self.percentage_rank(self.top5)

    def percentage_top10(self) -> float:
        """
        Calculate the percentage of matches within the top 10.

        Returns:
            float: The percentage of matches within the top 10 compared to the total count.
        """
        return self.percentage_rank(self.top10)

    def percentage_found(self) -> float:
        """
        Calculate the percentage of matches found.

        Returns:
            float: The percentage of matches found compared to the total count.
        """
        return self.percentage_rank(self.found)

    @staticmethod
    def percentage_difference(percentage_value_1: float, percentage_value_2: float) -> float:
        """
        Calculate the percentage difference between two percentage values.

        Args:
            percentage_value_1 (float): The first percentage value.
            percentage_value_2 (float): The second percentage value.

        Returns:
            float: The difference between the two percentage values.
        """
        return percentage_value_1 - percentage_value_2

    def mean_reciprocal_rank(self) -> float:
        """
        Calculate the Mean Reciprocal Rank (MRR) for the stored ranks.

        The Mean Reciprocal Rank is computed as the mean of the reciprocal ranks
        for the found cases.

        If the total number of cases differs from the number of found cases,
        this method extends the reciprocal ranks list with zeroes for missing cases.

        Returns:
            float: The calculated Mean Reciprocal Rank.
        """
        if len(self.reciprocal_ranks) != self.total:
            missing_cases = self.total - self.found
            self.reciprocal_ranks.extend([0] * missing_cases)
            return mean(self.reciprocal_ranks)
        return mean(self.reciprocal_ranks)

    def return_mean_reciprocal_rank(self) -> float:
        """
        Retrieve or calculate the Mean Reciprocal Rank (MRR).

        If a pre-calculated MRR value exists (stored in the 'mrr' attribute), this method returns that value.
        Otherwise, it computes the Mean Reciprocal Rank using the 'mean_reciprocal_rank' method.

        Returns:
            float: The Mean Reciprocal Rank value.
        """
        if self.mrr is not None:
            return self.mrr
        else:
            return self.mean_reciprocal_rank()

    def precision_at_k(self, k: int) -> float:
        """
        Calculate the precision at k.
        Precision at k is the ratio of relevant items in the top-k predictions to the total number of predictions.
        It measures the accuracy of the top-k predictions made by a model.

        Args:
            k (int): The number of top predictions to consider.

        Returns:
            float: The precision at k, ranging from 0.0 to 1.0.
            A higher precision indicates a better performance in identifying relevant items in the top-k predictions.
        """
        k_attr = getattr(self, f"top{k}") if k > 1 else self.top
        return k_attr / (self.total * k)

    @staticmethod
    def _average_precision_at_k(
        number_of_relevant_entities_at_k: int, precision_at_k: float
    ) -> float:
        """
        Calculate the Average Precision at k.

        Average Precision at k (AP@k) is a metric used to evaluate the precision of a ranked retrieval system.
        It measures the precision at each relevant position up to k and takes the average.

        Args:
            number_of_relevant_entities_at_k (int): The count of relevant entities in the top-k predictions.
            precision_at_k (float): The precision at k - the sum of the precision values at each relevant position.

        Returns:
            float: The Average Precision at k, ranging from 0.0 to 1.0.
                   A higher value indicates better precision in the top-k predictions.
        """
        return (
            (1 / number_of_relevant_entities_at_k) * precision_at_k
            if number_of_relevant_entities_at_k > 0
            else 0.0
        )

    def mean_average_precision_at_k(self, k: int) -> float:
        """
        Calculate the Mean Average Precision at k.

        Mean Average Precision at k (MAP@k) is a performance metric for ranked data.
        It calculates the average precision at k for each result rank and then takes the mean across all queries.

        Args:
            k (int): The number of top predictions to consider for precision calculation.

        Returns:
            float: The Mean Average Precision at k, ranging from 0.0 to 1.0.
                   A higher value indicates better performance in ranking relevant entities higher in the predictions.
        """
        cumulative_average_precision_scores = 0
        for result_ranks in self.relevant_result_ranks:
            precision_at_k, number_of_relevant_entities_at_k = 0, 0
            for rank in result_ranks:
                if 0 < rank <= k:
                    number_of_relevant_entities_at_k += 1
                    precision_at_k += number_of_relevant_entities_at_k / rank
                cumulative_average_precision_scores += self._average_precision_at_k(
                    number_of_relevant_entities_at_k, precision_at_k
                )
        return (1 / self.total) * cumulative_average_precision_scores

    def f_beta_score_at_k(self, percentage_at_k: float, k: int) -> float:
        """
        Calculate the F-beta score at k.

        The F-beta score is a metric that combines precision and recall,
        with beta controlling the emphasis on precision.
        The Beta value is set to the value of 1 to allow for equal weighting for both precision and recall.
        This method computes the F-beta score at a specific percentage threshold within the top-k predictions.

        Args:
            percentage_at_k (float): The percentage of true positive predictions within the top-k.
            k (int): The number of top predictions to consider.

        Returns:
            float: The F-beta score at k, ranging from 0.0 to 1.0.
                   A higher score indicates better trade-off between precision and recall.
        """
        precision = self.precision_at_k(k)
        recall_at_k = percentage_at_k / 100
        return (
            (2 * precision * recall_at_k) / (precision + recall_at_k)
            if (precision + recall_at_k) > 0
            else 0
        )

    def mean_normalised_discounted_cumulative_gain(self, k: int) -> float:
        """
        Calculate the mean Normalised Discounted Cumulative Gain (NDCG) for a given rank cutoff.

        NDCG measures the effectiveness of a ranking by considering both the relevance and the order of items.

        Args:
            k (int): The rank cutoff for calculating NDCG.

        Returns:
            float: The mean NDCG score across all query results.
        """
        ndcg_scores = []
        for result_ranks in self.relevant_result_ranks:
            result_ranks = [rank for rank in result_ranks if rank <= k]
            result_ranks = [3 if i in result_ranks else 0 for i in range(k)]
            ideal_ranking = sorted(result_ranks, reverse=True)
            ndcg_scores.append(ndcg_score(np.asarray([ideal_ranking]), np.asarray([result_ranks])))
        return np.mean(ndcg_scores)


class RankStatsWriter:
    """Class for writing the rank stats to a file."""

    def __init__(self, benchmark_name: str, table_name: str):
        """
        Initialise the RankStatsWriter class
        Args:
            table_name (str): Name of table to add statistics.
        """

        self.table_name = table_name
        self.benchmark_name = benchmark_name
        conn = BenchmarkDBManager(benchmark_name).conn
        conn.execute(
            f'CREATE TABLE IF NOT EXISTS "{self.table_name}" ('
            f"results_directory_path VARCHAR,"
            f"top INT,"
            f"top3 INT,"
            f"top5 INT,"
            f"top10 INT,"
            f'"found" INT,'
            f"total INT,"
            f"mean_reciprocal_rank FLOAT,"
            f"percentage_top FLOAT,"
            f"percentage_top3 FLOAT,"
            f"percentage_top5 FLOAT,"
            f"percentage_top10 FLOAT,"
            f"percentage_found FLOAT,"
            f'"precision@1" FLOAT,'
            f'"precision@3" FLOAT,'
            f'"precision@5" FLOAT,'
            f'"precision@10" FLOAT,'
            f'"MAP@1" FLOAT,'
            f'"MAP@3" FLOAT,'
            f'"MAP@5" FLOAT,'
            f'"MAP@10" FLOAT,'
            f'"f_beta_score@1" FLOAT,'
            f'"f_beta_score@3"FLOAT,'
            f'"f_beta_score@5" FLOAT,'
            f'"f_beta_score@10" FLOAT,'
            f'"NDCG@3" FLOAT,'
            f'"NDCG@5" FLOAT,'
            f'"NDCG@10" FLOAT,'
            f"true_positives INT,"
            f"false_positives INT,"
            f"true_negatives INT,"
            f"false_negatives INT,"
            f"sensitivity FLOAT,"
            f"specificity FLOAT,"
            f'"precision" FLOAT,'
            f"negative_predictive_value FLOAT,"
            f"false_positive_rate FLOAT,"
            f"false_discovery_rate FLOAT,"
            f"false_negative_rate FLOAT,"
            f"accuracy FLOAT,"
            f"f1_score FLOAT,"
            f"matthews_correlation_coefficient FLOAT,                        )"
        )
        conn.close()

    def add_statistics_entry(
        self,
        run_identifier: str,
        rank_stats: RankStats,
        binary_classification: BinaryClassificationStats,
    ):
        """
        Add statistics row to table for a run.
        Args:
            run_identifier (str): The run identifier.
            rank_stats (RankStats): RankStats object for the run.
            binary_classification (BinaryClassificationStats): BinaryClassificationStats object for the run.
        """
        conn = BenchmarkDBManager(self.benchmark_name).conn
        conn.execute(
            f' INSERT INTO "{self.table_name}" VALUES ( '
            f"'{run_identifier}',"
            f"{rank_stats.top},"
            f"{rank_stats.top3},"
            f"{rank_stats.top5},"
            f"{rank_stats.top10},"
            f"{rank_stats.found},"
            f"{rank_stats.total},"
            f"{rank_stats.mean_reciprocal_rank()},"
            f"{rank_stats.percentage_top()},"
            f"{rank_stats.percentage_top3()},"
            f"{rank_stats.percentage_top5()},"
            f"{rank_stats.percentage_top10()},"
            f"{rank_stats.percentage_found()},"
            f"{rank_stats.precision_at_k(1)},"
            f"{rank_stats.precision_at_k(3)},"
            f"{rank_stats.precision_at_k(5)},"
            f"{rank_stats.precision_at_k(10)},"
            f"{rank_stats.mean_average_precision_at_k(1)},"
            f"{rank_stats.mean_average_precision_at_k(3)},"
            f"{rank_stats.mean_average_precision_at_k(5)},"
            f"{rank_stats.mean_average_precision_at_k(10)},"
            f"{rank_stats.f_beta_score_at_k(rank_stats.percentage_top(), 1)},"
            f"{rank_stats.f_beta_score_at_k(rank_stats.percentage_top(), 3)},"
            f"{rank_stats.f_beta_score_at_k(rank_stats.percentage_top(), 5)},"
            f"{rank_stats.f_beta_score_at_k(rank_stats.percentage_top(), 10)},"
            f"{rank_stats.mean_normalised_discounted_cumulative_gain(3)},"
            f"{rank_stats.mean_normalised_discounted_cumulative_gain(5)},"
            f"{rank_stats.mean_normalised_discounted_cumulative_gain(10)},"
            f"{binary_classification.true_positives},"
            f"{binary_classification.false_positives},"
            f"{binary_classification.true_negatives},"
            f"{binary_classification.false_negatives},"
            f"{binary_classification.sensitivity()},"
            f"{binary_classification.specificity()},"
            f"{binary_classification.precision()},"
            f"{binary_classification.negative_predictive_value()},"
            f"{binary_classification.false_positive_rate()},"
            f"{binary_classification.false_discovery_rate()},"
            f"{binary_classification.false_negative_rate()},"
            f"{binary_classification.accuracy()},"
            f"{binary_classification.f1_score()},"
            f"{binary_classification.matthews_correlation_coefficient()})"
        )

        conn.close()
