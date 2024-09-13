from typing import Union

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)


class AssessPrioritisationBase:
    def __init__(
        self,
        db_connection: BenchmarkDBManager,
        table_name: str,
        column: str,
        threshold: float,
        score_order: str,
    ):
        """
        Initialise AssessPrioritisationBase class

        Args:
            db_connection (BenchmarkDBManager): DB connection.
            table_name (str): Table name.
            column (str): Column name.
            threshold (float): Threshold for scores
            score_order (str): Score order for results, either ascending or descending

        """
        self.threshold = threshold
        self.score_order = score_order
        self.db_connection = db_connection
        self.conn = db_connection.conn
        self.column = column
        self.table_name = table_name
        db_connection.add_column_integer_default(
            table_name=table_name, column=self.column, default=0
        )

    def _assess_with_threshold_ascending_order(
        self,
        result_entry: Union[
            RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult
        ],
    ) -> int:
        """
        Record the prioritisation rank if it meets the ascending order threshold.


        Args:
            result_entry (Union[RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult]):
                Ranked PhEval result entry

        Returns:
            int: Recorded prioritisation rank
        """
        if float(self.threshold) > float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _assess_with_threshold(
        self,
        result_entry: Union[
            RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult
        ],
    ) -> int:
        """
        Record the prioritisation rank if it meets the score threshold.

        Args:
            result_entry (Union[RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult]):
                Ranked PhEval result entry

        Returns:
            int: Recorded prioritisation rank
        """
        if float(self.threshold) < float(result_entry.score):
            return result_entry.rank
        else:
            return 0

    def _record_matched_entity(
        self,
        standardised_result: Union[
            RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult
        ],
    ) -> int:
        """
        Return the rank result - handling the specification of a threshold.
        Args:
            standardised_result (Union[RankedPhEvalGeneResult, RankedPhEvalDiseaseResult, RankedPhEvalVariantResult]):
                Ranked PhEval disease result entry

        Returns:
            int: Recorded entity prioritisation rank
        """
        if float(self.threshold) == 0.0:
            return standardised_result.rank
        else:
            return (
                self._assess_with_threshold(standardised_result)
                if self.score_order != "ascending"
                else self._assess_with_threshold_ascending_order(
                    standardised_result,
                )
            )
