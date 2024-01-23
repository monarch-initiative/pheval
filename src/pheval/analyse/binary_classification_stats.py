from dataclasses import dataclass
from typing import List, Union

from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)


@dataclass
class BinaryClassificationStats:
    """
    A data class representing counts of different categories in binary classification.

    Attributes:
        true_positives (int): The count of true positive instances - i.e., the number of known entities
            ranked 1 in the results.
        true_negatives (int): The count of true negative instances - i.e., the number of non-relevant entities
            ranked at a position other than 1 in the results.
        false_positives (int): The count of false positive instances - i.e., the number of non-relevant entities
            ranked at position 1 in the results.
        false_negatives (int): The count of false negative instances - i.e., the number of known entities
            ranked at a position other than 1 in the results.
    """

    true_positives: int = 0
    true_negatives: int = 0
    false_positives: int = 0
    false_negatives: int = 0

    @staticmethod
    def remove_relevant_ranks(
        pheval_results: Union[
            List[RankedPhEvalGeneResult],
            List[RankedPhEvalVariantResult],
            List[RankedPhEvalDiseaseResult],
        ],
        relevant_ranks: List[int],
    ) -> List[int]:
        """
        Remove the relevant entity ranks from all result ranks
        Args:
            pheval_results:
                (Union[List[RankedPhEvalGeneResult], List[RankedPhEvalVariantResult], List[RankedPhEvalDiseaseResult]]):
                    The list of all pheval results.
            relevant_ranks (List[int]): A list of the ranks associated with the known entities.

        Returns:
            List[int]: A list of the ranks with the relevant entity ranks removed.

        """
        all_result_ranks = [pheval_result.rank for pheval_result in pheval_results]
        for rank in relevant_ranks:
            if rank in all_result_ranks:
                all_result_ranks.remove(rank)
                continue
        return all_result_ranks

    def add_classification_for_known_entities(self, relevant_ranks: List[int]) -> None:
        """
        Update binary classification metrics for known entities based on their ranking.

        Args:
            relevant_ranks (List[int]): A list of the ranks associated with the known entities.
        """
        for rank in relevant_ranks:
            if rank == 1:
                self.true_positives += 1
            elif rank != 1:
                self.false_negatives += 1

    def add_classification_for_other_entities(self, ranks: List[int]) -> None:
        """
        Update binary classification metrics for other entities based on their ranking.

        Args:
            ranks (List[int]): A list of the ranks for all other entities.
        """
        for rank in ranks:
            if rank == 1:
                self.false_positives += 1
            elif rank != 1:
                self.true_negatives += 1

    def add_classification(
        self,
        pheval_results: Union[
            List[RankedPhEvalGeneResult],
            List[RankedPhEvalVariantResult],
            List[RankedPhEvalDiseaseResult],
        ],
        relevant_ranks: List[int],
    ) -> None:
        """
        Update binary classification metrics for known and unknown entities based on their ranks.
        Args:
            pheval_results:
                (Union[List[RankedPhEvalGeneResult], List[RankedPhEvalVariantResult], List[RankedPhEvalDiseaseResult]]):
                    The list of all pheval results.
            relevant_ranks (List[int]): A list of the ranks associated with the known entities.
        """
        self.add_classification_for_known_entities(relevant_ranks)
        self.add_classification_for_other_entities(
            self.remove_relevant_ranks(pheval_results, relevant_ranks)
        )

    def sensitivity(self) -> float:
        """
        Calculate sensitivity.

        Sensitivity measures the proportion of actual positive instances correctly identified by the model.

        Returns:
            float: The sensitivity of the model, calculated as true positives divided by the sum of true positives
            and false negatives. Returns 0 if both true positives and false negatives are zero.
        """
        return (
            self.true_positives / (self.true_positives + self.false_negatives)
            if (self.true_positives + self.false_negatives) > 0
            else 0.0
        )

    def specificity(self) -> float:
        """
        Calculate specificity.

        Specificity measures the proportion of actual negative instances correctly identified by the model.

        Returns:
            float: The specificity of the model, calculated as true negatives divided by the sum of true negatives
            and false positives. Returns 0.0 if both true negatives and false positives are zero.
        """
        return self.true_negatives / (self.true_negatives + self.false_positives) if (self.true_negatives + self.false_positives) > 0 else 0.0
