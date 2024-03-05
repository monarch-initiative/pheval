from dataclasses import dataclass, field
from math import sqrt
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
    labels: List = field(default_factory=list)
    scores: List = field(default_factory=list)

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

    def add_labels_and_scores(
        self,
        pheval_results: Union[
            List[RankedPhEvalGeneResult],
            List[RankedPhEvalVariantResult],
            List[RankedPhEvalDiseaseResult],
        ],
        relevant_ranks: List[int],
    ):
        """
        Adds scores and labels from the PhEval results.

        Args:
            pheval_results (Union[List[RankedPhEvalGeneResult], List[RankedPhEvalVariantResult],
                                  List[RankedPhEvalDiseaseResult]]):
                List of all PhEval results
            relevant_ranks (List[int]): A list of the ranks associated with the known entities.
        """
        relevant_ranks_copy = relevant_ranks.copy()
        for result in pheval_results:
            self.scores.append(result.score)
            label = 1 if result.rank in relevant_ranks_copy else 0
            self.labels.append(label)
            relevant_ranks_copy.remove(result.rank) if label == 1 else None

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
        self.add_labels_and_scores(pheval_results, relevant_ranks)

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
        return (
            self.true_negatives / (self.true_negatives + self.false_positives)
            if (self.true_negatives + self.false_positives) > 0
            else 0.0
        )

    def precision(self) -> float:
        """
        Calculate precision.

        Precision measures the proportion of correctly predicted positive instances out of all instances
        predicted as positive.

        Returns:
            float: The precision of the model, calculated as true positives divided by the sum of true positives
            and false positives. Returns 0.0 if both true positives and false positives are zero.
        """
        return (
            self.true_positives / (self.true_positives + self.false_positives)
            if (self.true_positives + self.false_positives) > 0
            else 0.0
        )

    def negative_predictive_value(self) -> float:
        """
        Calculate Negative Predictive Value (NPV).

        NPV measures the proportion of correctly predicted negative instances out of all instances predicted negative.

        Returns:
            float: The Negative Predictive Value of the model, calculated as true negatives divided by the sum of
            true negatives and false negatives. Returns 0.0 if both true negatives and false negatives are zero.
        """
        return (
            self.true_negatives / (self.true_negatives + self.false_negatives)
            if (self.true_negatives + self.false_negatives) > 0
            else 0.0
        )

    def false_positive_rate(self) -> float:
        """
        Calculate False Positive Rate (FPR).

        FPR measures the proportion of instances predicted as positive that are actually negative.

        Returns:
            float: The False Positive Rate of the model, calculated as false positives divided by the sum of
            false positives and true negatives. Returns 0.0 if both false positives and true negatives are zero.
        """
        return (
            self.false_positives / (self.false_positives + self.true_negatives)
            if (self.false_positives + self.true_negatives) > 0
            else 0.0
        )

    def false_discovery_rate(self) -> float:
        """
        Calculate False Discovery Rate (FDR).

        FDR measures the proportion of instances predicted as positive that are actually negative.

        Returns:
            float: The False Discovery Rate of the model, calculated as false positives divided by the sum of
            false positives and true positives. Returns 0.0 if both false positives and true positives are zero.
        """
        return (
            self.false_positives / (self.false_positives + self.true_positives)
            if (self.false_positives + self.true_positives) > 0
            else 0.0
        )

    def false_negative_rate(self) -> float:
        """
        Calculate False Negative Rate (FNR).

        FNR measures the proportion of instances that are actually positive but predicted as negative.

        Returns:
            float: The False Negative Rate of the model, calculated as false negatives divided by the sum of
            false negatives and true positives. Returns 0.0 if both false negatives and true positives are zero.
        """
        return (
            self.false_negatives / (self.false_negatives + self.true_positives)
            if (self.false_negatives + self.true_positives) > 0
            else 0.0
        )

    def accuracy(self) -> float:
        """
        Calculate Accuracy.

        Accuracy measures the proportion of correctly predicted instances out of all instances.

        Returns:
            float: The Accuracy of the model, calculated as the sum of true positives and true negatives divided by
            the sum of true positives, false positives, true negatives, and false negatives.
            Returns 0.0 if the total sum of counts is zero.
        """
        return (
            (self.true_positives + self.true_negatives)
            / (
                self.true_positives
                + self.false_positives
                + self.true_negatives
                + self.false_negatives
            )
            if (
                self.true_positives
                + self.false_negatives
                + self.true_negatives
                + self.false_negatives
            )
            > 0
            else 0.0
        )

    def f1_score(self) -> float:
        """
        Calculate F1 Score.

        F1 Score is the harmonic mean of precision and recall, providing a balance between false positives
        and false negatives.

        Returns:
            float: The F1 Score of the model, calculated as 2 * TP / (2 * TP + FP + FN).
            Returns 0.0 if the denominator is zero.
        """
        return (
            (2 * self.true_positives)
            / ((2 * self.true_positives) + self.false_positives + self.false_negatives)
            if (self.true_positives + self.false_positives + self.false_negatives) > 0
            else 0.0
        )

    def matthews_correlation_coefficient(self) -> float:
        """
        Calculate Matthews Correlation Coefficient (MCC).

        MCC is a measure of the quality of binary classifications, accounting for imbalances in the data.

        Returns:
            float: The Matthews Correlation Coefficient of the model, calculated as
            ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)).
            Returns 0.0 if the denominator is zero.
        """
        return (
            (
                (self.true_positives * self.true_negatives)
                - (self.false_positives * self.false_negatives)
            )
            / (
                sqrt(
                    (self.true_positives + self.false_positives)
                    * (self.true_positives + self.false_negatives)
                    * (self.true_negatives + self.false_positives)
                    * (self.true_negatives + self.false_negatives)
                )
            )
            if (
                self.true_positives
                + self.false_negatives
                + self.true_negatives
                + self.false_negatives
            )
            > 0
            else 0.0
        )
