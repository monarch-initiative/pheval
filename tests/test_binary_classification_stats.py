import unittest

from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.post_processing.post_processing import RankedPhEvalDiseaseResult

disease_results = [
    RankedPhEvalDiseaseResult(
        disease_name="Glutaric aciduria type 1",
        disease_identifier="OMIM:231670",
        score=1.0,
        rank=1,
    ),
    RankedPhEvalDiseaseResult(
        disease_name="Glutaric aciduria type 2",
        disease_identifier="OMIM:231680",
        score=0.5,
        rank=2,
    ),
    RankedPhEvalDiseaseResult(
        disease_name="Glutaric aciduria type 3",
        disease_identifier="OMIM:231690",
        score=0.5,
        rank=2,
    ),
    RankedPhEvalDiseaseResult(
        disease_name="Glutaric aciduria type 4",
        disease_identifier="OMIM:231700",
        score=0.3,
        rank=4,
    ),
]


class TestBinaryClassificationStats(unittest.TestCase):
    def setUp(self):
        self.binary_classification_stats = BinaryClassificationStats()
        self.complete_binary_classification_stats = BinaryClassificationStats(10, 15, 40, 5)

    def test_remove_relevant_ranks(self):
        self.assertEqual(
            self.binary_classification_stats.remove_relevant_ranks(disease_results, [2, 4]), [1, 2]
        )

    def test_add_classification_for_known_entities_true_pos(self):
        self.binary_classification_stats.add_classification_for_known_entities([1])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=0, false_positives=0, false_negatives=0
            ),
        )

    def test_add_classification_for_known_entities_false_neg(self):
        self.binary_classification_stats.add_classification_for_known_entities([3])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=0, false_positives=0, false_negatives=1
            ),
        )

    def test_add_classification_for_known_entities_multiple_results(self):
        self.binary_classification_stats.add_classification_for_known_entities([1, 3, 100])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=0, false_positives=0, false_negatives=2
            ),
        )

    def test_add_classification_for_other_entities_false_pos(self):
        self.binary_classification_stats.add_classification_for_other_entities([1])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=0, false_positives=1, false_negatives=0
            ),
        )

    def test_add_classification_for_other_entities_true_neg(self):
        self.binary_classification_stats.add_classification_for_other_entities([10])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=1, false_positives=0, false_negatives=0
            ),
        )

    def test_add_classification_for_other_entities_multiple_results(self):
        self.binary_classification_stats.add_classification_for_other_entities([1, 2, 3, 4, 5])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=4, false_positives=1, false_negatives=0
            ),
        )

    def test_add_classification(self):
        self.binary_classification_stats.add_classification(disease_results, [1, 2])
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=2, false_positives=0, false_negatives=1
            ),
        )

    def test_sensitivity_0(self):
        self.assertEqual(self.binary_classification_stats.sensitivity(), 0)

    def test_sensitivity(self):
        self.assertEqual(
            self.complete_binary_classification_stats.sensitivity(), 0.6666666666666666
        )
