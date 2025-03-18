import unittest

import polars as pl

from pheval.analyse.binary_classification_stats import BinaryClassificationStats


class TestBinaryClassificationStats(unittest.TestCase):
    def setUp(self):
        """Set up a test confusion matrix DataFrame."""
        self.confusion_matrix_df = pl.DataFrame(
            {
                "true_positives": [3],
                "false_positives": [2],
                "true_negatives": [4],
                "false_negatives": [1],
            }
        )

    def test_classification_metrics(self):
        """Test the calculation of classification metrics."""
        result = self.confusion_matrix_df.select(
            [
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
        self.assertTrue(
            (
                result.equals(
                    pl.DataFrame(
                        [
                            {
                                "sensitivity": 0.75,
                                "specificity": 0.6666666666666666,
                                "precision": 0.6,
                                "negative_predictive_value": 0.8,
                                "false_positive_rate": 0.3333333333333333,
                                "false_discovery_rate": 0.4,
                                "false_negative_rate": 0.25,
                                "accuracy": 0.7,
                                "f1_score": 0.6666666666666666,
                                "matthews_correlation_coefficient": 0.408248290463863,
                            }
                        ]
                    )
                )
            )
        )
