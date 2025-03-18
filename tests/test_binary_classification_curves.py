import unittest

import numpy as np
import polars as pl

from pheval.analyse.binary_classification_curves import BinaryClassificationCurves


class TestBinaryClassificationCurves(unittest.TestCase):
    def setUp(self):
        self.test_data = pl.LazyFrame(
            {
                "score": [0.5, np.nan, 0.8, np.inf, -np.inf, 0.3],
                "true_positive": [1, 0, 1, 1, 0, None],
            }
        )

    def test_compute_finite_bounds(self):
        """Test computing finite bounds (min and max finite values)."""
        max_finite, min_finite = BinaryClassificationCurves._compute_finite_bounds(self.test_data)
        self.assertEqual(max_finite, 0.8)
        self.assertEqual(min_finite, 0.3)

    def test_clean_and_extract_data(self):
        """Test that scores are properly cleaned and NaNs/Inf are replaced."""
        max_finite, min_finite = BinaryClassificationCurves._compute_finite_bounds(self.test_data)
        cleaned = BinaryClassificationCurves._clean_and_extract_data(
            self.test_data, max_finite, min_finite
        ).collect()

        expected_cleaned_scores = [0.5, 0.0, 0.8, 0.8, 0.3, 0.3]
        expected_true_positive = [1, 0, 1, 1, 0, 0]

        self.assertListEqual(cleaned["score"].to_list(), expected_cleaned_scores)
        self.assertListEqual(cleaned["true_positive"].to_list(), expected_true_positive)

    def test_compute_roc_pr_curves(self):
        """Test that ROC and PR curves are computed correctly."""
        labels = np.array([1, 0, 1, 1, 0])
        scores = np.array([0.9, 0.1, 0.8, 0.4, 0.2])

        curves = BinaryClassificationCurves._compute_roc_pr_curves(
            "test_run", labels, scores
        ).collect()
        self.assertTrue(
            curves.equals(
                pl.DataFrame(
                    [
                        {
                            "run_identifier": "test_run",
                            "fpr": [0.0, 0.0, 0.0, 1.0],
                            "tpr": [0.0, 0.3333333333333333, 1.0, 1.0],
                            "threshold_roc": [np.inf, 0.9, 0.4, 0.1],
                            "precision": [0.6, 0.75, 1.0, 1.0, 1.0, 1.0],
                            "recall": [1.0, 1.0, 1.0, 0.6666666666666666, 0.3333333333333333, 0.0],
                            "threshold_pr": [0.1, 0.2, 0.4, 0.8, 0.9],
                        }
                    ]
                )
            )
        )

    def test_process_integration(self):
        """Test full process integration from raw data to final curves."""
        result = BinaryClassificationCurves.process(self.test_data, "test_run").collect()
        self.assertTrue(
            result.equals(
                pl.DataFrame(
                    [
                        {
                            "run_identifier": "test_run",
                            "fpr": [0.0, 0.0, 0.0, 0.6666666666666666, 1.0],
                            "tpr": [0.0, 0.6666666666666666, 1.0, 1.0, 1.0],
                            "threshold_roc": [np.inf, 0.8, 0.5, 0.3, 0.0],
                            "precision": [0.5, 0.6, 1.0, 1.0, 1.0],
                            "recall": [1.0, 1.0, 1.0, 0.6666666666666666, 0.0],
                            "threshold_pr": [0.0, 0.3, 0.5, 0.8],
                        }
                    ]
                )
            )
        )
