import unittest

import numpy as np
import polars as pl

from pheval.analyse.rank_stats import Ranks


class TestRanks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Set up test data."""
        cls.test_df = pl.LazyFrame(
            {
                "file_path": [
                    "file1",
                    "file1",
                    "file2",
                    "file3",
                    "file3",
                    "file4",
                    "file5",
                    "file6",
                ],
                "rank": [1, 3, 4, 1, 6, 10, 12, 0],
                "true_positive": [True, True, True, True, True, True, True, True],
            }
        )
        cls.result = pl.DataFrame(
            [
                {
                    "top1": 2,
                    "top3": 3,
                    "top5": 4,
                    "top10": 6,
                    "found": 7,
                    "total": 8,
                    "number_of_samples": 6,
                }
            ]
        )

    def test_top_k_counts(self):
        result = self.test_df.select(
            [
                Ranks.TOP_1,
                Ranks.TOP_3,
                Ranks.TOP_5,
                Ranks.TOP_10,
                Ranks.FOUND,
                Ranks.TOTAL,
                Ranks.NUMBER_OF_SAMPLES,
            ]
        ).collect()
        self.assertTrue((result.equals(self.result)))

    def test_mrr(self):
        """Test Mean Reciprocal Rank (MRR) calculation."""
        self.assertAlmostEqual(self.test_df.select(Ranks.MRR).collect().item(), 0.366667, places=6)

    def test_percentage_at_k(self):
        """Test percentage calculations at K."""
        result = self.result.select(
            [
                Ranks.percentage_at_k(1),
                Ranks.percentage_at_k(3),
                Ranks.percentage_at_k(5),
                Ranks.percentage_at_k(10),
            ]
        )

        self.assertTrue(
            result.equals(
                pl.DataFrame(
                    [
                        {
                            "percentage@1": 25.0,
                            "percentage@3": 37.5,
                            "percentage@5": 50.0,
                            "percentage@10": 75.0,
                        }
                    ]
                )
            )
        )

    def test_precision_at_k(self):
        """Test precision calculations at K."""
        result = self.result.select(
            [
                Ranks.precision_at_k(1),
                Ranks.precision_at_k(3),
                Ranks.precision_at_k(5),
                Ranks.precision_at_k(10),
            ]
        )

        self.assertTrue(
            result.equals(
                pl.DataFrame(
                    [
                        {
                            "precision@1": 0.3333333333333333,
                            "precision@3": 0.16666666666666666,
                            "precision@5": 0.13333333333333333,
                            "precision@10": 0.1,
                        }
                    ]
                )
            )
        )

    def test_f_beta_score_at_k(self):
        result = self.result.select(
            [
                Ranks.f_beta_score_at_k(1),
                Ranks.f_beta_score_at_k(3),
                Ranks.f_beta_score_at_k(5),
                Ranks.f_beta_score_at_k(10),
            ]
        )
        self.assertTrue(
            result.equals(
                pl.DataFrame(
                    [
                        {
                            "f_beta@1": 0.28571428571428575,
                            "f_beta@3": 0.23076923076923078,
                            "f_beta@5": 0.2105263157894737,
                            "f_beta@10": 0.17647058823529416,
                        }
                    ]
                )
            )
        )

    def test__average_precision_at_k(self):
        self.assertTrue(
            Ranks._average_precision_at_k(self.test_df, 1)
            .collect()
            .sort("file_path")
            .equals(
                pl.DataFrame(
                    [{"file_path": "file1", "ap@1": 1.0}, {"file_path": "file3", "ap@1": 1.0}]
                )
            )
        )
        self.assertTrue(
            Ranks._average_precision_at_k(self.test_df, 3)
            .collect()
            .sort("file_path")
            .equals(
                pl.DataFrame(
                    [
                        {"file_path": "file1", "ap@3": 0.8333333333333333},
                        {"file_path": "file3", "ap@3": 1.0},
                    ]
                )
            )
        )
        self.assertTrue(
            Ranks._average_precision_at_k(self.test_df, 5)
            .collect()
            .sort("file_path")
            .equals(
                pl.DataFrame(
                    [
                        {"file_path": "file1", "ap@5": 0.8333333333333333},
                        {"file_path": "file2", "ap@5": 0.25},
                        {"file_path": "file3", "ap@5": 1.0},
                    ]
                )
            )
        )
        self.assertTrue(
            Ranks._average_precision_at_k(self.test_df, 10)
            .collect()
            .sort("file_path")
            .equals(
                pl.DataFrame(
                    [
                        {"file_path": "file1", "ap@10": 0.8333333333333333},
                        {"file_path": "file2", "ap@10": 0.25},
                        {"file_path": "file3", "ap@10": 0.6666666666666666},
                        {"file_path": "file4", "ap@10": 0.1},
                    ]
                )
            )
        )

    def test__compute_ap_k(self):
        self.assertEqual(Ranks._compute_ap_k(np.array([1])), 1)
        self.assertEqual(Ranks._compute_ap_k(np.array([1, 5])), 0.7)

    def test_mean_average_precision_at_k(self):
        self.assertAlmostEqual(Ranks.mean_average_precision_at_k(self.test_df, 1), 0.333, places=3)
        self.assertAlmostEqual(Ranks.mean_average_precision_at_k(self.test_df, 3), 0.306, places=3)
        self.assertAlmostEqual(Ranks.mean_average_precision_at_k(self.test_df, 5), 0.347, places=3)
        self.assertAlmostEqual(Ranks.mean_average_precision_at_k(self.test_df, 10), 0.308, places=3)

    def test__calculate_ndcg_at_k(self):
        self.assertEqual(Ranks._calculate_ndcg_at_k([1], 3), 1)
        self.assertAlmostEqual(Ranks._calculate_ndcg_at_k([1, 2, 4], 5), 0.858, places=3)

    def test_mean_normalised_discounted_cumulative_gain(self):
        self.assertAlmostEqual(
            Ranks.mean_normalised_discounted_cumulative_gain(self.test_df, 3), 0.301, places=3
        )
