import unittest

from pheval.analyse.rank_stats import RankStats


class TestRankStats(unittest.TestCase):
    def setUp(self) -> None:
        self.rank_stats = RankStats()
        self.complete_rank_stats = RankStats(top=0, top3=8, top5=10, top10=15, found=20, total=40)

    def test_add_rank(self):
        self.rank_stats.add_rank(1)
        self.rank_stats.add_rank(3)
        self.rank_stats.add_rank(5)
        self.rank_stats.add_rank(7)
        self.rank_stats.add_rank(10)
        self.assertEqual(
            self.rank_stats,
            RankStats(
                top=1,
                top3=2,
                top5=3,
                top10=5,
                found=5,
                total=0,
                reciprocal_ranks=[1.0, 0.3333333333333333, 0.2, 0.14285714285714285, 0.1],
            ),
        )

    def test_percentage_rank(self):
        self.rank_stats.total = 10
        self.assertTrue(self.rank_stats.percentage_rank(3) == 30)

    def test_percentage_top(self):
        self.rank_stats.top, self.rank_stats.total = 10, 20
        self.assertEqual(self.rank_stats.percentage_top(), 50)

    def test_percentage_top3(self):
        self.rank_stats.top3, self.rank_stats.total = 30, 50
        self.assertEqual(self.rank_stats.percentage_top3(), 60)

    def test_percentage_top5(self):
        self.rank_stats.top5, self.rank_stats.total = 70, 160
        self.assertEqual(self.rank_stats.percentage_top5(), 43.75)

    def test_percentage_top10(self):
        self.rank_stats.top10, self.rank_stats.total = 100, 160
        self.assertEqual(self.rank_stats.percentage_top10(), 62.5)

    def test_percentage_found(self):
        self.rank_stats.found, self.rank_stats.total = 100, 125
        self.assertEqual(self.rank_stats.percentage_found(), 80)

    def test_percentage_difference(self):
        self.assertEqual(self.rank_stats.percentage_difference(54, 23), 31)

    def test_mean_reciprocal_rank(self):
        self.rank_stats.reciprocal_ranks = [0.2, 0.4, 0.5, 0.6, 0.8]
        self.rank_stats.total = 5
        self.assertEqual(self.rank_stats.mean_reciprocal_rank(), 0.5)

    def test_mean_reciprocal_rank_missing_cases(self):
        self.rank_stats.reciprocal_ranks = [0.2, 0.4, 0.5, 0.6, 0.8]
        self.rank_stats.total = 20
        self.assertEqual(self.rank_stats.mean_reciprocal_rank(), 0.1)

    def test_return_mean_reciprocal_rank(self):
        self.rank_stats.reciprocal_ranks = [0.2, 0.4, 0.5, 0.6, 0.8]
        self.assertEqual(self.rank_stats.return_mean_reciprocal_rank(), 0.5)

    def test_return_mean_reciprocal_rank_from_mrr_variable(self):
        self.rank_stats.mrr = 0.1
        self.assertEqual(self.rank_stats.return_mean_reciprocal_rank(), 0.1)

    def test_precision_at_k_1(self):
        self.assertEqual(self.complete_rank_stats.precision_at_k(1), 0)

    def test_precision_at_k_3(self):
        self.assertEqual(self.complete_rank_stats.precision_at_k(3), 0.06666666666666667)

    def test_precision_at_k_5(self):
        self.assertEqual(self.complete_rank_stats.precision_at_k(5), 0.05)

    def test_precision_at_k_10(self):
        self.assertEqual(self.complete_rank_stats.precision_at_k(10), 0.0375)

    def test_f_beta_score_at_k_1(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(0, 1), 0)

    def test_f_beta_score_at_k_3(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(20, 3), 0.1)

    def test_f_beta_score_at_k_5(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(25, 5), 0.08333333333333334)

    def test_f_beta_score_at_k_10(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(37.5, 10), 0.06818181818181818)
