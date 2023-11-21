import unittest

from pheval.analyse.rank_stats import RankStats


class TestRankStats(unittest.TestCase):
    def setUp(self) -> None:
        self.rank_stats = RankStats()

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
        self.assertEqual(self.rank_stats.mean_reciprocal_rank(), 0.5)
