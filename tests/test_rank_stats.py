import unittest
from unittest.mock import patch

import duckdb

from pheval.analyse.rank_stats import RankStats


class TestRankStats(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db_connection = duckdb.connect(":memory:")
        cls.db_connection.execute(
            "CREATE TABLE test_table_gene (identifier VARCHAR(255) PRIMARY KEY, "
            "phenopacket VARCHAR, gene_symbol VARCHAR, gene_identifier VARCHAR, results_dir_1 INTEGER)"
        )
        cls.db_connection.execute(
            "INSERT INTO test_table_gene (identifier, phenopacket, gene_symbol, gene_identifier, results_dir_1) VALUES "
            "('phenopacket_1.json-GENE1', 'phenopacket_1.json', 'GENE1', 'GENEID1', 1),"
            "('phenopacket_2.json-GENE2', 'phenopacket_2.json', 'GENE2', 'GENEID2', 3),"
            "('phenopacket_3.json-GENE3', 'phenopacket_3.json', 'GENE3', 'GENEID3', 5),"
            "('phenopacket_4.json-GENE4', 'phenopacket_4.json', 'GENE4', 'GENEID4', 7),"
            "('phenopacket_5.json-GENE5', 'phenopacket_5.json', 'GENE5', 'GENEID5', 10),"
            "('phenopacket_2.json-GENE6', 'phenopacket_6.json', 'GENE6', 'GENEID6', 20),"
        )

    @classmethod
    def tearDownClass(cls):
        cls.db_connection.close()

    def setUp(self) -> None:
        self.rank_stats = RankStats()
        self.complete_rank_stats = RankStats(
            top=0,
            top3=8,
            top5=10,
            top10=15,
            found=20,
            total=40,
            relevant_result_ranks=[[4], [3], [6, 7], [2], [9], [20]],
        )

    @patch(
        "pheval.analyse.benchmark_db_manager.BenchmarkDBManager.get_connection",
        return_value=duckdb.connect(":memory:"),
    )
    def test_add_ranks(self, mock_get_connection):
        mock_get_connection.return_value = self.db_connection
        self.rank_stats.add_ranks("None", "test_table_gene", "results_dir_1")
        self.assertEqual(self.rank_stats.top, 1)
        self.assertEqual(self.rank_stats.top3, 2)
        self.assertEqual(self.rank_stats.top5, 3)
        self.assertEqual(self.rank_stats.top10, 5)
        self.assertEqual(self.rank_stats.found, 6)
        self.assertEqual(self.rank_stats.total, 6)
        self.assertEqual(
            self.rank_stats.reciprocal_ranks,
            [1.0, 0.3333333333333333, 0.2, 0.14285714285714285, 0.1, 0.05],
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

    def test__calculate_average_precision(self):
        self.assertEqual(self.rank_stats._average_precision_at_k(3, 0.5), 0.16666666666666666)

    def test__calculate_average_precision_0(self):
        self.assertEqual(self.rank_stats._average_precision_at_k(0, 0), 0)

    def test_mean_average_precision_at_k_1(self):
        self.assertEqual(self.complete_rank_stats.mean_average_precision_at_k(1), 0.0)

    def test_mean_average_precision_at_k_3(self):
        self.assertEqual(
            self.complete_rank_stats.mean_average_precision_at_k(3), 0.020833333333333332
        )

    def test_mean_average_precision_at_k_5(self):
        self.assertEqual(
            self.complete_rank_stats.mean_average_precision_at_k(5), 0.027083333333333334
        )

    def test_mean_average_precision_at_k_10(self):
        self.assertEqual(
            self.complete_rank_stats.mean_average_precision_at_k(10), 0.03968253968253968
        )

    def test_f_beta_score_at_k_1(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(0, 1), 0)

    def test_f_beta_score_at_k_3(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(20, 3), 0.1)

    def test_f_beta_score_at_k_5(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(25, 5), 0.08333333333333334)

    def test_f_beta_score_at_k_10(self):
        self.assertEqual(self.complete_rank_stats.f_beta_score_at_k(37.5, 10), 0.06818181818181818)

    def test_mean_normalised_discounted_cumulative_gain_3(self):
        self.assertEqual(
            self.complete_rank_stats.mean_normalised_discounted_cumulative_gain(3),
            0.09424414613095478,
        )

    def test_mean_normalised_discounted_cumulative_gain_5(self):
        self.assertEqual(
            self.complete_rank_stats.mean_normalised_discounted_cumulative_gain(5),
            0.243557389859924,
        )

    def test_mean_normalised_discounted_cumulative_gain_10(self):
        self.assertEqual(
            self.complete_rank_stats.mean_normalised_discounted_cumulative_gain(10),
            0.3368971541167727,
        )
