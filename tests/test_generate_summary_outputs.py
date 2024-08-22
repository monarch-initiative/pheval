import unittest
from unittest.mock import patch

import duckdb

from pheval.analyse.generate_summary_outputs import create_comparison_table, get_new_table_name
from pheval.analyse.get_connection import DBConnector


class TestGetNewTableName(unittest.TestCase):
    def test_get_new_table_name(self):
        new_table_name = get_new_table_name("run_1", "run_2", "gene")
        self.assertEqual(new_table_name, "run_1_vs_run_2_gene_rank_comparison")


class TestCreateComparisonTable(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.db_connection = duckdb.connect(":memory:")
        cls.db_connection.execute(
            "CREATE TABLE test_table_gene (identifier VARCHAR(255) PRIMARY KEY, "
            "phenopacket VARCHAR, gene_symbol VARCHAR, gene_identifier VARCHAR, "
            '"/path/to/result_dir_1" INTEGER, '
            '"/path/to/result_dir_2" INTEGER, '
            '"/path/to/result_dir_3" INTEGER)'
        )
        cls.db_connection.execute(
            "INSERT INTO test_table_gene (identifier, phenopacket, gene_symbol, gene_identifier, "
            '"/path/to/result_dir_1", "/path/to/result_dir_2", "/path/to/result_dir_3") VALUES '
            "('phenopacket_1.json-PLXNA1', 'phenopacket_1.json', 'PLXNA1', 'ENSG00000114554', 1, 0, 5),"
            "('phenopacket_1.json-LARGE1', 'phenopacket_1.json', 'LARGE1', 'ENSG00000133424', 2, 9, 0),"
        )

    @classmethod
    def tearDownClass(cls):
        cls.db_connection.close()

    def setUp(self):
        patcher = patch(
            "pheval.analyse.get_connection.DBConnector.get_connection",
            return_value=self.db_connection,
        )
        self.mock_get_connection = patcher.start()
        self.addCleanup(patcher.stop)
        self.db_connector = DBConnector("None")

    def test_create_comparison_table(self):
        create_comparison_table(
            "comparison_table_1",
            self.db_connector,
            ["/path/to/result_dir_1"],
            "/path/to/result_dir_2",
            "/path/to/result_dir_3",
            "test_table_gene",
        )
        self.db_connector.conn.execute("SELECT * FROM comparison_table_1")
        self.assertEqual(
            self.db_connector.conn.fetchall(),
            [
                ("phenopacket_1.json", "PLXNA1", "ENSG00000114554", 0, 5, "GAINED"),
                ("phenopacket_1.json", "LARGE1", "ENSG00000133424", 9, 0, "LOST"),
            ],
        )
