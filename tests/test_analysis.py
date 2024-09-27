import unittest
from copy import copy
from pathlib import Path
from unittest.mock import patch

import duckdb

from pheval.analyse.benchmark_db_manager import BenchmarkDBManager
from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.disease_prioritisation_analysis import AssessDiseasePrioritisation
from pheval.analyse.gene_prioritisation_analysis import AssessGenePrioritisation
from pheval.analyse.variant_prioritisation_analysis import AssessVariantPrioritisation
from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)


class TestAssessGenePrioritisation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db_connection = duckdb.connect(":memory:")
        cls.db_connection.execute(
            "CREATE TABLE test_table_gene (identifier VARCHAR(255) PRIMARY KEY, "
            "phenopacket VARCHAR, gene_symbol VARCHAR, gene_identifier VARCHAR)"
        )
        cls.db_connection.execute(
            "INSERT INTO test_table_gene (identifier, phenopacket, gene_symbol, gene_identifier) VALUES "
            "('phenopacket_1.json-PLXNA1', 'phenopacket_1.json', 'PLXNA1', 'ENSG00000114554'),"
            "('phenopacket_1.json-LARGE1', 'phenopacket_1.json', 'LARGE1', 'ENSG00000133424'),"
        )
        cls.db_connection.execute(
            "CREATE TABLE result (rank INTEGER, score DOUBLE, gene_symbol VARCHAR, gene_identifier VARCHAR)"
        )
        cls.db_connection.execute(
            "INSERT INTO result (rank, score, gene_symbol, gene_identifier) VALUES "
            "(1, 0.8764, 'PLXNA1', 'ENSG00000114554'),"
            "(2, 0.5777, 'ZNF804B', 'ENSG00000182348'),"
            "(2, 0.5777, 'SMCO2', 'ENSG00000165935'),"
            "(4, 0.3765, 'SPNS1', 'ENSG00000169682')"
        )

    @classmethod
    def tearDownClass(cls):
        cls.db_connection.close()

    def setUp(self):
        patcher = patch(
            "pheval.analyse.benchmark_db_manager.BenchmarkDBManager.get_connection",
            return_value=self.db_connection,
        )
        self.mock_get_connection = patcher.start()
        self.addCleanup(patcher.stop)
        self.db_connector = BenchmarkDBManager(None)
        self.assess_gene_prioritisation = AssessGenePrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_gene",
            column="run_1",
            threshold=0,
            score_order="descending",
        )
        self.assess_gene_prioritisation_ascending_order = AssessGenePrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_gene",
            column="run_1",
            threshold=0,
            score_order="ascending",
        )
        self.binary_classification_stats = BinaryClassificationStats()

    def test_assess_gene_with_ascending_order_threshold_fails_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.1
        self.assertEqual(
            assess_ascending_order_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_gene_with_ascending_order_threshold_meets_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.9
        self.assertEqual(
            assess_ascending_order_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_gene_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_gene_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.5
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_gene_prioritisation_no_threshold(self):
        self.db_connector.add_contains_function()
        self.assess_gene_prioritisation.assess_gene_prioritisation(
            "result",
            Path("/path/to/phenopacket_1.json"),
            self.binary_classification_stats,
        )
        self.db_connector.conn.execute("SELECT * FROM test_table_gene")
        self.assertEqual(
            self.db_connector.conn.fetchall(),
            [
                ("phenopacket_1.json-PLXNA1", "phenopacket_1.json", "PLXNA1", "ENSG00000114554", 1),
                ("phenopacket_1.json-LARGE1", "phenopacket_1.json", "LARGE1", "ENSG00000133424", 0),
            ],
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1,
                true_negatives=3,
                false_positives=0,
                false_negatives=0,
                labels=[1, 0, 0, 0],
                scores=[0.8764, 0.5777, 0.5777, 0.3765],
            ),
        )


class TestAssessVariantPrioritisation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db_connection = duckdb.connect(":memory:")
        cls.db_connection.execute(
            "CREATE TABLE test_table_variant (identifier VARCHAR(255) PRIMARY KEY,"
            "phenopacket VARCHAR, chrom VARCHAR, pos INTEGER, ref VARCHAR, alt VARCHAR)"
        )
        cls.db_connection.execute(
            "INSERT INTO test_table_variant (identifier, phenopacket, chrom, pos, ref, alt) VALUES "
            "('phenopacket_1.json-3-126741108-G-C', 'phenopacket_1.json', '3', 126741108, 'G', 'C'),"
            "('phenopacket_1.json-16-133564345-C-T', 'phenopacket_1.json', '16', 133564345, 'C', 'T'),"
        )
        cls.db_connection.execute(
            "CREATE TABLE result (rank INTEGER, score DOUBLE,"
            'chromosome VARCHAR, start INTEGER, "end" INTEGER, ref VARCHAR, alt VARCHAR)'
        )
        cls.db_connection.execute(
            'INSERT INTO result (rank, score, chromosome, start, "end", ref, alt) VALUES '
            "(1, 0.0484, '3', 126730873, 126730873 ,'G', 'A'),"
            "(1, 0.0484, '3', 126730873, 126730873 ,'G', 'T'),"
            "(3, 0.0484, '3', 126741108, 126741108 ,'G', 'C'),"
        )

    @classmethod
    def tearDownClass(cls):
        cls.db_connection.close()

    def setUp(self):
        patcher = patch(
            "pheval.analyse.benchmark_db_manager.BenchmarkDBManager.get_connection",
            return_value=self.db_connection,
        )
        self.mock_get_connection = patcher.start()
        self.addCleanup(patcher.stop)
        self.db_connector = BenchmarkDBManager("None")
        self.assess_variant_prioritisation = AssessVariantPrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_variant",
            column="run_1",
            threshold=0,
            score_order="descending",
        )
        self.assess_variant_prioritisation_ascending_order = AssessVariantPrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_variant",
            column="run_1",
            threshold=0,
            score_order="ascending",
        )
        self.binary_classification_stats = BinaryClassificationStats()

    def test_assess_variant_with_ascending_order_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_variant_with_ascending_order_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_variant_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_variant_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_variant_prioritisation(self):
        self.db_connector.add_contains_function()
        self.assess_variant_prioritisation.assess_variant_prioritisation(
            "result",
            Path("/path/to/phenopacket_1.json"),
            self.binary_classification_stats,
        )
        self.db_connector.conn.execute("SELECT * FROM test_table_variant")
        self.assertEqual(
            self.db_connector.conn.fetchall(),
            [
                (
                    "phenopacket_1.json-3-126741108-G-C",
                    "phenopacket_1.json",
                    "3",
                    126741108,
                    "G",
                    "C",
                    3,
                ),
                (
                    "phenopacket_1.json-16-133564345-C-T",
                    "phenopacket_1.json",
                    "16",
                    133564345,
                    "C",
                    "T",
                    0,
                ),
            ],
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0,
                true_negatives=0,
                false_positives=2,
                false_negatives=1,
                labels=[0, 0, 1],
                scores=[0.0484, 0.0484, 0.0484],
            ),
        )


class TestAssessDiseasePrioritisation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.db_connection = duckdb.connect(":memory:")
        cls.db_connection.execute(
            "CREATE TABLE test_table_disease (identifier VARCHAR(255) PRIMARY KEY, "
            "phenopacket VARCHAR, disease_identifier VARCHAR, disease_name VARCHAR)"
        )
        cls.db_connection.execute(
            "INSERT INTO test_table_disease (identifier, phenopacket, disease_identifier, disease_name) VALUES "
            "('phenopacket_1.json-OMIM:231670', 'phenopacket_1.json', 'OMIM:231670', 'Glutaric aciduria type 1'),"
        )
        cls.db_connection.execute(
            "CREATE TABLE result (rank INTEGER, score DOUBLE,"
            "disease_identifier VARCHAR, disease_name VARCHAR)"
        )
        cls.db_connection.execute(
            "INSERT INTO result (rank, score, disease_identifier, disease_name) VALUES "
            "(1, 1.0, 'OMIM:231670', 'Glutaric aciduria type 1'),"
            "(2, 0.5, 'OMIM:231680', 'Glutaric aciduria type 2'),"
            "(2, 0.5, 'OMIM:231690', 'Glutaric aciduria type 3'),"
            "(4, 0.3, 'OMIM:231700', 'Glutaric aciduria type 4'),"
        )

    @classmethod
    def tearDownClass(cls):
        cls.db_connection.close()

    def setUp(self):
        patcher = patch(
            "pheval.analyse.benchmark_db_manager.BenchmarkDBManager.get_connection",
            return_value=self.db_connection,
        )
        self.mock_get_connection = patcher.start()
        self.addCleanup(patcher.stop)
        self.db_connector = BenchmarkDBManager("None")
        self.assess_disease_prioritisation = AssessDiseasePrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_disease",
            column="run_1",
            threshold=0,
            score_order="descending",
        )
        self.assess_disease_prioritisation_ascending_order = AssessDiseasePrioritisation(
            db_connection=self.db_connector,
            table_name="test_table_disease",
            column="run_1",
            threshold=0,
            score_order="ascending",
        )
        self.binary_classification_stats = BinaryClassificationStats()

    def test_assess_disease_with_ascending_order_threshold_fails_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.1
        self.assertEqual(
            assess_ascending_order_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_disease_with_ascending_order_threshold_meets_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.9
        self.assertEqual(
            assess_ascending_order_threshold._assess_with_threshold_ascending_order(
                result_entry=RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_disease_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalDiseaseResult(
                    disease_identifier="OMIM:231670",
                    disease_name="Glutaric aciduria type 1",
                    score=0.8764,
                    rank=1,
                ),
            ),
            0,
        )

    def test_assess_disease_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 0.5
        self.assertEqual(
            assess_with_threshold._assess_with_threshold(
                result_entry=RankedPhEvalDiseaseResult(
                    disease_identifier="OMIM:231670",
                    disease_name="Glutaric aciduria type 1",
                    score=0.8764,
                    rank=1,
                ),
            ),
            1,
        )

    def test_assess_disease_prioritisation(self):
        self.db_connector.add_contains_function()
        self.assess_disease_prioritisation.assess_disease_prioritisation(
            "result",
            Path("/path/to/phenopacket_1.json"),
            self.binary_classification_stats,
        )
        self.db_connector.conn.execute("SELECT * FROM test_table_disease")
        self.assertEqual(
            self.db_connector.conn.fetchall(),
            [
                (
                    "phenopacket_1.json-OMIM:231670",
                    "phenopacket_1.json",
                    "OMIM:231670",
                    "Glutaric aciduria type 1",
                    1,
                )
            ],
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1,
                true_negatives=3,
                false_positives=0,
                false_negatives=0,
                labels=[1, 0, 0, 0],
                scores=[1.0, 0.5, 0.5, 0.3],
            ),
        )
