import unittest
from collections import defaultdict
from pathlib import Path, PosixPath

from pheval.post_process.post_processing_analysis import (
    PrioritisationRankRecorder,
    GenePrioritisationResultData,
    RankStats,
    VariantPrioritisationResultData,
)
from pheval.utils.phenopacket_utils import VariantData


class TestComparePrioritisationForRuns(unittest.TestCase):
    def setUp(self) -> None:
        self.add_new_phenopacket_variant_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            VariantPrioritisationResultData(
                Path("/path/to/phenopacket-2.json"), VariantData("1", 4896347, "C", "T"), 9
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12_120434_A_G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_directory_variant_record = PrioritisationRankRecorder(
            0,
            Path("directory2"),
            VariantPrioritisationResultData(
                Path("/path/to/phenopacket-1.json"), VariantData("12", 120434, "A", "G"), 9
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12_120434_A_G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_phenopacket_gene_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            GenePrioritisationResultData(Path("/path/to/phenopacket-2.json"), "GENE", 7),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Gene": "LARGE1",
                        PosixPath("directory1"): 4,
                    }
                },
            ),
        )
        self.add_new_directory_gene_record = PrioritisationRankRecorder(
            0,
            Path("directory2"),
            GenePrioritisationResultData(Path("/path/to/phenopacket-1.json"), "LARGE1", 1),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Gene": "LARGE1",
                        PosixPath("directory1"): 4,
                    }
                },
            ),
        )

    def test_record_rank(self):
        self.assertEqual(len(self.add_new_phenopacket_variant_record.run_comparison), 1)
        self.add_new_phenopacket_variant_record.record_rank()
        self.assertEqual(len(self.add_new_phenopacket_variant_record.run_comparison), 2)
        for i in list(self.add_new_directory_variant_record.run_comparison.values()):
            self.assertFalse(Path("directory2") in i)
        self.assertEqual(len(self.add_new_directory_variant_record.run_comparison), 1)
        self.add_new_directory_variant_record.record_rank()
        self.assertEqual(len(self.add_new_directory_variant_record.run_comparison), 1)
        for i in list(self.add_new_directory_variant_record.run_comparison.values()):
            self.assertTrue(Path("directory2") in i)
        self.assertEqual(len(self.add_new_phenopacket_gene_record.run_comparison), 1)
        self.add_new_phenopacket_gene_record.record_rank()
        self.assertEqual(len(self.add_new_phenopacket_gene_record.run_comparison), 2)
        for i in list(self.add_new_directory_gene_record.run_comparison.values()):
            self.assertFalse(Path("directory2") in i)
        self.add_new_directory_gene_record.record_rank()
        self.assertEqual(len(self.add_new_directory_gene_record.run_comparison), 1)
        for i in list(self.add_new_directory_gene_record.run_comparison.values()):
            self.assertTrue(Path("directory2") in i)


class TestRankStats(unittest.TestCase):
    def setUp(self) -> None:
        self.rank_stats = RankStats()

    def test_add_rank(self):
        self.rank_stats.add_rank(1)
        self.rank_stats.add_rank(3)
        self.rank_stats.add_rank(5)
        self.rank_stats.add_rank(7)
        self.rank_stats.add_rank(10)
        self.assertTrue(
            self.rank_stats.top == 1 and self.rank_stats.top3 == 2 and self.rank_stats.top5 == 3,
            self.rank_stats.total == 5
            and self.rank_stats.found == 5
            and len(self.rank_stats.reciprocal_ranks) == 5,
        )

    def test_percentage_rank(self):
        self.rank_stats.found = 10
        self.assertTrue(self.rank_stats.percentage_rank(3) == 30)

    def test_percentage_top(self):
        self.rank_stats.top, self.rank_stats.found = 10, 20
        self.assertEqual(self.rank_stats.percentage_top(), 50)

    def test_percentage_top3(self):
        self.rank_stats.top3, self.rank_stats.found = 30, 50
        self.assertEqual(self.rank_stats.percentage_top3(), 60)

    def test_percentage_top5(self):
        self.rank_stats.top5, self.rank_stats.found = 70, 160
        self.assertEqual(self.rank_stats.percentage_top5(), 43.75)

    def test_percentage_found(self):
        self.rank_stats.found, self.rank_stats.total = 100, 125
        self.assertEqual(self.rank_stats.percentage_found(), 80)

    def test_mean_reciprocal_rank(self):
        self.rank_stats.reciprocal_ranks = [0.2, 0.4, 0.5, 0.6, 0.8]
        self.assertEqual(self.rank_stats.mean_reciprocal_rank(), 0.5)
