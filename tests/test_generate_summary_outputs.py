import unittest
from collections import defaultdict

import pandas as pd

from pheval.analyse.generate_summary_outputs import RankComparisonGenerator, merge_results


class TestMergeResults(unittest.TestCase):
    def setUp(self) -> None:
        self.result_1 = {
            1: {
                "Phenopacket": "phenopacket1.json",
                "Gene": "GCDH",
                "/path/to/results_directory1": 1,
            }
        }
        self.result_2 = {
            1: {
                "Phenopacket": "phenopacket1.json",
                "Gene": "GCDH",
                "/path/to/results_directory2": 5,
            }
        }

    def test_merge_results(self):
        self.assertEqual(
            merge_results(self.result_1, self.result_2),
            {
                1: {
                    "Phenopacket": "phenopacket1.json",
                    "Gene": "GCDH",
                    "/path/to/results_directory1": 1,
                    "/path/to/results_directory2": 5,
                }
            },
        )


class TestRankComparisonGenerator(unittest.TestCase):
    def setUp(self) -> None:
        self.gene_rank_comparisons = RankComparisonGenerator(
            defaultdict(
                dict,
                {
                    1: {
                        "Phenopacket": "phenopacket1.json",
                        "Gene": "GCDH",
                        "/path/to/results_directory1": 1,
                        "/path/to/results_directory2": 5,
                    }
                },
            )
        )
        self.variant_rank_comparisons = RankComparisonGenerator(
            defaultdict(
                dict,
                {
                    1: {
                        "Phenopacket": "phenopacket1.json",
                        "Variant": "3-12563453454-C-T",
                        "/path/to/results_directory1": 9,
                        "/path/to/results_directory2": 3,
                    }
                },
            )
        )

    def test_generate_gene_dataframe(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Gene": "GCDH",
                    "/path/to/results_directory1": 1,
                    "/path/to/results_directory2": 5,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.gene_rank_comparisons._generate_dataframe()))

    def test_generate_variant_dataframe(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Variant": "3-12563453454-C-T",
                    "/path/to/results_directory1": 9,
                    "/path/to/results_directory2": 3,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.variant_rank_comparisons._generate_dataframe()))

    def test_calculate_gene_rank_difference(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Gene": "GCDH",
                    "/path/to/results_directory1": 1,
                    "/path/to/results_directory2": 5,
                    "rank_change": -4,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.gene_rank_comparisons._calculate_rank_difference()))

    def test_calculate_variant_rank_difference(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Variant": "3-12563453454-C-T",
                    "/path/to/results_directory1": 9,
                    "/path/to/results_directory2": 3,
                    "rank_change": 6,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.variant_rank_comparisons._calculate_rank_difference()))
