import unittest

import polars as pl

from pheval.post_processing.phenopacket_truth_set import calculate_end_pos
from pheval.post_processing.post_processing import SortOrder, _rank_results

gene_results = pl.DataFrame(
    [
        {"gene_symbol": "A4GNT", "gene_identifier": "ENSG00000118017", "score": 0.6529},
        {"gene_symbol": "MAP3K14", "gene_identifier": "ENSG00000006062", "score": 0.9234},
        {"gene_symbol": "OR14J1", "gene_identifier": "ENSG00000204695", "score": 0.6529},
        {"gene_symbol": "PAGE1", "gene_identifier": "ENSG00000068985", "score": 0.5235},
    ]
)

variant_results = pl.DataFrame(
    [
        {
            "chromosome": "5",
            "start": 23457444233,
            "end": 23457444234,
            "ref": "A",
            "alt": "C",
            "score": 0.9348,
            "grouping_id": "ABYEIF",
        },
        {
            "chromosome": "12",
            "start": 12754332,
            "end": 12754333,
            "ref": "T",
            "alt": "G",
            "score": 0.9999,
            "grouping_id": "JBDFUIBWIX",
        },
        {
            "chromosome": "12",
            "start": 12754332,
            "end": 12754333,
            "ref": "C",
            "alt": "A",
            "score": 0.9999,
            "grouping_id": "JBDFUIBWIX",
        },
        {
            "chromosome": "X",
            "start": 93473023,
            "end": 93473024,
            "ref": "A",
            "alt": "G",
            "score": 0.1245,
            "grouping_id": "DCBOASX",
        },
        {
            "chromosome": "8",
            "start": 532356,
            "end": 532357,
            "ref": "A",
            "alt": "C",
            "score": 0.4578,
            "grouping_id": "WIUDJCNOP",
        },
        {
            "chromosome": "8",
            "start": 532356,
            "end": 532357,
            "ref": "G",
            "alt": "T",
            "score": 0.4578,
            "grouping_id": "OIDJCXQP",
        },
    ]
)

disease_results = pl.DataFrame(
    [
        {"disease_identifier": "OMIM:604131", "mondo_identifier": "MONDO:0015264", "score": 0.9123},
        {"disease_identifier": "ORPHA:84", "mondo_identifier": "MONDO:0009825", "score": 0.8021},
        {
            "disease_identifier": "ORPHA:604131",
            "mondo_identifier": "MONDO:0015264",
            "score": 0.9123,
        },
        {"disease_identifier": "DOID:14330", "mondo_identifier": "MONDO:0007254", "score": 0.6734},
        {"disease_identifier": "OMIM:114480", "mondo_identifier": "MONDO:0011783", "score": 0.4312},
    ]
)


class TestCalculateEndPos(unittest.TestCase):
    def test_calculate_end_pos_single_ref_length(self):
        self.assertEqual(calculate_end_pos(variant_start=13007113, variant_ref="G"), 13007113)

    def test_calculate_end_pos(self):
        self.assertEqual(calculate_end_pos(variant_start=114269997, variant_ref="ACAG"), 114270000)


class TestRankResults(unittest.TestCase):
    def test__rank_results(self):
        self.assertTrue(
            _rank_results(gene_results, SortOrder.DESCENDING).equals(
                pl.DataFrame(
                    [
                        {
                            "gene_symbol": "MAP3K14",
                            "gene_identifier": "ENSG00000006062",
                            "score": 0.9234,
                            "rank": 1,
                        },
                        {
                            "gene_symbol": "A4GNT",
                            "gene_identifier": "ENSG00000118017",
                            "score": 0.6529,
                            "rank": 3,
                        },
                        {
                            "gene_symbol": "OR14J1",
                            "gene_identifier": "ENSG00000204695",
                            "score": 0.6529,
                            "rank": 3,
                        },
                        {
                            "gene_symbol": "PAGE1",
                            "gene_identifier": "ENSG00000068985",
                            "score": 0.5235,
                            "rank": 4,
                        },
                    ]
                )
            )
        )

    def test__rank_results_grouping_id(self):
        self.assertTrue(
            _rank_results(variant_results, SortOrder.DESCENDING).equals(
                pl.DataFrame(
                    [
                        {
                            "chromosome": "12",
                            "start": 12754332,
                            "end": 12754333,
                            "ref": "T",
                            "alt": "G",
                            "score": 0.9999,
                            "grouping_id": "JBDFUIBWIX",
                            "min_rank": 1,
                            "rank": 1,
                        },
                        {
                            "chromosome": "12",
                            "start": 12754332,
                            "end": 12754333,
                            "ref": "C",
                            "alt": "A",
                            "score": 0.9999,
                            "grouping_id": "JBDFUIBWIX",
                            "min_rank": 1,
                            "rank": 1,
                        },
                        {
                            "chromosome": "5",
                            "start": 23457444233,
                            "end": 23457444234,
                            "ref": "A",
                            "alt": "C",
                            "score": 0.9348,
                            "grouping_id": "ABYEIF",
                            "min_rank": 2,
                            "rank": 2,
                        },
                        {
                            "chromosome": "8",
                            "start": 532356,
                            "end": 532357,
                            "ref": "A",
                            "alt": "C",
                            "score": 0.4578,
                            "grouping_id": "WIUDJCNOP",
                            "min_rank": 3,
                            "rank": 4,
                        },
                        {
                            "chromosome": "8",
                            "start": 532356,
                            "end": 532357,
                            "ref": "G",
                            "alt": "T",
                            "score": 0.4578,
                            "grouping_id": "OIDJCXQP",
                            "min_rank": 4,
                            "rank": 4,
                        },
                        {
                            "chromosome": "X",
                            "start": 93473023,
                            "end": 93473024,
                            "ref": "A",
                            "alt": "G",
                            "score": 0.1245,
                            "grouping_id": "DCBOASX",
                            "min_rank": 5,
                            "rank": 5,
                        },
                    ]
                )
            )
        )

    def test__rank_results_mondo_id(self):
        self.assertTrue(
            _rank_results(disease_results, SortOrder.DESCENDING).equals(
                pl.DataFrame(
                    [
                        {
                            "disease_identifier": "OMIM:604131",
                            "mondo_identifier": "MONDO:0015264",
                            "score": 0.9123,
                            "min_rank": 1,
                            "rank": 1,
                        },
                        {
                            "disease_identifier": "ORPHA:604131",
                            "mondo_identifier": "MONDO:0015264",
                            "score": 0.9123,
                            "min_rank": 1,
                            "rank": 1,
                        },
                        {
                            "disease_identifier": "ORPHA:84",
                            "mondo_identifier": "MONDO:0009825",
                            "score": 0.8021,
                            "min_rank": 2,
                            "rank": 2,
                        },
                        {
                            "disease_identifier": "DOID:14330",
                            "mondo_identifier": "MONDO:0007254",
                            "score": 0.6734,
                            "min_rank": 3,
                            "rank": 3,
                        },
                        {
                            "disease_identifier": "OMIM:114480",
                            "mondo_identifier": "MONDO:0011783",
                            "score": 0.4312,
                            "min_rank": 4,
                            "rank": 4,
                        },
                    ]
                )
            )
        )
