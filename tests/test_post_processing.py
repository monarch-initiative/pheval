import unittest

import pandas as pd

from pheval.post_processing.post_processing import (
    PhEvalDiseaseResult,
    PhEvalGeneResult,
    PhEvalVariantResult,
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
    ResultRanker,
    ResultSorter,
    SortOrder,
    _create_pheval_result,
    _return_sort_order,
    calculate_end_pos,
)

pheval_gene_result = [
    PhEvalGeneResult(gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529),
    PhEvalGeneResult(gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234),
    PhEvalGeneResult(gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529),
    PhEvalGeneResult(gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235),
]
pheval_variant_result = [
    PhEvalVariantResult(
        chromosome="5",
        start=23457444233,
        end=23457444234,
        ref="A",
        alt="C",
        score=0.9348,
    ),
    PhEvalVariantResult(
        chromosome="12", start=12754332, end=12754333, ref="T", alt="G", score=0.9999
    ),
    PhEvalVariantResult(
        chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
    ),
    PhEvalVariantResult(chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578),
]
pheval_disease_result = [
    PhEvalDiseaseResult(
        disease_name="Glutaric acidemia I", disease_identifier="OMIM:231670", score=4.284
    ),
    PhEvalDiseaseResult(
        disease_name="Diencephalic-mesencephalic junction dysplasia syndrome 2",
        disease_identifier="OMIM:618646",
        score=4.284,
    ),
    PhEvalDiseaseResult(
        disease_name=" Brain small vessel disease 2", disease_identifier="OMIM:614483", score=-1.871
    ),
]


class TestRankedPhEvalGeneResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_gene_result = PhEvalGeneResult(
            gene_symbol="A4GNT",
            gene_identifier="ENSG00000118017",
            score=0.6529,
        )

    def test_from_gene_result(self):
        self.assertEqual(
            RankedPhEvalGeneResult.from_gene_result(self.pheval_gene_result, 1),
            RankedPhEvalGeneResult(
                gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529, rank=1
            ),
        )


class TestRankedPhEvalVariantResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_variant_result = PhEvalVariantResult(
            chromosome="12",
            start=12754332,
            end=12754333,
            ref="T",
            alt="G",
            score=0.9999,
        )

    def test_from_variant_result(self):
        self.assertEqual(
            RankedPhEvalVariantResult.from_variant_result(self.pheval_variant_result, 3),
            RankedPhEvalVariantResult(
                chromosome="12",
                start=12754332,
                end=12754333,
                ref="T",
                alt="G",
                score=0.9999,
                rank=3,
            ),
        )


class TestRankedPhEvalDiseaseResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_disease_result = PhEvalDiseaseResult(
            disease_name="Bardet-Biedl Syndrome",
            disease_identifier="OMIM:209900",
            score=0.9999,
        )

    def test_from_disease_result(self):
        self.assertEqual(
            RankedPhEvalDiseaseResult.from_disease_result(self.pheval_disease_result, 1),
            RankedPhEvalDiseaseResult(
                disease_name="Bardet-Biedl Syndrome",
                disease_identifier="OMIM:209900",
                score=0.9999,
                rank=1,
            ),
        )


class TestResultSorter(unittest.TestCase):
    def setUp(self) -> None:
        self.gene_results = ResultSorter(
            pheval_results=pheval_gene_result,
            sort_order=SortOrder.DESCENDING,
        )
        self.variant_results = ResultSorter(
            pheval_results=pheval_variant_result,
            sort_order=SortOrder.ASCENDING,
        )

    def test_sort_by_decreasing_score(self):
        self.assertEqual(
            self.gene_results._sort_by_decreasing_score(),
            [
                PhEvalGeneResult(
                    gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
                ),
                PhEvalGeneResult(
                    gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235
                ),
            ],
        )

    def test_sort_by_increasing_score(self):
        self.assertEqual(
            self.variant_results._sort_by_increasing_score(),
            [
                PhEvalVariantResult(
                    chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
                ),
                PhEvalVariantResult(
                    chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578
                ),
                PhEvalVariantResult(
                    chromosome="5",
                    start=23457444233,
                    end=23457444234,
                    ref="A",
                    alt="C",
                    score=0.9348,
                ),
                PhEvalVariantResult(
                    chromosome="12", start=12754332, end=12754333, ref="T", alt="G", score=0.9999
                ),
            ],
        )

    def test_sort_pheval_results_not_pvalue(self):
        print(self.gene_results.sort_order)
        self.assertEqual(
            self.gene_results.sort_pheval_results(),
            [
                PhEvalGeneResult(
                    gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
                ),
                PhEvalGeneResult(
                    gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235
                ),
            ],
        )

    def test_sort_pheval_results_pvalue(self):
        self.assertEqual(
            self.variant_results.sort_pheval_results(),
            [
                PhEvalVariantResult(
                    chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
                ),
                PhEvalVariantResult(
                    chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578
                ),
                PhEvalVariantResult(
                    chromosome="5",
                    start=23457444233,
                    end=23457444234,
                    ref="A",
                    alt="C",
                    score=0.9348,
                ),
                PhEvalVariantResult(
                    chromosome="12", start=12754332, end=12754333, ref="T", alt="G", score=0.9999
                ),
            ],
        )


class TestResultRanker(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gene_result_ranker = ResultRanker(
            [
                PhEvalGeneResult(
                    gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
                ),
                PhEvalGeneResult(
                    gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529
                ),
                PhEvalGeneResult(
                    gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235
                ),
            ],
            SortOrder.DESCENDING,
        )
        cls.variant_result_ranker = ResultRanker(
            [
                PhEvalVariantResult(
                    chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
                ),
                PhEvalVariantResult(
                    chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578
                ),
                PhEvalVariantResult(
                    chromosome="5",
                    start=23457444233,
                    end=23457444234,
                    ref="A",
                    alt="C",
                    score=0.9999,
                ),
                PhEvalVariantResult(
                    chromosome="12", start=12754332, end=12754333, ref="T", alt="G", score=0.9999
                ),
            ],
            sort_order=SortOrder.ASCENDING,
        )
        cls.variant_result_ranker_grouping_id = ResultRanker(
            [
                PhEvalVariantResult(
                    chromosome="X",
                    start=93473023,
                    end=93473024,
                    ref="A",
                    alt="G",
                    score=0.1245,
                    grouping_id="4567",
                ),
                PhEvalVariantResult(
                    chromosome="8",
                    start=532356,
                    end=532357,
                    ref="A",
                    alt="C",
                    score=0.4578,
                    grouping_id="789",
                ),
                PhEvalVariantResult(
                    chromosome="5",
                    start=23457444233,
                    end=23457444234,
                    ref="A",
                    alt="C",
                    score=0.9999,
                    grouping_id="12345",
                ),
                PhEvalVariantResult(
                    chromosome="12",
                    start=12754332,
                    end=12754333,
                    ref="T",
                    alt="G",
                    score=0.9999,
                    grouping_id="12345",
                ),
            ],
            sort_order=SortOrder.DESCENDING,
        )
        cls.disease_result_ranker = ResultRanker(
            pheval_disease_result, sort_order=SortOrder.DESCENDING
        )

    def test__has_valid_grouping_id(self):
        df = pd.DataFrame({"score": [0.5, 0.7, 0.3], "grouping_id": ["A", "B", "C"]})
        self.assertTrue(self.variant_result_ranker._has_valid_grouping_id(df))

    def test__has_valid_grouping_id_present_with_none(self):
        df = pd.DataFrame({"score": [0.5, 0.7, 0.3], "grouping_id": ["A", None, "C"]})
        self.assertFalse(self.variant_result_ranker._has_valid_grouping_id(df))

    def test__has_valid_grouping_id_not_present(self):
        df = pd.DataFrame({"score": [0.5, 0.7, 0.3]})
        self.assertFalse(self.variant_result_ranker._has_valid_grouping_id(df))

    def test__rank_with_grouping_id(self):
        df = pd.DataFrame(
            {
                "score": [0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5, 0.4],
                "grouping_id": ["A", "A", "B", "C", "D", "E", "F", "F", "G"],
            }
        )
        self.assertTrue(
            self.gene_result_ranker._rank_with_grouping_id(df).equals(
                pd.DataFrame(
                    {
                        "score": [0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5, 0.4],
                        "grouping_id": ["A", "A", "B", "C", "D", "E", "F", "F", "G"],
                        "min_rank": [1.0, 1.0, 3.0, 5.0, 4.0, 6.0, 7.0, 7.0, 9.0],
                        "rank": [1.0, 1.0, 3.0, 5.0, 5.0, 6.0, 7.0, 7.0, 9.0],
                    }
                )
            )
        )

    def test__rank_without_grouping_id(self):
        df = pd.DataFrame(
            {
                "score": [0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5, 0.4],
            }
        )
        self.assertTrue(
            self.gene_result_ranker._rank_without_grouping_id(df).equals(
                pd.DataFrame(
                    {
                        "score": [0.9, 0.9, 0.8, 0.7, 0.7, 0.6, 0.5, 0.5, 0.4],
                        "rank": [2.0, 2.0, 3.0, 5.0, 5.0, 6.0, 8.0, 8.0, 9.0],
                    }
                )
            )
        )

    def test_rank_pheval_results_gene(self):
        self.assertTrue(
            self.gene_result_ranker.rank().equals(
                pd.DataFrame(
                    {
                        "gene_symbol": ["MAP3K14", "A4GNT", "OR14J1", "PAGE1"],
                        "gene_identifier": [
                            "ENSG00000006062",
                            "ENSG00000118017",
                            "ENSG00000204695",
                            "ENSG00000068985",
                        ],
                        "score": [0.9234, 0.6529, 0.6529, 0.5235],
                        "rank": [1.0, 3.0, 3.0, 4.0],
                    }
                )
            )
        )

    def test_rank_pheval_results_variant(self):
        self.assertTrue(
            self.variant_result_ranker.rank().equals(
                pd.DataFrame(
                    {
                        "chromosome": ["X", "8", "5", "12"],
                        "start": [93473023, 532356, 23457444233, 12754332],
                        "end": [93473024, 532357, 23457444234, 12754333],
                        "ref": ["A", "A", "A", "T"],
                        "alt": ["G", "C", "C", "G"],
                        "score": [0.1245, 0.4578, 0.9999, 0.9999],
                        "rank": [1.0, 2.0, 4.0, 4.0],
                    }
                )
            )
        )

    def test_rank_pheval_results_variant_grouping_id(self):
        self.assertTrue(
            self.variant_result_ranker_grouping_id.rank().equals(
                pd.DataFrame(
                    {
                        "chromosome": ["X", "8", "5", "12"],
                        "start": [93473023, 532356, 23457444233, 12754332],
                        "end": [93473024, 532357, 23457444234, 12754333],
                        "ref": ["A", "A", "A", "T"],
                        "alt": ["G", "C", "C", "G"],
                        "score": [0.1245, 0.4578, 0.9999, 0.9999],
                        "rank": [4.0, 3.0, 1.0, 1.0],
                    }
                )
            )
        )

    def test_rank_pheval_results_disease(self):
        self.assertTrue(
            self.disease_result_ranker.rank().equals(
                pd.DataFrame(
                    {
                        "disease_name": {
                            0: "Glutaric acidemia I",
                            1: "Diencephalic-mesencephalic junction dysplasia syndrome 2",
                            2: " Brain small vessel disease 2",
                        },
                        "disease_identifier": {
                            0: "OMIM:231670",
                            1: "OMIM:618646",
                            2: "OMIM:614483",
                        },
                        "score": {0: 4.284, 1: 4.284, 2: -1.871},
                        "rank": {0: 2.0, 1: 2.0, 2: 3.0},
                    }
                )
            )
        )


class TestCreatePhEvalResult(unittest.TestCase):
    def test_create_pheval_result_gene(self):
        self.assertTrue(
            _create_pheval_result(pheval_gene_result, "descending").equals(
                pd.DataFrame(
                    {
                        "gene_symbol": ["MAP3K14", "A4GNT", "OR14J1", "PAGE1"],
                        "gene_identifier": [
                            "ENSG00000006062",
                            "ENSG00000118017",
                            "ENSG00000204695",
                            "ENSG00000068985",
                        ],
                        "score": [0.9234, 0.6529, 0.6529, 0.5235],
                        "rank": [1.0, 3.0, 3.0, 4.0],
                    }
                )
            )
        )

    def test_create_pheval_result_variant(self):
        self.assertTrue(
            _create_pheval_result(pheval_variant_result, "ascending").equals(
                pd.DataFrame(
                    {
                        "chromosome": ["X", "8", "5", "12"],
                        "start": [93473023, 532356, 23457444233, 12754332],
                        "end": [93473024, 532357, 23457444234, 12754333],
                        "ref": ["A", "A", "A", "T"],
                        "alt": ["G", "C", "C", "G"],
                        "score": [0.1245, 0.4578, 0.9348, 0.9999],
                        "rank": [1.0, 2.0, 3.0, 4.0],
                    }
                )
            )
        )

    def test_create_pheval_result_disease(self):
        self.assertTrue(
            _create_pheval_result(pheval_disease_result, "descending").equals(
                pd.DataFrame(
                    {
                        "disease_name": {
                            0: "Glutaric acidemia I",
                            1: "Diencephalic-mesencephalic junction dysplasia syndrome 2",
                            2: " Brain small vessel disease 2",
                        },
                        "disease_identifier": {
                            0: "OMIM:231670",
                            1: "OMIM:618646",
                            2: "OMIM:614483",
                        },
                        "score": {0: 4.284, 1: 4.284, 2: -1.871},
                        "rank": {0: 2.0, 1: 2.0, 2: 3.0},
                    }
                )
            )
        )


class TestReturnSortOrder(unittest.TestCase):
    def test_return_sort_order(self):
        self.assertEqual(_return_sort_order("ascending"), SortOrder.ASCENDING)
        self.assertEqual(_return_sort_order("descending"), SortOrder.DESCENDING)

    def test_return_sort_order_incompatible(self):
        with self.assertRaises(ValueError):
            _return_sort_order("null")


class TestCalculateEndPos(unittest.TestCase):
    def test_calculate_end_pos_single_ref_length(self):
        self.assertEqual(calculate_end_pos(variant_start=13007113, variant_ref="G"), 13007113)

    def test_calculate_end_pos(self):
        self.assertEqual(calculate_end_pos(variant_start=114269997, variant_ref="ACAG"), 114270000)
