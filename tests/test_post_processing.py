import unittest

import pandas as pd

from pheval.post_processing.post_processing import (
    PhEvalDiseaseResult,
    PhEvalGeneResult,
    PhEvalVariantResult,
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
    ResultSorter,
    SortOrder,
    _create_pheval_result,
    _rank_pheval_result,
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


class TestRankPhEvalResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.sorted_gene_result = [
            PhEvalGeneResult(
                gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
            ),
            PhEvalGeneResult(gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529),
            PhEvalGeneResult(gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529),
            PhEvalGeneResult(gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235),
        ]
        cls.sorted_variant_result = [
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
        ]
        cls.sorted_disease_result = pheval_disease_result

    def test_rank_pheval_results_gene(self):
        self.assertTrue(
            _rank_pheval_result(self.sorted_gene_result, SortOrder.DESCENDING).equals(
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
            _rank_pheval_result(self.sorted_variant_result, SortOrder.ASCENDING).equals(
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

    def test_rank_pheval_results_disease(self):
        self.assertTrue(
            _rank_pheval_result(self.sorted_disease_result, SortOrder.DESCENDING).equals(
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
