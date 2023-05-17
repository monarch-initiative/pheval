import unittest

from pheval.post_processing.post_processing import (
    PhEvalGeneResult,
    PhEvalVariantResult, PhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult, RankedPhEvalDiseaseResult,
    ResultSorter,
    ScoreRanker,
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
    PhEvalDiseaseResult(disease_name="Glutaric acidemia I", disease_identifier="OMIM:231670", score=4.284),
    PhEvalDiseaseResult(disease_name="Diencephalic-mesencephalic junction dysplasia syndrome 2",
                        disease_identifier="OMIM:618646", score=-1.378),
    PhEvalDiseaseResult(disease_name=" Brain small vessel disease 2", disease_identifier="OMIM:614483", score=-1.871)
]


class TestRankedPhEvalGeneResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_gene_result = RankedPhEvalGeneResult(
            pheval_gene_result=PhEvalGeneResult(
                gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
            ),
            rank=1,
        )

    def test_as_dict(self):
        self.assertEqual(
            self.pheval_gene_result.as_dict(),
            {
                "gene_symbol": "A4GNT",
                "gene_identifier": "ENSG00000118017",
                "score": 0.6529,
                "rank": 1,
            },
        )


class TestRankedPhEvalVariantResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_variant_result = RankedPhEvalVariantResult(
            pheval_variant_result=PhEvalVariantResult(
                chromosome="12", start=12754332, end=12754333, ref="T", alt="G", score=0.9999
            ),
            rank=3,
        )

    def test_as_dict(self):
        self.assertEqual(
            self.pheval_variant_result.as_dict(),
            {
                "chromosome": "12",
                "start": 12754332,
                "end": 12754333,
                "ref": "T",
                "alt": "G",
                "score": 0.9999,
                "rank": 3,
            },
        )


class TestRankedPhEvalDiseaseResult(unittest.TestCase):
    def setUp(self) -> None:
        self.pheval_disease_result = RankedPhEvalDiseaseResult(
            pheval_disease_result=PhEvalDiseaseResult(disease_name="Glutaric acidemia I",
                                                      disease_identifier="OMIM:231670", score=4.284),
            rank=1,
        )

    def test_as_dict(self):
        self.assertEqual(
            self.pheval_disease_result.as_dict(),
            {
                "disease_name": "Glutaric acidemia I",
                "disease_identifier": "OMIM:231670",
                "score": 4.284,
                "rank": 1,
            },
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


class TestScoreRanker(unittest.TestCase):
    def setUp(self) -> None:
        self.score_ranker_descending = ScoreRanker(SortOrder.DESCENDING)
        self.score_ranker_ascending = ScoreRanker(SortOrder.ASCENDING)

    def test_check_rank_order_descending(self):
        self.score_ranker_descending.rank_scores(0.9)
        self.score_ranker_descending.rank_scores(0.8)
        with self.assertRaises(ValueError):
            self.score_ranker_descending.rank_scores(0.9)

    def test_check_rank_order_ascending(self):
        self.score_ranker_ascending.rank_scores(0.1)
        self.score_ranker_ascending.rank_scores(0.2)
        with self.assertRaises(ValueError):
            self.score_ranker_ascending.rank_scores(0.1)

    def test_rank_scores_first_rank(self):
        self.assertEqual(self.score_ranker_descending.rank_scores(0.7342), 1)

    def test_rank_scores_increase_rank(self):
        self.assertEqual(self.score_ranker_descending.rank_scores(0.7342), 1)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.3452), 2)

    def test_rank_scores_same_rank(self):
        self.assertEqual(self.score_ranker_descending.rank_scores(0.7342), 1)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.3452), 2)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.3452), 2)

    def test_rank_scores_count_increase(self):
        self.assertEqual(self.score_ranker_descending.rank_scores(0.7342), 1)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.3452), 2)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.3452), 2)
        self.assertEqual(self.score_ranker_descending.rank_scores(0.1234), 4)


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
            _rank_pheval_result(self.sorted_gene_result, SortOrder.DESCENDING),
            [
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
                    ),
                    rank=1,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
                    ),
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529
                    ),
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235
                    ),
                    rank=4,
                ),
            ],
        )

    def test_rank_pheval_results_variant(self):
        self.assertEqual(
            _rank_pheval_result(self.sorted_variant_result, SortOrder.ASCENDING),
            [
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
                    ),
                    rank=1,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578
                    ),
                    rank=2,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="5",
                        start=23457444233,
                        end=23457444234,
                        ref="A",
                        alt="C",
                        score=0.9348,
                    ),
                    rank=3,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="12",
                        start=12754332,
                        end=12754333,
                        ref="T",
                        alt="G",
                        score=0.9999,
                    ),
                    rank=4,
                ),
            ],
        )

    def test_rank_pheval_results_disease(self):
        self.assertEqual(_rank_pheval_result(self.sorted_disease_result, SortOrder.DESCENDING),
                         [RankedPhEvalDiseaseResult(
                             pheval_disease_result=PhEvalDiseaseResult(disease_name='Glutaric acidemia I',
                                                                       disease_identifier='OMIM:231670', score=4.284),
                             rank=1),
                             RankedPhEvalDiseaseResult(
                                 pheval_disease_result=PhEvalDiseaseResult(
                                     disease_name='Diencephalic-mesencephalic junction dysplasia syndrome 2',
                                     disease_identifier='OMIM:618646',
                                     score=-1.378), rank=2),
                             RankedPhEvalDiseaseResult(
                                 pheval_disease_result=PhEvalDiseaseResult(
                                     disease_name=' Brain small vessel disease 2',
                                     disease_identifier='OMIM:614483',
                                     score=-1.871), rank=3)]
                         )


class TestCreatePhEvalResult(unittest.TestCase):
    def test_create_pheval_result_gene(self):
        self.assertEqual(
            _create_pheval_result(pheval_gene_result, "descending"),
            [
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="MAP3K14", gene_identifier="ENSG00000006062", score=0.9234
                    ),
                    rank=1,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="A4GNT", gene_identifier="ENSG00000118017", score=0.6529
                    ),
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="OR14J1", gene_identifier="ENSG00000204695", score=0.6529
                    ),
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    pheval_gene_result=PhEvalGeneResult(
                        gene_symbol="PAGE1", gene_identifier="ENSG00000068985", score=0.5235
                    ),
                    rank=4,
                ),
            ],
        )

    def test_create_pheval_result_variant(self):
        self.assertEqual(
            _create_pheval_result(pheval_variant_result, "ascending"),
            [
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="X", start=93473023, end=93473024, ref="A", alt="G", score=0.1245
                    ),
                    rank=1,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="8", start=532356, end=532357, ref="A", alt="C", score=0.4578
                    ),
                    rank=2,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="5",
                        start=23457444233,
                        end=23457444234,
                        ref="A",
                        alt="C",
                        score=0.9348,
                    ),
                    rank=3,
                ),
                RankedPhEvalVariantResult(
                    pheval_variant_result=PhEvalVariantResult(
                        chromosome="12",
                        start=12754332,
                        end=12754333,
                        ref="T",
                        alt="G",
                        score=0.9999,
                    ),
                    rank=4,
                ),
            ],
        )

    def test_create_pheval_result_disease(self):
        self.assertEqual(_create_pheval_result(pheval_disease_result, "descending"), [RankedPhEvalDiseaseResult(
            pheval_disease_result=PhEvalDiseaseResult(disease_name='Glutaric acidemia I',
                                                      disease_identifier='OMIM:231670', score=4.284), rank=1),
            RankedPhEvalDiseaseResult(
                pheval_disease_result=PhEvalDiseaseResult(
                    disease_name='Diencephalic-mesencephalic junction dysplasia syndrome 2',
                    disease_identifier='OMIM:618646',
                    score=-1.378), rank=2),
            RankedPhEvalDiseaseResult(
                pheval_disease_result=PhEvalDiseaseResult(
                    disease_name=' Brain small vessel disease 2',
                    disease_identifier='OMIM:614483',
                    score=-1.871), rank=3)]
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
