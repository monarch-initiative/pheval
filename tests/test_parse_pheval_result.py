from unittest import TestCase

from pheval.analyse.parse_pheval_result import parse_pheval_result
from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)


class TestParsePhEvalResult(TestCase):
    def test_parse_pheval_gene_result(self):
        self.assertEqual(
            parse_pheval_result(
                RankedPhEvalGeneResult,
                [
                    {
                        "gene_symbol": "PLXNA1",
                        "gene_identifier": "ENSG00000114554",
                        "score": 0.8764,
                        "rank": 1,
                    },
                    {
                        "gene_symbol": "ZNF804B",
                        "gene_identifier": "ENSG00000182348",
                        "score": 0.5777,
                        "rank": 2,
                    },
                    {
                        "gene_symbol": "SMCO2",
                        "gene_identifier": "ENSG00000165935",
                        "score": 0.5777,
                        "rank": 2,
                    },
                    {
                        "gene_symbol": "SPNS1",
                        "gene_identifier": "ENSG00000169682",
                        "score": 0.3765,
                        "rank": 4,
                    },
                ],
            ),
            [
                RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                RankedPhEvalGeneResult(
                    gene_symbol="ZNF804B",
                    gene_identifier="ENSG00000182348",
                    score=0.5777,
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    gene_symbol="SMCO2",
                    gene_identifier="ENSG00000165935",
                    score=0.5777,
                    rank=2,
                ),
                RankedPhEvalGeneResult(
                    gene_symbol="SPNS1",
                    gene_identifier="ENSG00000169682",
                    score=0.3765,
                    rank=4,
                ),
            ],
        )

    def test_parse_pheval_variant_result(self):
        self.assertEqual(
            parse_pheval_result(
                RankedPhEvalVariantResult,
                [
                    {
                        "chromosome": "3",
                        "start": 126730873,
                        "end": 126730873,
                        "ref": "G",
                        "alt": "A",
                        "score": 0.0484,
                        "rank": 1,
                    },
                    {
                        "chromosome": "3",
                        "start": 126730873,
                        "end": 126730873,
                        "ref": "G",
                        "alt": "A",
                        "score": 0.0484,
                        "rank": 1,
                    },
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "score": 0.0484,
                        "rank": 1,
                    },
                ],
            ),
            [
                RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126730873,
                    end=126730873,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126730873,
                    end=126730873,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
            ],
        )

    def test_parse_pheval_disease_result(self):
        self.assertEqual(
            parse_pheval_result(
                RankedPhEvalDiseaseResult,
                [
                    {
                        "disease_name": "Glutaric aciduria type 1",
                        "disease_identifier": "OMIM:231670",
                        "score": 1.0,
                        "rank": 1,
                    },
                    {
                        "disease_name": "Glutaric aciduria type 2",
                        "disease_identifier": "OMIM:231680",
                        "score": 0.8,
                        "rank": 2,
                    },
                    {
                        "disease_name": "Glutaric aciduria type 3",
                        "disease_identifier": "OMIM:231690",
                        "score": 0.6,
                        "rank": 3,
                    },
                ],
            ),
            [
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=1.0,
                    rank=1,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 2",
                    disease_identifier="OMIM:231680",
                    score=0.8,
                    rank=2,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 3",
                    disease_identifier="OMIM:231690",
                    score=0.6,
                    rank=3,
                ),
            ],
        )
