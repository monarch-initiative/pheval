import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl

from pheval.post_processing.mondo_mapping import parse_mondo_mapping_table
from pheval.post_processing.phenopacket_truth_set import PhenopacketTruthSet
from pheval.utils.phenopacket_utils import GenomicVariant, ProbandCausativeGene, ProbandDisease

diagnosed_genes = [
    ProbandCausativeGene(gene_symbol="FGD1", gene_identifier="ENSG00000102302"),
    ProbandCausativeGene(gene_symbol="RTTN", gene_identifier="ENSG00000176225"),
]

variants = [
    GenomicVariant(chrom="X", pos=54492285, ref="C", alt="T"),
    GenomicVariant(chrom="18", pos=67691994, ref="G", alt="A"),
]

diseases = [ProbandDisease(disease_name="Cystic Fibrosis", disease_identifier="OMIM:219700")]

mondo_mapping_table = parse_mondo_mapping_table()


class TestPhenopacketTruthSet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.phenopacket_truth_set = PhenopacketTruthSet(Path("/path/to/phenopacket_dir"))
        cls.phenopacket_truth_set._get_causative_genes = MagicMock(return_value=diagnosed_genes)
        cls.phenopacket_truth_set._get_causative_variants = MagicMock(return_value=variants)
        cls.phenopacket_truth_set._get_causative_diseases = MagicMock(return_value=diseases)
        cls.mock_gene_classified_results = pl.DataFrame(
            [
                {
                    "gene_symbol": "RTTN",
                    "gene_identifier": "ENSG00000176225",
                    "score": 0.0,
                    "rank": 0,
                    "true_positive": True,
                },
                {
                    "gene_symbol": "FGD1",
                    "gene_identifier": "ENSG00000102302",
                    "score": 0.0,
                    "rank": 0,
                    "true_positive": True,
                },
            ]
        )
        cls.mock_gene_ranked_results = pl.DataFrame(
            {
                "gene_symbol": ["FGD1", "BRCA1"],
                "gene_identifier": ["NCBIGene:2245", "ENSG00000012048"],
                "score": [0.9, 0.8],
                "rank": [1, 2],
            }
        )
        cls.mock_variant_classified_results = pl.DataFrame(
            [
                {
                    "chrom": "18",
                    "start": 67691994,
                    "end": 67691994,
                    "ref": "G",
                    "alt": "A",
                    "variant_id": "18-67691994-G-A",
                    "score": 0.0,
                    "rank": 0,
                    "true_positive": True,
                },
                {
                    "chrom": "X",
                    "start": 54492285,
                    "end": 54492285,
                    "ref": "C",
                    "alt": "T",
                    "variant_id": "X-54492285-C-T",
                    "score": 0.0,
                    "rank": 0,
                    "true_positive": True,
                },
            ]
        )
        cls.mock_variant_ranked_results = pl.DataFrame(
            [
                {
                    "chrom": "X",
                    "start": 54492285,
                    "end": 54492285,
                    "ref": "C",
                    "alt": "T",
                    "variant_id": "X-54492285-C-T",
                    "score": 1.0,
                    "rank": 1,
                },
                {
                    "chrom": "19",
                    "start": 67691994,
                    "end": 67691994,
                    "ref": "G",
                    "alt": "A",
                    "variant_id": "18-67691994-G-A",
                    "score": 0.9,
                    "rank": 2,
                },
            ]
        )
        cls.mock_disease_classified_results = pl.DataFrame(
            [
                {
                    "disease_identifier": "OMIM:219700",
                    "score": 0.0,
                    "rank": 0,
                    "true_positive": True,
                    "mondo_identifier": "MONDO:0009061",
                }
            ]
        )
        cls.mock_disease_ranked_results = pl.DataFrame(
            {
                "disease_identifier": ["OMIM:12345", "OMIM:6789"],
                "score": [0.9, 0.8],
                "rank": [1, 2],
            }
        )

    def test_classified_gene(self):
        self.assertTrue(
            self.phenopacket_truth_set.classified_gene("dummy_result_name")
            .sort("gene_identifier")
            .equals(self.mock_gene_classified_results.sort("gene_identifier"))
        )

    @patch("polars.read_parquet")
    def test_merge_gene_results(self, mock_read_parquet):
        mock_read_parquet.return_value = self.mock_gene_classified_results
        self.assertTrue(
            self.phenopacket_truth_set.merge_gene_results(
                self.mock_gene_ranked_results, "output_file"
            )
            .sort("gene_identifier")
            .equals(
                pl.DataFrame(
                    [
                        {
                            "gene_symbol": "FGD1",
                            "gene_identifier": "NCBIGene:2245",
                            "score": 0.9,
                            "rank": 1,
                            "true_positive": True,
                        },
                        {
                            "gene_symbol": "BRCA1",
                            "gene_identifier": "ENSG00000012048",
                            "score": 0.8,
                            "rank": 2,
                            "true_positive": False,
                        },
                        {
                            "gene_symbol": "RTTN",
                            "gene_identifier": "ENSG00000176225",
                            "score": 0.0,
                            "rank": 0,
                            "true_positive": True,
                        },
                    ]
                ).sort("gene_identifier")
            )
        )

    def test_classified_variant(self):
        self.assertTrue(
            self.phenopacket_truth_set.classified_variant("dummy_result_name")
            .sort("chrom")
            .equals(self.mock_variant_classified_results.sort("chrom"))
        )

    @patch("polars.read_parquet")
    def test_merge_variant_results(self, mock_read_parquet):
        mock_read_parquet.return_value = self.mock_variant_classified_results
        self.assertTrue(
            self.phenopacket_truth_set.merge_variant_results(
                self.mock_variant_ranked_results, "output_file"
            )
            .sort("chrom")
            .equals(
                pl.DataFrame(
                    [
                        {
                            "chrom": "18",
                            "start": 67691994,
                            "end": 67691994,
                            "ref": "G",
                            "alt": "A",
                            "variant_id": "18-67691994-G-A",
                            "score": 0.0,
                            "rank": 0,
                            "true_positive": True,
                        },
                        {
                            "chrom": "19",
                            "start": 67691994,
                            "end": 67691994,
                            "ref": "G",
                            "alt": "A",
                            "variant_id": "18-67691994-G-A",
                            "score": 0.9,
                            "rank": 2,
                            "true_positive": False,
                        },
                        {
                            "chrom": "X",
                            "start": 54492285,
                            "end": 54492285,
                            "ref": "C",
                            "alt": "T",
                            "variant_id": "X-54492285-C-T",
                            "score": 1.0,
                            "rank": 1,
                            "true_positive": True,
                        },
                    ]
                ).sort("chrom")
            )
        )

    def test_classified_disease(self):
        self.assertTrue(
            self.phenopacket_truth_set.classified_disease("dummy_result_name", mondo_mapping_table)
            .sort("disease_identifier")
            .equals(self.mock_disease_classified_results.sort("disease_identifier"))
        )

    @patch("polars.read_parquet")
    def test_merge_disease_results(self, mock_read_parquet):
        mock_read_parquet.return_value = self.mock_disease_classified_results
        print(
            self.phenopacket_truth_set.merge_disease_results(
                self.mock_disease_ranked_results, "output_file", mondo_mapping_table
            ).sort("disease_identifier")
        )
        self.assertTrue(
            self.phenopacket_truth_set.merge_disease_results(
                self.mock_disease_ranked_results, "output_file", mondo_mapping_table
            )
            .sort("disease_identifier")
            .equals(
                pl.DataFrame(
                    [
                        {
                            "disease_identifier": "OMIM:12345",
                            "score": 0.9,
                            "rank": 1,
                            "true_positive": False,
                            "mondo_identifier": "OMIM:12345",
                        },
                        {
                            "disease_identifier": "OMIM:6789",
                            "score": 0.8,
                            "rank": 2,
                            "true_positive": False,
                            "mondo_identifier": "OMIM:6789",
                        },
                        {
                            "disease_identifier": "OMIM:219700",
                            "score": 0.0,
                            "rank": 0,
                            "true_positive": True,
                            "mondo_identifier": "MONDO:0009061",
                        },
                    ]
                ).sort("disease_identifier")
            )
        )
