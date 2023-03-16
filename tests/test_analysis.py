import unittest
from collections import defaultdict
from copy import copy
from pathlib import Path, PosixPath

import pandas as pd

from pheval.analyse.analysis import (
    AssessGenePrioritisation,
    AssessVariantPrioritisation,
    GenePrioritisationResultData,
    PrioritisationRankRecorder,
    RankComparisonGenerator,
    RankStats,
    VariantPrioritisationResultData,
    merge_results,
)
from pheval.utils.phenopacket_utils import GenomicVariant, ProbandCausativeGene


class TestPrioritisationRankRecorder(unittest.TestCase):
    def setUp(self) -> None:
        self.add_new_phenopacket_variant_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            VariantPrioritisationResultData(
                Path("/path/to/phenopacket-2.json"), GenomicVariant("1", 4896347, "C", "T"), 9
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
                Path("/path/to/phenopacket-1.json"), GenomicVariant("12", 120434, "A", "G"), 9
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


class TestAssessGenePrioritisation(unittest.TestCase):
    def setUp(self) -> None:
        self.assess_gene_prioritisation = AssessGenePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_gene_results=pd.DataFrame(
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
                ]
            ),
            threshold=0.0,
            ranking_method="combinedScore",
            proband_causative_genes=[
                ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                ProbandCausativeGene(gene_symbol="LARGE1", gene_identifier="ENSG00000133424"),
            ],
        )
        self.assess_gene_prioritisation_pvalue = AssessGenePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_gene_results=pd.DataFrame(
                [
                    {
                        "gene_symbol": "SPNS1",
                        "gene_identifier": "ENSG00000169682",
                        "score": 0.3765,
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
                        "gene_symbol": "PLXNA1",
                        "gene_identifier": "ENSG00000114554",
                        "score": 0.8764,
                        "rank": 4,
                    },
                ]
            ),
            threshold=0.0,
            ranking_method="pValue",
            proband_causative_genes=[
                ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                ProbandCausativeGene(gene_symbol="LARGE1", gene_identifier="ENSG00000133424"),
            ],
        )
        self.gene_rank_stats = RankStats(0, 0, 0, 0)
        self.gene_rank_records = defaultdict(dict)

    def test_record_gene_prioritisation_match(self):
        self.assertEqual(
            self.assess_gene_prioritisation.record_gene_prioritisation_match(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry={
                    "gene_symbol": "PLXNA1",
                    "gene_identifier": "ENSG00000114554",
                    "score": 0.8764,
                    "rank": 1,
                },
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )

    def test_assess_gene_with_pvalue_threshold_fails_cutoff(self):
        assess_pvalue_threshold = copy(self.assess_gene_prioritisation_pvalue)
        assess_pvalue_threshold.threshold = 0.1
        self.assertEqual(
            assess_pvalue_threshold.assess_gene_with_pvalue_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry={
                    "gene_symbol": "PLXNA1",
                    "gene_identifier": "ENSG00000114554",
                    "score": 0.8764,
                    "rank": 1,
                },
                rank_stats=self.gene_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_gene_with_pvalue_threshold_meets_cutoff(self):
        assess_pvalue_threshold = copy(self.assess_gene_prioritisation_pvalue)
        assess_pvalue_threshold.threshold = 0.9
        self.assertEqual(
            assess_pvalue_threshold.assess_gene_with_pvalue_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry={
                    "gene_symbol": "PLXNA1",
                    "gene_identifier": "ENSG00000114554",
                    "score": 0.8764,
                    "rank": 1,
                },
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_gene_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold.assess_gene_with_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry={
                    "gene_symbol": "PLXNA1",
                    "gene_identifier": "ENSG00000114554",
                    "score": 0.8764,
                    "rank": 1,
                },
                rank_stats=self.gene_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_gene_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.5
        self.assertEqual(
            assess_with_threshold.assess_gene_with_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry={
                    "gene_symbol": "PLXNA1",
                    "gene_identifier": "ENSG00000114554",
                    "score": 0.8764,
                    "rank": 1,
                },
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_gene_prioritisation_no_threshold(self):
        self.assess_gene_prioritisation.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.gene_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "PLXNA1",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "LARGE1",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_gene_prioritisation_threshold_fails_pvalue_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.01
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.gene_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "PLXNA1",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "LARGE1",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_gene_prioritisation_threshold_meets_pvalue_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=1, found=1, total=2, reciprocal_ranks=[0.25]),
        )
        self.assertEqual(
            self.gene_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "PLXNA1",
                    Path("/path/to/results_dir"): 4,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "LARGE1",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_gene_prioritisation_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.gene_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "PLXNA1",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "LARGE1",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_gene_prioritisation_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.1
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.gene_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "PLXNA1",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Gene": "LARGE1",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )


class TestAssessVariantPrioritisation(unittest.TestCase):
    def setUp(self) -> None:
        variant_results = pd.DataFrame(
            [
                {
                    "variant": {
                        "chromosome": "3",
                        "start": 126730873,
                        "end": 126730873,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                    },
                    "score": 0.0484,
                    "rank": 1,
                },
                {
                    "variant": {
                        "chromosome": "3",
                        "start": 126730873,
                        "end": 126730873,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                    },
                    "score": 0.0484,
                    "rank": 1,
                },
                {
                    "variant": {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                    },
                    "score": 0.0484,
                    "rank": 1,
                },
            ]
        )
        variant_result_format = variant_results.drop("variant", axis=1).join(
            variant_results.variant.apply(pd.Series)
        )

        self.assess_variant_prioritisation = AssessVariantPrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_variant_results=variant_result_format,
            threshold=0.0,
            ranking_method="combinedScore",
            proband_causative_variants=[
                GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                GenomicVariant(chrom="16", pos=133564345, ref="C", alt="T"),
            ],
        )
        self.assess_variant_prioritisation_pvalue = AssessVariantPrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_variant_results=variant_result_format,
            threshold=0.0,
            ranking_method="pValue",
            proband_causative_variants=[
                GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                GenomicVariant(chrom="16", pos=133564345, ref="C", alt="T"),
            ],
        )
        self.variant_rank_stats = RankStats()
        self.variant_rank_records = defaultdict(dict)

    def test_record_variant_prioritisation_match(self):
        self.assertEqual(
            self.assess_variant_prioritisation.record_variant_prioritisation_match(
                result_entry=pd.Series(
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                        "score": 0.0484,
                        "rank": 1,
                    }
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )

    def test_assess_variant_with_pvalue_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold.assess_variant_with_pvalue_threshold(
                result_entry=pd.Series(
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                        "score": 0.0484,
                        "rank": 1,
                    }
                ),
                rank_stats=self.variant_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_variant_with_pvalue_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold.assess_variant_with_pvalue_threshold(
                result_entry=pd.Series(
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                        "score": 0.0484,
                        "rank": 1,
                    }
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_variant_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold.assess_variant_with_threshold(
                result_entry=pd.Series(
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                        "score": 0.0484,
                        "rank": 1,
                    }
                ),
                rank_stats=self.variant_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_variant_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold.assess_variant_with_threshold(
                result_entry=pd.Series(
                    {
                        "chromosome": "3",
                        "start": 126741108,
                        "end": 126741108,
                        "ref": "G",
                        "alt": "A",
                        "gene": "PLXNA1",
                        "score": 0.0484,
                        "rank": 1,
                    }
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResultData(
                phenopacket=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_variant_prioritisation_no_threshold(self):
        self.assess_variant_prioritisation.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records
        )

        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3_126741108_G_A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16_133564345_C_T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_variant_prioritisation_fails_pvalue_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.01
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3_126741108_G_A",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16_133564345_C_T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_variant_prioritisation_meets_pvalue_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3_126741108_G_A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16_133564345_C_T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_variant_prioritisation_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3_126741108_G_A",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16_133564345_C_T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )

    def test_assess_variant_prioritisation_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_pvalue)
        assess_with_threshold.threshold = 0.1
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3_126741108_G_A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16_133564345_C_T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )


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
                        "Variant": "3_12563453454_C_T",
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
        self.assertTrue(result.equals(self.gene_rank_comparisons.generate_dataframe()))

    def test_generate_variant_dataframe(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Variant": "3_12563453454_C_T",
                    "/path/to/results_directory1": 9,
                    "/path/to/results_directory2": 3,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.variant_rank_comparisons.generate_dataframe()))

    def test_calculate_gene_rank_difference(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Gene": "GCDH",
                    "/path/to/results_directory1": 1,
                    "/path/to/results_directory2": 5,
                    "rank_decrease": 4,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.gene_rank_comparisons.calculate_rank_difference()))

    def test_calculate_variant_rank_difference(self):
        result = pd.DataFrame(
            [
                {
                    "Phenopacket": "phenopacket1.json",
                    "Variant": "3_12563453454_C_T",
                    "/path/to/results_directory1": 9,
                    "/path/to/results_directory2": 3,
                    "rank_decrease": -6,
                }
            ]
        )
        result.index += 1
        self.assertTrue(result.equals(self.variant_rank_comparisons.calculate_rank_difference()))
