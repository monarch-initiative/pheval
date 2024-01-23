import unittest
from collections import defaultdict
from copy import copy
from pathlib import Path, PosixPath

from pheval.analyse.binary_classification_stats import BinaryClassificationStats
from pheval.analyse.disease_prioritisation_analysis import AssessDiseasePrioritisation
from pheval.analyse.gene_prioritisation_analysis import AssessGenePrioritisation
from pheval.analyse.prioritisation_rank_recorder import PrioritisationRankRecorder
from pheval.analyse.prioritisation_result_types import (
    DiseasePrioritisationResult,
    GenePrioritisationResult,
    VariantPrioritisationResult,
)
from pheval.analyse.rank_stats import RankStats
from pheval.analyse.variant_prioritisation_analysis import AssessVariantPrioritisation
from pheval.post_processing.post_processing import (
    RankedPhEvalDiseaseResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
)
from pheval.utils.phenopacket_utils import GenomicVariant, ProbandCausativeGene, ProbandDisease


class TestPrioritisationRankRecorder(unittest.TestCase):
    def setUp(self) -> None:
        self.add_new_phenopacket_variant_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            VariantPrioritisationResult(
                Path("/path/to/phenopacket-2.json"), GenomicVariant("1", 4896347, "C", "T"), 9
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_directory_variant_record = PrioritisationRankRecorder(
            0,
            Path("directory2"),
            VariantPrioritisationResult(
                Path("/path/to/phenopacket-1.json"), GenomicVariant("12", 120434, "A", "G"), 9
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_phenopacket_gene_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            GenePrioritisationResult(Path("/path/to/phenopacket-2.json"), "GENE", 7),
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
            GenePrioritisationResult(Path("/path/to/phenopacket-1.json"), "LARGE1", 1),
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
        self.add_new_directory_disease_record = PrioritisationRankRecorder(
            0,
            Path("directory2"),
            DiseasePrioritisationResult(
                Path("/path/to/phenopacket-1.json"),
                ProbandDisease(disease_name="DISEASE1", disease_identifier="OMIM:12345"),
                1,
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 4,
                    }
                },
            ),
        )
        self.add_new_phenopacket_disease_record = PrioritisationRankRecorder(
            1,
            Path("directory1"),
            DiseasePrioritisationResult(
                Path("/path/to/phenopacket-2.json"),
                ProbandDisease(disease_name="DISEASE1", disease_identifier="OMIM:12345"),
                7,
            ),
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 4,
                    }
                },
            ),
        )

    def test__record_gene_rank_new_directory(self):
        self.assertEqual(
            self.add_new_directory_gene_record.run_comparison,
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
        self.add_new_directory_gene_record.record_rank()
        self.assertEqual(
            self.add_new_directory_gene_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Gene": "LARGE1",
                        PosixPath("directory1"): 4,
                        PosixPath("directory2"): 1,
                    }
                },
            ),
        )

    def test__record_gene_rank_new_phenopacket(self):
        self.assertEqual(
            self.add_new_phenopacket_gene_record.run_comparison,
            defaultdict(
                dict,
                {0: {"Phenopacket": "phenopacket-1.json", "Gene": "LARGE1", Path("directory1"): 4}},
            ),
        )
        self.add_new_phenopacket_gene_record.record_rank()
        self.assertEqual(
            self.add_new_phenopacket_gene_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Gene": "LARGE1",
                        Path("directory1"): 4,
                    },
                    1: {"Phenopacket": "phenopacket-2.json", "Gene": "GENE", Path("directory1"): 7},
                },
            ),
        )

    def test__variant_rank_new_directory(self):
        self.assertEqual(
            self.add_new_directory_variant_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_directory_variant_record.record_rank()
        self.assertEqual(
            self.add_new_directory_variant_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        PosixPath("directory1"): 3,
                        PosixPath("directory2"): 9,
                    }
                },
            ),
        )

    def test__variant_rank_new_phenopacket(self):
        self.assertEqual(
            self.add_new_phenopacket_variant_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        Path("directory1"): 3,
                    }
                },
            ),
        )
        self.add_new_phenopacket_variant_record.record_rank()
        self.assertEqual(
            self.add_new_phenopacket_variant_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        PosixPath("directory1"): 3,
                    },
                    1: {
                        "Phenopacket": "phenopacket-2.json",
                        "Variant": "1-4896347-C-T",
                        PosixPath("directory1"): 9,
                    },
                },
            ),
        )

    def test__disease_rank_new_directory(self):
        self.assertEqual(
            self.add_new_directory_disease_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        Path("directory1"): 4,
                    }
                },
            ),
        )
        self.add_new_directory_disease_record.record_rank()
        self.assertEqual(
            self.add_new_directory_disease_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 4,
                        PosixPath("directory2"): 1,
                    }
                },
            ),
        )

    def test__disease_rank_new_phenopacket(self):
        self.assertEqual(
            self.add_new_phenopacket_disease_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        Path("directory1"): 4,
                    }
                },
            ),
        )
        self.add_new_phenopacket_disease_record.record_rank()
        self.assertEqual(
            self.add_new_phenopacket_disease_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 4,
                    },
                    1: {
                        "Phenopacket": "phenopacket-2.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 7,
                    },
                },
            ),
        )

    def test_record_rank_gene(self):
        self.add_new_directory_gene_record.record_rank()
        self.assertEqual(
            self.add_new_directory_gene_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Gene": "LARGE1",
                        PosixPath("directory1"): 4,
                        PosixPath("directory2"): 1,
                    }
                },
            ),
        )

    def test_record_rank_variant(self):
        self.add_new_directory_variant_record.record_rank()
        self.assertEqual(
            self.add_new_directory_variant_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Variant": "12-120434-A-G",
                        PosixPath("directory1"): 3,
                        PosixPath("directory2"): 9,
                    }
                },
            ),
        )

    def test_record_rank_disease(self):
        self.add_new_directory_disease_record.record_rank()
        self.assertEqual(
            self.add_new_directory_disease_record.run_comparison,
            defaultdict(
                dict,
                {
                    0: {
                        "Phenopacket": "phenopacket-1.json",
                        "Disease": "OMIM:12345",
                        PosixPath("directory1"): 4,
                        PosixPath("directory2"): 1,
                    }
                },
            ),
        )


class TestAssessGenePrioritisation(unittest.TestCase):
    def setUp(self) -> None:
        self.assess_gene_prioritisation = AssessGenePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_gene_results=[
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
            threshold=0.0,
            score_order="descending",
            proband_causative_genes=[
                ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                ProbandCausativeGene(gene_symbol="LARGE1", gene_identifier="ENSG00000133424"),
            ],
        )
        self.assess_gene_prioritisation_ascending_order = AssessGenePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_gene_results=[
                RankedPhEvalGeneResult(
                    gene_symbol="SPNS1",
                    gene_identifier="ENSG00000169682",
                    score=0.3765,
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
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=4,
                ),
            ],
            threshold=0.0,
            score_order="ascending",
            proband_causative_genes=[
                ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                ProbandCausativeGene(gene_symbol="LARGE1", gene_identifier="ENSG00000133424"),
            ],
        )
        self.gene_rank_stats = RankStats(0, 0, 0, 0, 0)
        self.gene_rank_records = defaultdict(dict)
        self.binary_classification_stats = BinaryClassificationStats()

    def test_record_gene_prioritisation_match(self):
        self.assertEqual(
            self.assess_gene_prioritisation._record_gene_prioritisation_match(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )

    def test_assess_gene_with_ascending_order_threshold_fails_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.1
        self.assertEqual(
            assess_ascending_order_threshold._assess_gene_with_threshold_ascending_order(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.gene_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_gene_with_ascending_order_threshold_meets_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.9
        self.assertEqual(
            assess_ascending_order_threshold._assess_gene_with_threshold_ascending_order(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_gene_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_gene_with_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.gene_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_gene_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.5
        self.assertEqual(
            assess_with_threshold._assess_gene_with_threshold(
                gene=ProbandCausativeGene(gene_symbol="PLXNA1", gene_identifier="ENSG00000114554"),
                result_entry=RankedPhEvalGeneResult(
                    gene_symbol="PLXNA1",
                    gene_identifier="ENSG00000114554",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.gene_rank_stats,
            ),
            GenePrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"), gene="PLXNA1", rank=1
            ),
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_gene_prioritisation_no_threshold(self):
        self.assess_gene_prioritisation.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[1.0]),
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
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=3, false_positives=0, false_negatives=0
            ),
        )

    def test_assess_gene_prioritisation_threshold_fails_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.01
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=2, reciprocal_ranks=[]),
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
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=3, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_gene_prioritisation_threshold_meets_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[0.25]),
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
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=2, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_gene_prioritisation_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=2, reciprocal_ranks=[]),
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
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=3, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_gene_prioritisation_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_gene_prioritisation)
        assess_with_threshold.threshold = 0.1
        assess_with_threshold.assess_gene_prioritisation(
            self.gene_rank_stats, self.gene_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.gene_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[1.0]),
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
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=3, false_positives=0, false_negatives=0
            ),
        )


class TestAssessVariantPrioritisation(unittest.TestCase):
    def setUp(self) -> None:
        variant_results = [
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
        ]
        self.assess_variant_prioritisation = AssessVariantPrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_variant_results=variant_results,
            threshold=0.0,
            score_order="descending",
            proband_causative_variants=[
                GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                GenomicVariant(chrom="16", pos=133564345, ref="C", alt="T"),
            ],
        )
        self.assess_variant_prioritisation_ascending_order = AssessVariantPrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_variant_results=variant_results,
            threshold=0.0,
            score_order="ascending",
            proband_causative_variants=[
                GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                GenomicVariant(chrom="16", pos=133564345, ref="C", alt="T"),
            ],
        )
        self.variant_rank_stats = RankStats()
        self.variant_rank_records = defaultdict(dict)
        self.binary_classification_stats = BinaryClassificationStats()

    def test_record_variant_prioritisation_match(self):
        self.assertEqual(
            self.assess_variant_prioritisation._record_variant_prioritisation_match(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )

    def test_assess_variant_with_ascending_order_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold._assess_variant_with_threshold_ascending_order(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                rank_stats=self.variant_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_variant_with_ascending_order_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_variant_with_threshold_ascending_order(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_variant_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_variant_with_threshold(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                rank_stats=self.variant_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_variant_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.01
        self.assertEqual(
            assess_with_threshold._assess_variant_with_threshold(
                result_entry=RankedPhEvalVariantResult(
                    chromosome="3",
                    start=126741108,
                    end=126741108,
                    ref="G",
                    alt="A",
                    score=0.0484,
                    rank=1,
                ),
                rank_stats=self.variant_rank_stats,
            ),
            VariantPrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"),
                variant=GenomicVariant(chrom="3", pos=126741108, ref="G", alt="A"),
                rank=1,
            ),
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_variant_prioritisation_no_threshold(self):
        self.assess_variant_prioritisation.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records, self.binary_classification_stats
        )

        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3-126741108-G-A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16-133564345-C-T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=0, false_positives=2, false_negatives=0
            ),
        )

    def test_assess_variant_prioritisation_fails_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.01
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3-126741108-G-A",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16-133564345-C-T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=0, false_positives=3, false_negatives=1
            ),
        )

    def test_assess_variant_prioritisation_meets_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3-126741108-G-A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16-133564345-C-T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=0, false_positives=2, false_negatives=0
            ),
        )

    def test_assess_variant_prioritisation_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=2, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3-126741108-G-A",
                    Path("/path/to/results_dir"): 0,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16-133564345-C-T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=0, false_positives=3, false_negatives=1
            ),
        )

    def test_assess_variant_prioritisation_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_variant_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.1
        assess_with_threshold.assess_variant_prioritisation(
            self.variant_rank_stats, self.variant_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.variant_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=2, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.variant_rank_records,
            {
                1: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "3-126741108-G-A",
                    Path("/path/to/results_dir"): 1,
                },
                2: {
                    "Phenopacket": "phenopacket.json",
                    "Variant": "16-133564345-C-T",
                    Path("/path/to/results_dir"): 0,
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=0, false_positives=2, false_negatives=0
            ),
        )


class TestAssessDiseasePrioritisation(unittest.TestCase):
    def setUp(self) -> None:
        self.assess_disease_prioritisation = AssessDiseasePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_disease_results=[
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=1.0,
                    rank=1,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 2",
                    disease_identifier="OMIM:231680",
                    score=0.5,
                    rank=2,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 3",
                    disease_identifier="OMIM:231690",
                    score=0.5,
                    rank=2,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 4",
                    disease_identifier="OMIM:231700",
                    score=0.3,
                    rank=4,
                ),
            ],
            threshold=0.0,
            score_order="descending",
            proband_diseases=[
                ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                )
            ],
        )
        self.assess_disease_prioritisation_ascending_order = AssessDiseasePrioritisation(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            results_dir=Path("/path/to/results_dir"),
            standardised_disease_results=[
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 4",
                    disease_identifier="OMIM:231690",
                    score=0.3765,
                    rank=1,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 2",
                    disease_identifier="OMIM:231680",
                    score=0.5777,
                    rank=2,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 3",
                    disease_identifier="OMIM:231690",
                    score=0.5777,
                    rank=2,
                ),
                RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=4,
                ),
            ],
            threshold=0.0,
            score_order="ascending",
            proband_diseases=[
                ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                )
            ],
        )
        self.disease_rank_stats = RankStats(0, 0, 0, 0, 0)
        self.disease_rank_records = defaultdict(dict)
        self.binary_classification_stats = BinaryClassificationStats()

    def test_record_disease_prioritisation_match(self):
        self.assertEqual(
            self.assess_disease_prioritisation._record_disease_prioritisation_match(
                disease=ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                ),
                result_entry=RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.disease_rank_stats,
            ),
            DiseasePrioritisationResult(
                phenopacket_path=PosixPath("/path/to/phenopacket.json"),
                disease=ProbandDisease(
                    disease_name="Glutaric aciduria type 1", disease_identifier="OMIM:231670"
                ),
                rank=1,
            ),
        )

    def test_assess_disease_with_ascending_order_threshold_fails_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.1
        self.assertEqual(
            assess_ascending_order_threshold._assess_disease_with_threshold_ascending_order(
                disease=ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                ),
                result_entry=RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.disease_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_disease_with_ascending_order_threshold_meets_cutoff(self):
        assess_ascending_order_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_ascending_order_threshold.threshold = 0.9
        self.assertEqual(
            assess_ascending_order_threshold._assess_disease_with_threshold_ascending_order(
                disease=ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                ),
                result_entry=RankedPhEvalDiseaseResult(
                    disease_name="Glutaric aciduria type 1",
                    disease_identifier="OMIM:231670",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.disease_rank_stats,
            ),
            DiseasePrioritisationResult(
                phenopacket_path=PosixPath("/path/to/phenopacket.json"),
                disease=ProbandDisease(
                    disease_name="Glutaric aciduria type 1", disease_identifier="OMIM:231670"
                ),
                rank=1,
            ),
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_disease_with_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 0.9
        self.assertEqual(
            assess_with_threshold._assess_disease_with_threshold(
                disease=ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                ),
                result_entry=RankedPhEvalDiseaseResult(
                    disease_identifier="OMIM:231670",
                    disease_name="Glutaric aciduria type 1",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.disease_rank_stats,
            ),
            None,
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=0, reciprocal_ranks=[]),
        )

    def test_assess_disease_with_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 0.5
        self.assertEqual(
            assess_with_threshold._assess_disease_with_threshold(
                disease=ProbandDisease(
                    disease_identifier="OMIM:231670", disease_name="Glutaric aciduria type 1"
                ),
                result_entry=RankedPhEvalDiseaseResult(
                    disease_identifier="OMIM:231670",
                    disease_name="Glutaric aciduria type 1",
                    score=0.8764,
                    rank=1,
                ),
                rank_stats=self.disease_rank_stats,
            ),
            DiseasePrioritisationResult(
                phenopacket_path=Path("/path/to/phenopacket.json"),
                disease=ProbandDisease(
                    disease_name="Glutaric aciduria type 1", disease_identifier="OMIM:231670"
                ),
                rank=1,
            ),
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=0, reciprocal_ranks=[1.0]),
        )

    def test_assess_disease_prioritisation_no_threshold(self):
        self.assess_disease_prioritisation.assess_disease_prioritisation(
            self.disease_rank_stats, self.disease_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=1, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.disease_rank_records,
            {
                1: {
                    Path("/path/to/results_dir"): 1,
                    "Disease": "OMIM:231670",
                    "Phenopacket": "phenopacket.json",
                }
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=3, false_positives=0, false_negatives=0
            ),
        )

    def test_assess_disease_prioritisation_threshold_fails_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.01
        assess_with_threshold.assess_disease_prioritisation(
            self.disease_rank_stats, self.disease_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=1, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.disease_rank_records,
            {
                1: {
                    Path("/path/to/results_dir"): 0,
                    "Disease": "OMIM:231670",
                    "Phenopacket": "phenopacket.json",
                }
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=3, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_disease_prioritisation_threshold_meets_ascending_order_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation_ascending_order)
        assess_with_threshold.threshold = 0.9
        assess_with_threshold.assess_disease_prioritisation(
            self.disease_rank_stats, self.disease_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=0, top3=0, top5=1, top10=1, found=1, total=1, reciprocal_ranks=[0.25]),
        )
        self.assertEqual(
            self.disease_rank_records,
            {
                1: {
                    Path("/path/to/results_dir"): 4,
                    "Disease": "OMIM:231670",
                    "Phenopacket": "phenopacket.json",
                },
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=2, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_disease_prioritisation_threshold_fails_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 1.0
        assess_with_threshold.assess_disease_prioritisation(
            self.disease_rank_stats, self.disease_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=0, top3=0, top5=0, top10=0, found=0, total=1, reciprocal_ranks=[]),
        )
        self.assertEqual(
            self.disease_rank_records,
            {
                1: {
                    Path("/path/to/results_dir"): 0,
                    "Disease": "OMIM:231670",
                    "Phenopacket": "phenopacket.json",
                }
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=0, true_negatives=3, false_positives=1, false_negatives=1
            ),
        )

    def test_assess_disease_prioritisation_threshold_meets_cutoff(self):
        assess_with_threshold = copy(self.assess_disease_prioritisation)
        assess_with_threshold.threshold = 0.1
        assess_with_threshold.assess_disease_prioritisation(
            self.disease_rank_stats, self.disease_rank_records, self.binary_classification_stats
        )
        self.assertEqual(
            self.disease_rank_stats,
            RankStats(top=1, top3=1, top5=1, top10=1, found=1, total=1, reciprocal_ranks=[1.0]),
        )
        self.assertEqual(
            self.disease_rank_records,
            {
                1: {
                    Path("/path/to/results_dir"): 1,
                    "Disease": "OMIM:231670",
                    "Phenopacket": "phenopacket.json",
                }
            },
        )
        self.assertEqual(
            self.binary_classification_stats,
            BinaryClassificationStats(
                true_positives=1, true_negatives=3, false_positives=0, false_negatives=0
            ),
        )
