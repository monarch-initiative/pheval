import copy
import unittest
from pathlib import Path

from phenopackets import (
    Diagnosis,
    Family,
    File,
    GeneDescriptor,
    GenomicInterpretation,
    Individual,
    Interpretation,
    MetaData,
    OntologyClass,
    Pedigree,
    Phenopacket,
    PhenotypicFeature,
    Resource,
    VariantInterpretation,
    VariationDescriptor,
    VcfRecord,
)

from pheval.prepare.custom_exceptions import IncorrectFileFormatError
from pheval.utils.phenopacket_utils import (
    IncompatibleGenomeAssemblyError,
    PhenopacketRebuilder,
    PhenopacketUtil,
    ProbandCausativeVariant,
)

proband = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=[
        PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
        PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
        PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
        PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
        PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
    ],
    interpretations=[
        Interpretation(
            id="test-subject-1-int",
            progress_status="SOLVED",
            diagnosis=Diagnosis(
                genomic_interpretations=[
                    GenomicInterpretation(
                        subject_or_biosample_id="test-subject-1",
                        interpretation_status=4,
                        variant_interpretation=VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=VariationDescriptor(
                                gene_context=GeneDescriptor(
                                    value_id="NCBIGene:2245", symbol="FGD1"
                                ),
                                vcf_record=VcfRecord(
                                    genome_assembly="GRCh37",
                                    chrom="X",
                                    pos=54492285,
                                    ref="C",
                                    alt="T",
                                ),
                                allelic_state=OntologyClass(
                                    id="GENO:0000134",
                                    label="hemizygous",
                                ),
                            ),
                        ),
                    ),
                    GenomicInterpretation(
                        subject_or_biosample_id="test-subject-1",
                        interpretation_status=4,
                        variant_interpretation=VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=VariationDescriptor(
                                gene_context=GeneDescriptor(value_id="HGNC:18654", symbol="RTTN"),
                                vcf_record=VcfRecord(
                                    genome_assembly="GRCh37",
                                    chrom="18",
                                    pos=67691994,
                                    ref="G",
                                    alt="A",
                                ),
                                allelic_state=OntologyClass(
                                    id="GENO:0000402", label="compound heterozygous"
                                ),
                            ),
                        ),
                    ),
                ]
            ),
        )
    ],
)

phenopacket = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=[
        PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
        PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
        PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
        PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
        PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
        PhenotypicFeature(
            type=OntologyClass(id="HP:0008494", label="Inferior lens subluxation"), excluded=True
        ),
    ],
    interpretations=[
        Interpretation(
            id="test-subject-1-int",
            progress_status="SOLVED",
            diagnosis=Diagnosis(
                genomic_interpretations=[
                    GenomicInterpretation(
                        subject_or_biosample_id="test-subject-1",
                        interpretation_status=4,
                        variant_interpretation=VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=VariationDescriptor(
                                gene_context=GeneDescriptor(
                                    value_id="NCBIGene:2245", symbol="FGD1"
                                ),
                                vcf_record=VcfRecord(
                                    genome_assembly="GRCh37",
                                    chrom="X",
                                    pos=54492285,
                                    ref="C",
                                    alt="T",
                                ),
                                allelic_state=OntologyClass(
                                    id="GENO:0000134",
                                    label="hemizygous",
                                ),
                            ),
                        ),
                    ),
                    GenomicInterpretation(
                        subject_or_biosample_id="test-subject-1",
                        interpretation_status=4,
                        variant_interpretation=VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=VariationDescriptor(
                                gene_context=GeneDescriptor(value_id="HGNC:18654", symbol="RTTN"),
                                vcf_record=VcfRecord(
                                    genome_assembly="GRCh37",
                                    chrom="18",
                                    pos=67691994,
                                    ref="G",
                                    alt="A",
                                ),
                                allelic_state=OntologyClass(
                                    id="GENO:0000402", label="compound heterozygous"
                                ),
                            ),
                        ),
                    ),
                ]
            ),
        )
    ],
    files=[
        File(
            uri="test/path/to/test_1.vcf",
            file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
        ),
        File(
            uri="test_1.ped",
            file_attributes={"fileFormat": "PED", "genomeAssembly": "GRCh37"},
        ),
    ],
    meta_data=MetaData(
        created_by="pheval-converter",
        resources=[
            Resource(
                id="hp",
                name="human phenotype ontology",
                url="http://purl.obolibrary.org/obo/hp.owl",
                version="hp/releases/2019-11-08",
                namespace_prefix="HP",
                iri_prefix="http://purl.obolibrary.org/obo/HP_",
            )
        ],
        phenopacket_schema_version="2.0",
    ),
)

family = Family(
    id="test-family-1",
    proband=proband,
    pedigree=Pedigree(
        persons=[
            Pedigree.Person(
                family_id="test-family-1",
                individual_id="test-subject-1",
                paternal_id="MOTHER",
                maternal_id="FATHER",
                sex=1,
                affected_status=1,
            )
        ]
    ),
    files=phenopacket.files,
    meta_data=phenopacket.meta_data,
)


class TestPhenopacketUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket = PhenopacketUtil(phenopacket)
        cls.family = PhenopacketUtil(family)

    def test_phenotypic_features(self):
        self.assertTrue(
            "HP:" in pheno_phenotypic_feature.type.id
            for pheno_phenotypic_feature in self.phenopacket.phenotypic_features()
        )
        self.assertTrue(
            "HP:" in family_phenotypic_feature.type.id
            for family_phenotypic_feature in self.family.phenotypic_features()
        )
        self.assertEqual(len(self.phenopacket.phenotypic_features()), 6)
        self.assertEqual(len(self.family.phenotypic_features()), 5)

    def test_remove_excluded_phenotypic_features(self):
        self.assertEqual(len(self.phenopacket.remove_excluded_phenotypic_features()), 5)
        self.assertEqual(len(self.family.remove_excluded_phenotypic_features()), 5)

    def test_interpretations(self):
        for interpretation in self.phenopacket.interpretations() and self.family.interpretations():
            self.assertTrue(hasattr(interpretation, "progress_status"))
            self.assertTrue(hasattr(interpretation, "diagnosis"))
            self.assertTrue(hasattr(interpretation.diagnosis, "genomic_interpretations"))
            for genomic_interpretation in interpretation.diagnosis.genomic_interpretations:
                self.assertTrue(hasattr(genomic_interpretation, "variant_interpretation"))
                self.assertTrue(
                    hasattr(genomic_interpretation.variant_interpretation, "variation_descriptor")
                )
                self.assertTrue(
                    hasattr(
                        genomic_interpretation.variant_interpretation.variation_descriptor,
                        "vcf_record",
                    )
                )

    def test_causative_variants(self):
        self.assertEqual(
            len(self.phenopacket.causative_variants(Path("test-phenopacket-1.json"))), 2
        )
        self.assertFalse(
            self.phenopacket.causative_variants(Path("test-phenopacket-1.json"))[0]
            == self.phenopacket.causative_variants(Path("test-phenopacket-1.json"))[1]
        )
        for causative_variant in self.phenopacket.causative_variants(
            Path("test-phenopacket-1.json")
        ):
            self.assertEqual(type(causative_variant), ProbandCausativeVariant)

    def test_files(self):
        self.assertEqual(len(self.phenopacket.files()), 2)
        self.assertEqual(len(self.phenopacket.files()), 2)
        for file in self.phenopacket.files():
            self.assertTrue(type(file) is File)

    def test_vcf_file_data(self):
        vcf_file_data = self.phenopacket.vcf_file_data(
            Path("test-phenopacket-1.json"), Path("input_dir")
        )
        self.assertTrue(vcf_file_data.uri == "input_dir/test_1.vcf")
        self.assertTrue(vcf_file_data.file_attributes["fileFormat"] == "VCF")
        incorrect_genome_assembly_ppacket = copy.copy(phenopacket)
        del incorrect_genome_assembly_ppacket.files[:]
        incorrect_genome_assembly = [
            File(
                uri="test/path/to/test_1.vcf",
                file_attributes={"fileFormat": "VCF", "genomeAssembly": "hg10"},
            )
            if file.file_attributes["fileFormat"] == "VCF"
            else file
            for file in phenopacket.files
        ]
        incorrect_genome_assembly_ppacket.files.extend(incorrect_genome_assembly)
        with self.assertRaises(IncompatibleGenomeAssemblyError):
            PhenopacketUtil(incorrect_genome_assembly_ppacket).vcf_file_data(
                Path("test-phenopacket-1.json"), Path("input_dir")
            )
        incorrect_file_format_ppacket = copy.copy(phenopacket)
        del incorrect_file_format_ppacket.files[:]
        incorrect_file_format = [
            File(
                uri="test/path/to/test_1.ped",
                file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
            )
            if file.file_attributes["fileFormat"] == "VCF"
            else file
            for file in phenopacket.files
        ]
        incorrect_file_format_ppacket.files.extend(incorrect_file_format)
        with self.assertRaises(IncorrectFileFormatError):
            PhenopacketUtil(incorrect_file_format_ppacket).vcf_file_data(
                Path("test-phenopacket-1.json"), Path("input_dir")
            )

    def test_diagnosed_genes(self):
        self.assertEqual(len(self.phenopacket.diagnosed_genes()), 2)
        self.assertTrue(
            all(gene in ["RTTN", "FGD1"] for gene in self.phenopacket.diagnosed_genes())
        )

    def test_diagnosed_variants(self):
        self.assertEqual(len(self.phenopacket.diagnosed_variants()), 2)
        for diagnosed_variant in self.phenopacket.diagnosed_variants():
            self.assertTrue(
                diagnosed_variant.gene == "RTTN"
                and diagnosed_variant.chrom == "18"
                and diagnosed_variant.pos == 67691994
                and diagnosed_variant.ref == "G"
                and diagnosed_variant.alt == "A"
                or diagnosed_variant.gene == "FGD1"
                and diagnosed_variant.chrom == "X"
                and diagnosed_variant.pos == 54492285
                and diagnosed_variant.ref == "C"
                and diagnosed_variant.alt == "T"
            )


class TestPhenopacketRebuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket_rebuilder = PhenopacketRebuilder(phenopacket)
        cls.family_rebuilder = PhenopacketRebuilder(family)

    def test_add_randomised_hpo(self):
        self.phenopacket_rebuilder.add_randomised_hpo(
            [
                PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="Hypertelorism")),
                PhenotypicFeature(type=OntologyClass(id="HP:RANDOM", label="RANDOM")),
            ]
        )
        self.family_rebuilder.add_randomised_hpo(
            [
                PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="Hypertelorism")),
                PhenotypicFeature(type=OntologyClass(id="HP:RANDOM", label="RANDOM")),
            ]
        )
        self.assertEqual(
            len(self.phenopacket_rebuilder.phenopacket_contents.phenotypic_features), 2
        )
        self.assertEqual(
            len(self.family_rebuilder.phenopacket_contents.proband.phenotypic_features), 2
        )
        for p_f in self.phenopacket_rebuilder.phenopacket_contents.phenotypic_features:
            self.assertTrue(
                p_f.type.id == "HP:0000316"
                and p_f.type.label == "Hypertelorism"
                or p_f.type.id == "HP:RANDOM"
                and p_f.type.label == "RANDOM"
            )

    def test_add_created_vcf_path(self):
        self.phenopacket_rebuilder.add_created_vcf_path(
            Path("input_dir/test_vcf_dir/test_1.vcf"), "GRCh37"
        )
        vcf_file = [
            file
            for file in self.phenopacket_rebuilder.phenopacket_contents.files
            if file.file_attributes["fileFormat"] == "VCF"
        ][0]
        self.assertEqual(vcf_file.uri, str(Path("input_dir/test_vcf_dir/test_1.vcf").absolute()))
