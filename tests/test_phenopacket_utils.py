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
from pheval.utils.phenopacket_utils import (  # GeneIdentifierUpdater,; create_hgnc_dict,
    GeneIdentifierUpdater,
    IncompatibleGenomeAssemblyError,
    PhenopacketRebuilder,
    PhenopacketUtil,
    ProbandCausativeVariant,
    VariantData,
    create_hgnc_dict,
)

interpretations = [
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
                            gene_context=GeneDescriptor(value_id="NCBIGene:2245", symbol="FGD1"),
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
]
updated_interpretations = [
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
                                value_id="ENSG00000102302",
                                symbol="FGD1",
                                alternate_ids=[
                                    "HGNC:3663",
                                    "ncbigene:2245",
                                    "ensembl:ENSG00000102302",
                                    "symbol:FGD1",
                                ],
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
                            gene_context=GeneDescriptor(
                                value_id="ENSG00000176225",
                                symbol="RTTN",
                                alternate_ids=[
                                    "HGNC:18654",
                                    "ncbigene:25914",
                                    "ensembl:ENSG00000176225",
                                    "symbol:RTTN",
                                ],
                            ),
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
]

phenotypic_features_none_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
]
phenotypic_features_with_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0008494", label="Inferior lens subluxation"), excluded=True
    ),
]

phenotypic_features_all_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly"), excluded=True),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy"), excluded=True),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage"), excluded=True
    ),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0003150", label="Glutaric aciduria"), excluded=True
    ),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia"), excluded=True),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0008494", label="Inferior lens subluxation"), excluded=True
    ),
]

proband = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_none_excluded,
    interpretations=interpretations,
)

phenopacket_files = [
    File(
        uri="test/path/to/test_1.vcf",
        file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
    ),
    File(
        uri="test_1.ped",
        file_attributes={"fileFormat": "PED", "genomeAssembly": "GRCh37"},
    ),
]
incorrect_genome_assembly = [
    File(
        uri="test/path/to/test_1.vcf",
        file_attributes={"fileFormat": "VCF", "genomeAssembly": "hg10"},
    ),
    File(
        uri="test_1.ped",
        file_attributes={"fileFormat": "PED", "genomeAssembly": "hg10"},
    ),
]
incorrect_file_format = [
    File(
        uri="test/path/to/test_1.ped",
        file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
    ),
    File(
        uri="test_1.vcf",
        file_attributes={"fileFormat": "PED", "genomeAssembly": "GRCh37"},
    ),
]
phenopacket_metadata = MetaData(
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
)

phenopacket = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_with_excluded,
    interpretations=interpretations,
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
)
phenopacket_with_all_excluded_terms = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_all_excluded,
    interpretations=interpretations,
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
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
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
)
family_incorrect_files = Family(
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
    files=incorrect_genome_assembly,
    meta_data=phenopacket_metadata,
)

family_incorrect_file_format = Family(
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
    files=incorrect_file_format,
    meta_data=phenopacket_metadata,
)


class TestPhenopacketUtil(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket = PhenopacketUtil(phenopacket)
        cls.phenopacket_excluded_pf = PhenopacketUtil(phenopacket_with_all_excluded_terms)
        cls.family = PhenopacketUtil(family)
        cls.family_incorrect_files = PhenopacketUtil(family_incorrect_files)
        cls.family_incorrect_file_format = PhenopacketUtil(family_incorrect_file_format)

    def test_phenotypic_features_phenopacket(self):
        self.assertEqual(
            list(self.phenopacket.phenotypic_features()), phenotypic_features_with_excluded
        )
        self.assertEqual(
            list(self.phenopacket_excluded_pf.phenotypic_features()),
            phenotypic_features_all_excluded,
        )

    def test_phenotypic_features_family(self):
        self.assertEqual(list(self.family.phenotypic_features()), phenotypic_features_none_excluded)

    def test_observed_phenotypic_features(self):
        self.assertEqual(list(self.phenopacket_excluded_pf.observed_phenotypic_features()), [])
        self.assertEqual(
            list(self.phenopacket.observed_phenotypic_features()), phenotypic_features_none_excluded
        )

    def test_interpretations_phenopacket(self):
        self.assertEqual(list(self.phenopacket.interpretations()), interpretations)

    def test_interpretations_family(self):
        self.assertEqual(list(self.family.interpretations()), interpretations)

    def test_causative_variants(self):
        for causative_variant in self.phenopacket.causative_variants():
            self.assertEqual(type(causative_variant), ProbandCausativeVariant)

    def test_files_phenopacket(self):
        self.assertEqual(list(self.phenopacket.files()), phenopacket_files)

    def test_files_family(self):
        self.assertEqual(list(self.family.files()), phenopacket_files)

    def test_vcf_file_data(self):
        vcf_file_data = self.phenopacket.vcf_file_data(
            Path("test-phenopacket-1.json"), Path("input_dir")
        )
        self.assertEqual(
            vcf_file_data,
            File(
                uri="input_dir/test_1.vcf",
                file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
            ),
        )
        with self.assertRaises(IncompatibleGenomeAssemblyError):
            self.family_incorrect_files.vcf_file_data(
                Path("test-phenopacket-1.json"), Path("input_dir")
            )
        with self.assertRaises(IncorrectFileFormatError):
            self.family_incorrect_file_format.vcf_file_data(
                Path("test-phenopacket-1.json"), Path("input_dir")
            )

    def test_diagnosed_genes(self):
        self.assertEqual(set((["RTTN", "FGD1"])), set(self.phenopacket.diagnosed_genes()))

    def test_diagnosed_variants(self):
        self.assertEqual(
            list(self.phenopacket.diagnosed_variants()),
            [
                VariantData(chrom="X", pos=54492285, ref="C", alt="T", gene="FGD1"),
                VariantData(chrom="18", pos=67691994, ref="G", alt="A", gene="RTTN"),
            ],
        )


class TestPhenopacketRebuilder(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket_rebuilder = PhenopacketRebuilder(phenopacket)
        cls.family_rebuilder = PhenopacketRebuilder(family)
        cls.randomised_phenotype = [
            PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="Hypertelorism")),
            PhenotypicFeature(type=OntologyClass(id="HP:RANDOM", label="RANDOM")),
        ]

    def test_update_interpretations_phenopacket(self):
        phenopacket_updated_interpretations = self.phenopacket_rebuilder.update_interpretations(
            updated_interpretations
        )
        self.assertNotEqual(
            list(updated_interpretations),
            list(self.phenopacket_rebuilder.phenopacket.interpretations),
        )
        self.assertEqual(
            list(phenopacket_updated_interpretations.interpretations), updated_interpretations
        )

    def test_update_interpretations_family(self):
        family_updated_interpretations = self.family_rebuilder.update_interpretations(
            updated_interpretations
        )
        self.assertNotEqual(
            list(updated_interpretations),
            list(self.family_rebuilder.phenopacket.proband.interpretations),
        )
        self.assertEqual(
            list(family_updated_interpretations.proband.interpretations), updated_interpretations
        )

    def test_add_randomised_hpo_phenopacket(self):
        random_phenopacket = self.phenopacket_rebuilder.add_randomised_hpo(
            self.randomised_phenotype
        )
        self.assertNotEqual(
            random_phenopacket.phenotypic_features,
            self.phenopacket_rebuilder.phenopacket.phenotypic_features,
        )
        self.assertEqual(list(random_phenopacket.phenotypic_features), self.randomised_phenotype)

    def test_add_randomised_hpo_family(self):
        random_family = self.family_rebuilder.add_randomised_hpo(self.randomised_phenotype)
        self.assertNotEqual(
            random_family.proband.phenotypic_features,
            self.family_rebuilder.phenopacket.proband.phenotypic_features,
        )
        self.assertEqual(list(random_family.proband.phenotypic_features), self.randomised_phenotype)

    def test_add_created_vcf_path(self):
        updated_phenopacket = self.phenopacket_rebuilder.add_created_vcf_path(
            Path("input_dir/test_vcf_dir/test_1.vcf"), "GRCh37"
        )
        vcf_file = [
            file
            for file in updated_phenopacket.files
            if file.file_attributes["fileFormat"] == "VCF"
        ][0]
        self.assertEqual(vcf_file.uri, str(Path("input_dir/test_vcf_dir/test_1.vcf").absolute()))


class TestGeneIdentifierUpdater(unittest.TestCase):
    @classmethod
    def setUpClass(cls, hgnc_dict=None) -> None:
        if hgnc_dict is None:
            hgnc_dict = create_hgnc_dict()
        cls.gene_identifier_updater_ens = GeneIdentifierUpdater(hgnc_dict, "ensembl_id")
        cls.gene_identifier_updater_entrez = GeneIdentifierUpdater(hgnc_dict, "entrez_id")

    def test_find_identifier(self):
        self.assertEqual(self.gene_identifier_updater_ens.find_identifier("A2M"), "ENSG00000175899")
        self.assertEqual(self.gene_identifier_updater_ens.find_identifier("GBA"), "ENSG00000177628")
        self.assertEqual(self.gene_identifier_updater_entrez.find_identifier("A2M"), "2")
        self.assertEqual(self.gene_identifier_updater_entrez.find_identifier("GBA"), "2629")

    def test_find_alternate_ids(self):
        self.assertEqual(
            ["HGNC:4177", "ncbigene:2629", "ensembl:ENSG00000177628", "symbol:GBA1"],
            self.gene_identifier_updater_ens.find_alternate_ids("GBA"),
        )
        self.assertEqual(
            ["HGNC:7", "ncbigene:2", "ensembl:ENSG00000175899", "symbol:A2M"],
            self.gene_identifier_updater_entrez.find_alternate_ids("A2M"),
        )

    def test_update_genomic_interpretations_gene_identifier(self):
        self.assertEqual(
            list(
                self.gene_identifier_updater_ens.update_genomic_interpretations_gene_identifier(
                    phenopacket.interpretations
                )
            ),
            updated_interpretations,
        )
