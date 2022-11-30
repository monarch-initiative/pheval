import os.path
import pathlib
import unittest
from pathlib import Path

from phenopackets import OntologyClass, PhenotypicFeature

from pheval.utils import phenopacket_utils

current_path = str(pathlib.Path(__file__).parent.resolve())
test_phenopacket = Path(current_path + "/input_dir/test_phenopacket_1.json")
test_incorrect_assembly_phenopacket = Path(
    current_path + "/input_dir/test_phenopacket_incorrect_assembly_1.json"
)
test_incorrect_filetype_phenopacket = Path(
    current_path + "/input_dir/test_phenopacket_incorrect_filetype_1.json"
)


class TestPhenopacketReader(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket = phenopacket_utils.PhenopacketReader(Path(test_phenopacket))
        cls.incorrect_assembly_phenopacket = phenopacket_utils.PhenopacketReader(
            Path(test_incorrect_assembly_phenopacket)
        )
        cls.incorrect_filetype_phenopacket = phenopacket_utils.PhenopacketReader(
            Path(test_incorrect_filetype_phenopacket)
        )

    def test_phenotypic_features(self):
        self.assertTrue(
            "HP:" in phenotypic_feature.type.id
            for phenotypic_feature in self.phenopacket.phenotypic_features()
        )
        self.assertEqual(len(self.phenopacket.phenotypic_features()), 10)

    def test_remove_excluded_phenotypic_features(self):
        self.assertEqual(len(self.phenopacket.remove_excluded_phenotypic_features()), 9)

    def test_interpretations(self):
        for interpretation in self.phenopacket.interpretations():
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
        self.assertEqual(len(self.phenopacket.causative_variants()), 2)
        self.assertFalse(
            self.phenopacket.causative_variants()[0] == self.phenopacket.causative_variants()[1]
        )
        for causative_variant in self.phenopacket.causative_variants():
            self.assertEqual(type(causative_variant), phenopacket_utils.CausativeVariant)

    def test_files(self):
        self.assertEqual(len(self.phenopacket.files()), 1)
        self.assertTrue("test_phenopacket_1.vcf" in file.uri for file in self.phenopacket.files())
        self.assertEqual(
            [{"genomeAssembly": "GRCh37", "fileFormat": "VCF"}],
            [file.file_attributes for file in self.phenopacket.files()],
        )

    def test_vcf_file_data(self):
        self.assertEqual(
            self.phenopacket.vcf_file_data(Path("input_dir")),
            ("input_dir/test_phenopacket_1.vcf", "GRCh37"),
        )
        with self.assertRaises(phenopacket_utils.IncompatibleGenomeAssemblyError):
            self.incorrect_assembly_phenopacket.vcf_file_data(Path("input"))
        with self.assertRaises(phenopacket_utils.IncorrectFileFormatError):
            self.incorrect_filetype_phenopacket.vcf_file_data(Path("input"))

    def test_diagnosed_genes(self):
        self.assertEqual(len(self.phenopacket.diagnosed_genes()), 2)
        self.assertTrue(
            all(gene in ["RTTN", "FGD1"] for gene in self.phenopacket.diagnosed_genes())
        )

    def test_diagnosed_variants(self):
        self.assertEqual(len(self.phenopacket.diagnosed_variants()), 2)
        for diagnosed_variant in self.phenopacket.diagnosed_variants():
            self.assertTrue(
                diagnosed_variant["geneSymbol"] == "RTTN"
                or diagnosed_variant["geneSymbol"] == "FGD1"
            )
            self.assertTrue(
                hasattr(diagnosed_variant["variant"], "genome_assembly")
                and diagnosed_variant["variant"].genome_assembly == "GRCh37"
            )
            self.assertTrue(
                hasattr(diagnosed_variant["variant"], "chrom")
                and diagnosed_variant["variant"].chrom == "X"
                or diagnosed_variant["variant"].chrom == "18"
            )
            self.assertTrue(
                hasattr(diagnosed_variant["variant"], "pos")
                and diagnosed_variant["variant"].pos == 54492285
                or diagnosed_variant["variant"].pos == 67691994
            )
            self.assertTrue(
                hasattr(diagnosed_variant["variant"], "ref")
                and diagnosed_variant["variant"].ref == "C"
                or diagnosed_variant["variant"].ref == "G"
            )
            self.assertTrue(
                hasattr(diagnosed_variant["variant"], "alt")
                and diagnosed_variant["variant"].alt == "T"
                or diagnosed_variant["variant"].alt == "A"
            )


class TestPhenopacketRebuilder(unittest.TestCase):
    phenopacket_contents = phenopacket_utils.PhenopacketReader(Path(test_phenopacket))

    @classmethod
    def setUpClass(cls) -> None:
        cls.phenopacket_rebuilder = phenopacket_utils.PhenopacketRebuilder(cls.phenopacket_contents)
        cls.hpo_list = [
            PhenotypicFeature(type=OntologyClass(id="HP:0000316", label="Hypertelorism")),
            PhenotypicFeature(type=OntologyClass(id="HP:RANDOM", label="RANDOM")),
        ]

    def test_add_randomised_hpo(self):
        randomised_hpo_phenopacket = self.phenopacket_rebuilder.add_randomised_hpo(self.hpo_list)
        self.assertEqual(
            self.phenopacket_contents.pheno.subject, randomised_hpo_phenopacket.subject
        )
        self.assertNotEqual(
            self.phenopacket_contents.pheno.phenotypic_features,
            randomised_hpo_phenopacket.phenotypic_features,
        )
        for p_f in randomised_hpo_phenopacket.phenotypic_features:
            self.assertTrue(p_f.type.id == "HP:0000316" or p_f.type.id == "HP:RANDOM")
            self.assertTrue(p_f.type.label == "Hypertelorism" or p_f.type.label == "RANDOM")

    def test_create_json_message(self):
        phenopacket_to_create = self.phenopacket_rebuilder.add_randomised_hpo(self.hpo_list)
        message = self.phenopacket_rebuilder.create_json_message(phenopacket_to_create)
        self.assertEqual(type(message), str)


class TestPhenopacketWriter(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        ppacket = phenopacket_utils.PhenopacketReader(Path(test_phenopacket))
        ppacket_message = phenopacket_utils.PhenopacketRebuilder(ppacket).create_json_message(
            ppacket.pheno
        )
        cls.phenopacket_writer = phenopacket_utils.PhenopacketWriter(
            ppacket_message, Path("TEST-OUTPUT.json")
        )

    @classmethod
    def tearDownClass(cls) -> None:
        Path.unlink(Path("TEST-OUTPUT.json"))

    def test_write_file(self):
        self.phenopacket_writer.write_file()
        self.assertTrue(os.path.exists("TEST-OUTPUT.json"))
