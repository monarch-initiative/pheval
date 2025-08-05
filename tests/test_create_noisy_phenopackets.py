import os
import unittest
from copy import copy

from phenopackets import OntologyClass, PhenotypicFeature

from pheval.prepare.create_noisy_phenopackets import HpoRandomiser, load_ontology

test_phenopacket = os.path.abspath("tests/input_dir/test_phenopacket_1.json")

phenotypic_features_single_term = [PhenotypicFeature(type=OntologyClass(id="HP:0000965", label="Cutis marmorata"))]
phenotypic_features = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
]
many_phenotypic_features = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
    PhenotypicFeature(type=OntologyClass(id="HP:0000965", label="Cutis marmorata")),
    PhenotypicFeature(type=OntologyClass(id="HP:0007843", label="Attenuation of retinal blood vessels")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001513", label="Obesity")),
    PhenotypicFeature(type=OntologyClass(id="HP:0000608", label="Macular degeneration")),
    PhenotypicFeature(type=OntologyClass(id="HP:0000486", label="Strabismus")),
]


class TestHpoRandomiser(unittest.TestCase):
    hpo_randomiser_0_75 = None
    hpo_randomiser_0_25 = None
    hpo_randomiser_default = None
    hpo_ontology = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.hpo_ontology = load_ontology()
        cls.hpo_randomiser_default = HpoRandomiser(cls.hpo_ontology, 0.5)
        cls.hpo_randomiser_0_25 = copy(cls.hpo_randomiser_default)
        cls.hpo_randomiser_0_25.scramble_factor = 0.25
        cls.hpo_randomiser_0_75 = copy(cls.hpo_randomiser_default)
        cls.hpo_randomiser_0_75.scramble_factor = 0.75

    def test_scramble_factor_proportions_default(self):
        self.assertEqual(
            self.hpo_randomiser_default.scramble_factor_proportions(phenotypic_features_single_term),
            1,
        )
        self.assertEqual(self.hpo_randomiser_default.scramble_factor_proportions(phenotypic_features), 2)
        self.assertEqual(self.hpo_randomiser_default.scramble_factor_proportions(many_phenotypic_features), 5)

    def test_scramble_factor_0_25(self):
        self.assertEqual(self.hpo_randomiser_0_25.scramble_factor_proportions(phenotypic_features_single_term), 1)
        self.assertEqual(self.hpo_randomiser_0_25.scramble_factor_proportions(phenotypic_features), 1)
        self.assertEqual(self.hpo_randomiser_0_25.scramble_factor_proportions(many_phenotypic_features), 2)

    def test_scramble_factor_0_75(self):
        self.assertEqual(self.hpo_randomiser_0_75.scramble_factor_proportions(phenotypic_features_single_term), 1)
        self.assertEqual(self.hpo_randomiser_0_75.scramble_factor_proportions(phenotypic_features), 4)
        self.assertEqual(self.hpo_randomiser_0_75.scramble_factor_proportions(many_phenotypic_features), 8)

    def test_retrieve_hpo_term(self):
        self.assertEqual(
            self.hpo_randomiser_default.retrieve_hpo_term("HP:0000256"),
            PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
        )

    def test_retain_real_patient_terms(self):
        self.assertEqual(
            len(self.hpo_randomiser_default.retain_real_patient_terms(phenotypic_features, 2)),
            3,
        )
        self.assertEqual(
            len(self.hpo_randomiser_default.retain_real_patient_terms(phenotypic_features_single_term, 2)),
            1,
        )

    def test_convert_patient_terms_to_parent(self):
        self.assertEqual(
            len(
                self.hpo_randomiser_default.convert_patient_terms_to_parent(
                    phenotypic_features,
                    self.hpo_randomiser_default.retain_real_patient_terms(phenotypic_features, 2),
                    2,
                )
            ),
            2,
        )

        self.assertEqual(
            len(
                self.hpo_randomiser_default.convert_patient_terms_to_parent(
                    phenotypic_features_single_term,
                    self.hpo_randomiser_default.retain_real_patient_terms(phenotypic_features_single_term, 1),
                    1,
                )
            ),
            0,
        )

    def test_create_random_hpo_terms(self):
        self.assertEqual(len(self.hpo_randomiser_default.create_random_hpo_terms(2)), 2)
        self.assertEqual(len(self.hpo_randomiser_default.create_random_hpo_terms(4)), 4)

    def test_randomise_hpo_terms(self):
        self.assertEqual(
            len(
                self.hpo_randomiser_default.randomise_hpo_terms(
                    phenotypic_features,
                )
            ),
            7,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser_default.randomise_hpo_terms(
                    many_phenotypic_features,
                )
            ),
            15,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser_default.randomise_hpo_terms(
                    phenotypic_features_single_term,
                )
            ),
            2,
        )
