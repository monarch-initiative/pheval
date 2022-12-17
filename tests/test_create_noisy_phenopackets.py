import os
import unittest

from phenopackets import OntologyClass, PhenotypicFeature

from pheval.prepare.create_noisy_phenopackets import HpoRandomiser, load_ontology

test_phenopacket = os.path.abspath("tests/input_dir/test_phenopacket_1.json")

phenotypic_features_single_term = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000965", label="Cutis marmorata"))
]
phenotypic_features = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
]


class TestHpoRandomiser(unittest.TestCase):
    hpo_ontology = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.hpo_ontology = load_ontology()
        cls.hpo_randomiser = HpoRandomiser(cls.hpo_ontology)

    def test_retrieve_hpo_term(self):
        self.assertEqual(
            self.hpo_randomiser.retrieve_hpo_term("HP:0000256"),
            PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
        )

    def test_retain_real_patient_terms(self):
        self.assertEqual(
            len(
                self.hpo_randomiser.retain_real_patient_terms(
                    phenotypic_features, number_of_real_id=2
                )
            ),
            2,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser.retain_real_patient_terms(
                    phenotypic_features, number_of_real_id=7
                )
            ),
            3,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser.retain_real_patient_terms(
                    phenotypic_features_single_term, number_of_real_id=7
                )
            ),
            1,
        )

    def test_convert_patient_terms_to_parent(self):
        self.assertEqual(
            len(
                self.hpo_randomiser.convert_patient_terms_to_parent(
                    phenotypic_features,
                    self.hpo_randomiser.retain_real_patient_terms(
                        phenotypic_features, number_of_real_id=2
                    ),
                    number_of_parent_terms=2,
                )
            ),
            2,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser.convert_patient_terms_to_parent(
                    phenotypic_features,
                    self.hpo_randomiser.retain_real_patient_terms(
                        phenotypic_features, number_of_real_id=7
                    ),
                    number_of_parent_terms=6,
                )
            ),
            2,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser.convert_patient_terms_to_parent(
                    phenotypic_features_single_term,
                    self.hpo_randomiser.retain_real_patient_terms(
                        phenotypic_features_single_term, number_of_real_id=7
                    ),
                    number_of_parent_terms=6,
                )
            ),
            0,
        )

    def test_create_random_hpo_terms(self):
        self.assertEqual(
            len(self.hpo_randomiser.create_random_hpo_terms(number_of_random_terms=3)), 3
        )
        self.assertEqual(
            len(self.hpo_randomiser.create_random_hpo_terms(number_of_random_terms=0)), 0
        )

    def test_randomise_hpo_terms(self):
        self.assertEqual(
            len(
                self.hpo_randomiser.randomise_hpo_terms(
                    phenotypic_features,
                    max_real_id=3,
                    number_of_parent_terms=4,
                    number_of_random_terms=3,
                )
            ),
            8,
        )
        self.assertEqual(
            len(
                self.hpo_randomiser.randomise_hpo_terms(
                    phenotypic_features_single_term,
                    max_real_id=3,
                    number_of_parent_terms=4,
                    number_of_random_terms=3,
                )
            ),
            4,
        )
