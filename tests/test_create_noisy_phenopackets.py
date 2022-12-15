import os
import unittest

from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import OntologyClass, PhenotypicFeature

from pheval.prepare import create_noisy_phenopackets
from pheval.prepare.create_noisy_phenopackets import phenotypic_abnormality_entities

test_phenopacket = os.path.abspath("tests/input_dir/test_phenopacket_1.json")


class TestLoadOntology(unittest.TestCase):
    def test_load_ontology(self):
        self.assertFalse(create_noisy_phenopackets.ontology_loaded)
        create_noisy_phenopackets.load_ontology()
        self.assertTrue(create_noisy_phenopackets.ontology_loaded)


class TestRandomisePhenopackets(unittest.TestCase):
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

    @classmethod
    def setUpClass(cls) -> None:
        resource = OntologyResource(slug="hp.obo", local=False)
        ontology = ProntoImplementation(resource)
        cls.randomise_phenopacket = create_noisy_phenopackets.RandomisePhenopackets(
            ontology, cls.phenotypic_features, 3, 2, 3
        )
        cls.randomise_phenopacket_single_term = create_noisy_phenopackets.RandomisePhenopackets(
            ontology, cls.phenotypic_features_single_term, 3, 2, 3
        )
        cls.phenotypic_abnormality_children = phenotypic_abnormality_entities(ontology)

    def test_max_real_patient_id(self):
        self.assertTrue(len(self.randomise_phenopacket.retain_real_patient_terms()), 3)
        self.assertTrue(
            all(
                hpo_term in self.phenotypic_features
                for hpo_term in self.randomise_phenopacket.retain_real_patient_terms()
            )
        )
        self.assertTrue(len(self.randomise_phenopacket_single_term.retain_real_patient_terms()), 1)
        self.assertEqual(
            self.phenotypic_features_single_term,
            self.randomise_phenopacket_single_term.retain_real_patient_terms(),
        )

    def test_change_to_parent_term(self):
        self.assertTrue(len(self.randomise_phenopacket.convert_patient_terms_to_parent()), 2)
        self.assertFalse(
            all(
                hpo_term in self.phenotypic_features
                for hpo_term in self.randomise_phenopacket.convert_patient_terms_to_parent()
            )
        )
        self.assertEqual(
            self.randomise_phenopacket_single_term.convert_patient_terms_to_parent(), []
        )

    def test_random_hpo_terms(self):
        self.assertEqual(
            len(
                self.randomise_phenopacket.create_random_hpo_terms(
                    self.phenotypic_abnormality_children
                )
            ),
            3,
        )
        self.assertEqual(
            len(
                self.randomise_phenopacket_single_term.create_random_hpo_terms(
                    self.phenotypic_abnormality_children
                )
            ),
            3,
        )

    def test_combine_hpo_terms(self):
        self.assertEqual(
            len(
                self.randomise_phenopacket.randomise_hpo_terms(self.phenotypic_abnormality_children)
            ),
            8,
        )
        self.assertEqual(
            len(
                self.randomise_phenopacket_single_term.randomise_hpo_terms(
                    self.phenotypic_abnormality_children
                )
            ),
            4,
        )
