import unittest

from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource

from pheval.prepare.create_noisy_phenopackets import RandomisePhenopackets


class TestRandomisePhenopackets(unittest.TestCase):
    phenotypic_features_single_term = {"HP:0000965": "Cutis marmorata"}
    phenotypic_features = {
        "HP:0000256": "Macrocephaly",
        "HP:0002059": "Cerebral atrophy",
        "HP:0100309": "Subdural hemorrhage",
        "HP:0003150": "Glutaric aciduria",
        "HP:0001332": "Dystonia",
    }

    @classmethod
    def setUpClass(cls) -> None:
        resource = OntologyResource(slug="hp.obo", local=False)
        ontology = ProntoImplementation(resource)
        cls.randomise_phenopacket = RandomisePhenopackets(
            ontology, cls.phenotypic_features, 3, 2, 3
        )
        cls.randomise_phenopacket_single_term = RandomisePhenopackets(
            ontology, cls.phenotypic_features_single_term, 3, 2, 3
        )

    def test_create_clean_entities(self):
        clean_entities = self.randomise_phenopacket.create_clean_entities()
        self.assertTrue("HPO:" in term for term in clean_entities)
        self.assertFalse("HP:0034334" in clean_entities)

    def test_max_real_patient_id(self):
        self.assertTrue(len(self.randomise_phenopacket.max_real_patient_id()), 3)
        self.assertTrue(
            all(
                hpo_term in self.phenotypic_features.keys()
                for hpo_term in self.randomise_phenopacket.max_real_patient_id().keys()
            )
        )
        self.assertTrue(len(self.randomise_phenopacket_single_term.max_real_patient_id()), 1)
        self.assertEqual(
            self.phenotypic_features_single_term.keys(),
            self.randomise_phenopacket_single_term.max_real_patient_id().keys(),
        )

    def test_change_to_parent_term(self):
        self.assertTrue(len(self.randomise_phenopacket.change_to_parent_term()), 2)
        self.assertFalse(
            all(
                hpo_term in self.phenotypic_features.keys()
                for hpo_term in self.randomise_phenopacket.change_to_parent_term().keys()
            )
        )
        self.assertEqual(self.randomise_phenopacket_single_term.change_to_parent_term(), {})

    def test_random_hpo_terms(self):
        self.assertEqual(len(self.randomise_phenopacket.random_hpo_terms()), 3)
        self.assertEqual(len(self.randomise_phenopacket_single_term.random_hpo_terms()), 3)

    def test_combine_hpo_terms(self):
        self.assertEqual(len(self.randomise_phenopacket.combine_hpo_terms()), 8)
        self.assertEqual(len(self.randomise_phenopacket_single_term.combine_hpo_terms()), 4)
