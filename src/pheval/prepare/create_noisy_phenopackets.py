import random
from pathlib import Path

from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import Family, OntologyClass, Phenopacket, PhenotypicFeature

from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import (
    PhenopacketRebuilder,
    PhenopacketUtil,
    phenopacket_reader,
    write_phenopacket,
)


def load_ontology():
    """Loads human phenotype ontology."""
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


class HpoRandomiser:
    """Randomises phenopacket phenotypic features."""

    def __init__(self, hpo_ontology, scramble_factor: float):
        self.hpo_ontology = hpo_ontology
        self.phenotypic_abnormalities = set(hpo_ontology.roots(predicates=["HP:0000118"]))
        self.scramble_factor = scramble_factor

    def scramble_factor_proportions(self, phenotypic_features: list[PhenotypicFeature]):
        """Calculate proportion of scrambled hpo terms from scramble factor."""
        if len(phenotypic_features) == 1:
            return 1
        else:
            return int(round(len(phenotypic_features) * self.scramble_factor, 0))

    def retrieve_hpo_term(self, hpo_id: str) -> PhenotypicFeature:
        """Retrieves term for hpo id."""
        rels = self.hpo_ontology.entity_alias_map(hpo_id)
        hpo_term = "".join(rels[(list(rels.keys())[0])])
        return PhenotypicFeature(type=OntologyClass(id=hpo_id, label=hpo_term))

    @staticmethod
    def retain_real_patient_terms(
        phenotypic_features: list[PhenotypicFeature],
        number_of_scrambled_terms: int,
    ) -> list[PhenotypicFeature]:
        """Returns a list of the maximum number of real patient HPO terms."""
        if len(phenotypic_features) > 1:
            number_of_real_id = len(phenotypic_features) - number_of_scrambled_terms
        else:
            number_of_real_id = 1
        return random.sample(phenotypic_features, number_of_real_id)

    def convert_patient_terms_to_parent(
        self,
        phenotypic_features: list[PhenotypicFeature],
        retained_phenotypic_features: list[PhenotypicFeature],
        number_of_scrambled_terms: int,
    ) -> list[PhenotypicFeature]:
        """Returns a list of the HPO terms that have been converted to a parent term."""
        remaining_hpo = [i for i in phenotypic_features if i not in retained_phenotypic_features]
        if len(remaining_hpo) == 0:
            number_of_scrambled_terms = 0
        hpo_terms_to_be_changed = list(random.sample(remaining_hpo, number_of_scrambled_terms))
        parent_terms = []
        for term in hpo_terms_to_be_changed:
            try:
                parent_terms.append(
                    self.retrieve_hpo_term(
                        self.hpo_ontology.hierararchical_parents(term.type.id)[0]
                    )
                )
            except IndexError:
                obsolete_term = self.hpo_ontology.entity_metadata_map(term.type.id)
                updated_term = list(obsolete_term.values())[0][0]
                parent_terms.append(
                    self.retrieve_hpo_term(
                        self.hpo_ontology.hierararchical_parents(updated_term)[0]
                    )
                )
        return parent_terms

    def create_random_hpo_terms(self, number_of_scrambled_terms: int) -> list[PhenotypicFeature]:
        """Returns a list of random HPO terms"""
        random_ids = list(
            random.sample(sorted(self.phenotypic_abnormalities), number_of_scrambled_terms)
        )
        return [self.retrieve_hpo_term(random_id) for random_id in random_ids]

    def randomise_hpo_terms(
        self,
        phenotypic_features: list[PhenotypicFeature],
    ) -> list[PhenotypicFeature]:
        """Returns a list of randomised HPO terms."""
        number_of_scrambled_terms = self.scramble_factor_proportions(phenotypic_features)
        retained_patient_terms = self.retain_real_patient_terms(
            phenotypic_features, number_of_scrambled_terms
        )
        return (
            retained_patient_terms
            + self.convert_patient_terms_to_parent(
                phenotypic_features, retained_patient_terms, number_of_scrambled_terms
            )
            + self.create_random_hpo_terms(number_of_scrambled_terms)
        )


def add_noise_to_phenotypic_profile(
    hpo_randomiser: HpoRandomiser,
    phenopacket: Phenopacket or Family,
) -> Phenopacket or Family:
    """Randomises the phenotypic profile of a phenopacket."""
    # phenopacket_util = PhenopacketUtil(phenopacket)
    # phenotypic_features = phenopacket_util.observed_phenotypic_features()
    phenotypic_features = PhenopacketUtil(phenopacket).observed_phenotypic_features()
    random_phenotypes = hpo_randomiser.randomise_hpo_terms(phenotypic_features)
    randomised_phenopacket = PhenopacketRebuilder(phenopacket).add_randomised_hpo(random_phenotypes)
    return randomised_phenopacket


def create_scrambled_phenopacket(
    output_dir: Path, phenopacket_path: Path, scramble_factor: float
) -> None:
    """Creates a scrambled phenopacket."""
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    ontology = load_ontology()
    hpo_randomiser = HpoRandomiser(ontology, scramble_factor)
    phenopacket = phenopacket_reader(phenopacket_path)
    created_noisy_phenopacket = add_noise_to_phenotypic_profile(
        hpo_randomiser,
        phenopacket,
    )
    write_phenopacket(
        created_noisy_phenopacket,
        output_dir.joinpath(phenopacket_path.name),
    )


def create_scrambled_phenopackets(
    output_dir: Path, phenopacket_dir: Path, scramble_factor: float
) -> None:
    """Creates scrambled phenopackets."""
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    ontology = load_ontology()
    hpo_randomiser = HpoRandomiser(ontology, scramble_factor)
    phenopacket_files = files_with_suffix(phenopacket_dir, ".json")
    for phenopacket_path in phenopacket_files:
        phenopacket = phenopacket_reader(phenopacket_path)
        created_noisy_phenopacket = add_noise_to_phenotypic_profile(hpo_randomiser, phenopacket)
        write_phenopacket(
            created_noisy_phenopacket,
            output_dir.joinpath(
                phenopacket_path.name,
            ),
        )


def scramble_phenopackets(
    output_dir: Path, phenopacket_path: Path, phenopacket_dir: Path, scramble_factor: float
) -> None:
    """Create scrambled phenopackets from either a single phenopacket or directory of phenopackets."""
    if phenopacket_path is not None:
        create_scrambled_phenopacket(output_dir, phenopacket_path, scramble_factor)
    elif phenopacket_dir is not None:
        create_scrambled_phenopackets(output_dir, phenopacket_dir, scramble_factor)
