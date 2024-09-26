import random
from pathlib import Path
from typing import List, Union

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


def load_ontology(local_cached_ontology: Path = None) -> ProntoImplementation:
    """
    Load the Human Phenotype Ontology (HPO).
    Args:
        local_cached_ontology(Path): Path to the local cached ontology.
    Returns:
        ProntoImplementation: An instance of ProntoImplementation containing the loaded HPO.
    """
    if local_cached_ontology is None:
        resource = OntologyResource(slug="hp.obo", local=False)
        return ProntoImplementation(resource)
    else:
        resource = OntologyResource(slug=local_cached_ontology, local=True)
        return ProntoImplementation(resource)


class HpoRandomiser:
    """Class for randomising phenopacket phenotypic features using Human Phenotype Ontology (HPO)."""

    def __init__(self, hpo_ontology: ProntoImplementation, scramble_factor: float):
        """
        Initialise the HpoRandomiser.

        Args:
            hpo_ontology (ProntoImplementation): The instance of the HPO ontology.
            scramble_factor (float): A factor for scrambling phenotypic features.
        """
        self.hpo_ontology = hpo_ontology
        self.phenotypic_abnormalities = set(hpo_ontology.roots(predicates=["HP:0000118"]))
        self.scramble_factor = scramble_factor

    def scramble_factor_proportions(self, phenotypic_features: list[PhenotypicFeature]) -> int:
        """
        Calculate the proportion of scrambled HPO terms based on the scramble factor.

        Args:
            phenotypic_features (list[PhenotypicFeature]): List of phenotypic features.

        Returns:
            int: The calculated number of phenotypic features to be scrambled.
        """
        if len(phenotypic_features) == 1:
            return 1
        else:
            return int(round(len(phenotypic_features) * self.scramble_factor, 0))

    def retrieve_hpo_term(self, hpo_id: str) -> PhenotypicFeature:
        """
        Retrieve an HPO term based on the provided HPO ID.

        Args:
            hpo_id (str): The HPO ID of the term to retrieve.

        Returns:
            PhenotypicFeature: The PhenotypicFeature object representing the retrieved HPO term.
        """
        rels = self.hpo_ontology.entity_alias_map(hpo_id)
        hpo_term = "".join(rels[(list(rels.keys())[0])])
        return PhenotypicFeature(type=OntologyClass(id=hpo_id, label=hpo_term))

    @staticmethod
    def retain_real_patient_terms(
        phenotypic_features: List[PhenotypicFeature],
        number_of_scrambled_terms: int,
    ) -> List[PhenotypicFeature]:
        """
        Return a list of real patient HPO terms, retaining a specific number of non-scrambled terms.

        Args:
            phenotypic_features (List[PhenotypicFeature]): List of phenotypic features.
            number_of_scrambled_terms (int): The count of scrambled HPO terms.

        Returns:
            List[PhenotypicFeature]: A list of non-scrambled (real patient) HPO terms.
        """
        if len(phenotypic_features) > 1:
            number_of_real_id = len(phenotypic_features) - number_of_scrambled_terms
        else:
            number_of_real_id = 1
        return random.sample(phenotypic_features, number_of_real_id)

    def convert_patient_terms_to_parent(
        self,
        phenotypic_features: List[PhenotypicFeature],
        retained_phenotypic_features: List[PhenotypicFeature],
        number_of_scrambled_terms: int,
    ) -> List[PhenotypicFeature]:
        """
        Convert a subset of patient HPO terms to their respective parent terms.

        Args:
            phenotypic_features (List[PhenotypicFeature]): List of all phenotypic features.
            retained_phenotypic_features (List[PhenotypicFeature]): List of retained non-scrambled phenotypic features.
            number_of_scrambled_terms (int): The count of scrambled HPO terms.

        Returns:
            List[PhenotypicFeature]: A list of HPO terms converted to their parent terms.

        Note:
            This method identifies a subset of patient HPO terms that are not retained among the
            non-scrambled phenotypic features and converts them to their respective parent terms.
            It then returns a list of parent HPO terms based on the provided scrambled terms count.
            If no remaining HPO terms are available for conversion, no parent terms are returned.
        """
        remaining_hpo = [i for i in phenotypic_features if i not in retained_phenotypic_features]
        if len(remaining_hpo) == 0:
            number_of_scrambled_terms = 0
        hpo_terms_to_be_changed = list(random.sample(remaining_hpo, number_of_scrambled_terms))
        parent_terms = []
        for term in hpo_terms_to_be_changed:
            if self.hpo_ontology.label(term.type.id).startswith("obsolete"):
                obsolete_term = self.hpo_ontology.entity_metadata_map(term.type.id)
                updated_term = list(obsolete_term.values())[0][0]
                parents = self.hpo_ontology.hierarchical_parents(updated_term)
            else:
                parents = self.hpo_ontology.hierarchical_parents(term.type.id)
            if not parents:
                parent_terms.append(term)
            else:
                parent_terms.append(self.retrieve_hpo_term(random.choice(parents)))
        return parent_terms

    def create_random_hpo_terms(self, number_of_scrambled_terms: int) -> List[PhenotypicFeature]:
        """
        Generate a list of random HPO terms.

        Args:
            number_of_scrambled_terms (int): The count of random HPO terms to be generated.

        Returns:
            List[PhenotypicFeature]: A list of randomly selected HPO terms.
        """
        random_ids = list(
            random.sample(sorted(self.phenotypic_abnormalities), number_of_scrambled_terms)
        )
        return [self.retrieve_hpo_term(random_id) for random_id in random_ids]

    def randomise_hpo_terms(
        self,
        phenotypic_features: List[PhenotypicFeature],
    ) -> List[PhenotypicFeature]:
        """
        Randomise the provided phenotypic features by combining retained, parent-converted, and random HPO terms.

        Args:
            phenotypic_features (List[PhenotypicFeature]): List of phenotypic features to be randomised.

        Returns:
            List[PhenotypicFeature]: A list of randomised HPO terms.

        Note:
            This method randomises the provided phenotypic features by incorporating three types of HPO terms:
            1. Retained Patient Terms: Non-scrambled (real patient) HPO terms retained based on the scramble factor.
            2. Converted to Parent Terms: Subset of HPO terms converted to their respective parent terms.
            3. Random HPO Terms: Newly generated random HPO terms based on the scramble factor.

            The method determines the count of terms for each category and combines them to form a final list
            of randomised HPO terms to be used in the phenotypic features.
        """
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
        self,
        phenopacket: Union[Phenopacket, Family],
    ) -> Union[Phenopacket, Family]:
        """
        Randomise the phenotypic profile of a Phenopacket or Family.

        Args:
            phenopacket (Union[Phenopacket, Family]): The Phenopacket or Family to be randomised.

        Returns:
            Union[Phenopacket, Family]: The randomised Phenopacket or Family.
        """
        phenotypic_features = PhenopacketUtil(phenopacket).observed_phenotypic_features()
        random_phenotypes = self.randomise_hpo_terms(phenotypic_features)
        randomised_phenopacket = PhenopacketRebuilder(phenopacket).add_randomised_hpo(
            random_phenotypes
        )
        return randomised_phenopacket

    def create_scrambled_phenopacket(
        self,
        output_dir: Path,
        phenopacket_path: Path,
    ) -> None:
        """
        Create a scrambled version of a Phenopacket.

        Args:
            output_dir (Path): The directory to store the output scrambled Phenopacket.
            phenopacket_path (Path): The path to the original Phenopacket file.
        """
        phenopacket = phenopacket_reader(phenopacket_path)
        created_noisy_phenopacket = self.add_noise_to_phenotypic_profile(
            phenopacket,
        )
        write_phenopacket(
            created_noisy_phenopacket,
            output_dir.joinpath(phenopacket_path.name),
        )

    def create_scrambled_phenopackets(
        self,
        output_dir: Path,
        phenopacket_dir: Path,
    ) -> None:
        """
        Create scrambled versions of Phenopackets within a directory.

        Args:
            output_dir (Path): The directory to store the output scrambled Phenopackets.
            phenopacket_dir (Path): The directory containing the original Phenopacket files.
        """
        phenopacket_files = files_with_suffix(phenopacket_dir, ".json")
        for phenopacket_path in phenopacket_files:
            phenopacket = phenopacket_reader(phenopacket_path)
            created_noisy_phenopacket = self.add_noise_to_phenotypic_profile(phenopacket)
            write_phenopacket(
                created_noisy_phenopacket,
                output_dir.joinpath(
                    phenopacket_path.name,
                ),
            )


def scramble_phenopackets(
    output_dir: Path,
    phenopacket_path: Path,
    phenopacket_dir: Path,
    scramble_factor: float,
    local_cached_ontology: Path,
) -> None:
    """
    Create scrambled phenopackets from either a single phenopacket or a directory of phenopackets.

    Args:
        output_dir (Path): The directory to store the output scrambled Phenopackets.
        phenopacket_path (Path): The path to a single Phenopacket file (if applicable).
        phenopacket_dir (Path): The directory containing multiple Phenopacket files (if applicable).
        scramble_factor (float): A factor determining the level of scrambling for phenotypic features.
        local_cached_ontology (Path): The path to the local cached ontology.
    """
    output_dir.mkdir(exist_ok=True)
    ontology = load_ontology(local_cached_ontology)
    if phenopacket_path is not None:
        HpoRandomiser(ontology, scramble_factor).create_scrambled_phenopacket(
            output_dir, phenopacket_path
        )
    elif phenopacket_dir is not None:
        HpoRandomiser(ontology, scramble_factor).create_scrambled_phenopackets(
            output_dir,
            phenopacket_dir,
        )
