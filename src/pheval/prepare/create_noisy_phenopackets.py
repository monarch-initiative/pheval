import os
import random
from pathlib import Path

import click
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import OntologyClass, PhenotypicFeature

from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import PhenopacketRebuilder, PhenopacketUtil, phenopacket_reader

ontology_loaded = False


def load_ontology():
    """Loads ontology once for use in other methods."""
    global ontology_loaded

    if ontology_loaded:
        pass

    ontology_loaded = True
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


class RandomisePhenopackets:
    """Randomises the Phenopacket phenotypic profile."""

    def __init__(
        self,
        ontology,
        phenotypic_features: list[PhenotypicFeature],
        number_of_real_id: int,
        number_of_changed_terms: int,
        number_of_random_terms: int,
    ):
        self.ontology = ontology
        self.phenotypic_features = phenotypic_features
        self.number_of_real_id = number_of_real_id
        self.number_of_changed_terms = number_of_changed_terms
        self.number_of_random_terms = number_of_random_terms

    def create_clean_entities(self) -> list:
        """Returns a list of HPO terms."""
        entities = list(self.ontology.entities())
        return [x for x in entities if "HP:" in x and "HP:0034334" not in x]

    def retrieve_hpo_label(self, hpo_id: str) -> PhenotypicFeature:
        rels = self.ontology.entity_alias_map(hpo_id)
        hpo_term = "".join(rels[(list(rels.keys())[0])])
        return PhenotypicFeature(type=OntologyClass(id=hpo_id, label=hpo_term))

    def retain_real_patient_terms(self) -> list[PhenotypicFeature]:
        """Returns a dictionary of the maximum number of real patient HPO terms."""
        if self.number_of_real_id > len(self.phenotypic_features):
            if len(self.phenotypic_features) - 2 > 0:
                self.number_of_real_id = len(self.phenotypic_features) - 2
            else:
                self.number_of_real_id = 1
        return random.sample(self.phenotypic_features, self.number_of_real_id)

    def convert_patient_terms_to_parent(self) -> list[PhenotypicFeature]:
        """Returns a dictionary of the HPO terms that have been converted to a parent term."""
        retained_hpo = self.retain_real_patient_terms()
        remaining_hpo = [i for i in self.phenotypic_features if i not in retained_hpo]
        if self.number_of_changed_terms > len(remaining_hpo):
            self.number_of_changed_terms = len(remaining_hpo)
        hpo_terms_to_be_changed = random.sample(remaining_hpo, self.number_of_changed_terms)
        parent_terms = []
        for term in hpo_terms_to_be_changed:
            try:
                parent_terms.append(
                    self.retrieve_hpo_label(self.ontology.hierararchical_parents(term.type.id)[0])
                )
            except IndexError:
                obsolete_term = self.ontology.entity_metadata_map(term.type.id)
                updated_term = list(obsolete_term.values())[0][0]
                parent_terms.append(
                    self.retrieve_hpo_label(self.ontology.hierararchical_parents(updated_term)[0])
                )
        return parent_terms

    def create_random_hpo_terms(self) -> list[PhenotypicFeature]:
        """Returns a dictionary of random HPO terms"""
        clean_entities = self.create_clean_entities()
        random_ids = list(random.sample(clean_entities, self.number_of_random_terms))
        return [self.retrieve_hpo_label(random_id) for random_id in random_ids]

    def randomise_hpo_terms(self) -> list[PhenotypicFeature]:
        """Combines real patient HPO terms, parent terms and randomised terms."""
        return (
            self.retain_real_patient_terms()
            + self.convert_patient_terms_to_parent()
            + self.create_random_hpo_terms()
        )


def noisy_phenopacket(
    phenopacket: Path,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    output_file_suffix: str,
    output_dir: Path,
):
    """Randomises a single phenopacket phenotypic profile, writing to a new .json file"""
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    ontology = load_ontology()
    phenopacket_contents = phenopacket_reader(phenopacket)
    phenotypic_features = PhenopacketUtil(
        phenopacket_contents
    ).remove_excluded_phenotypic_features()
    randomised_hpo_terms = RandomisePhenopackets(
        ontology,
        phenotypic_features,
        max_real_id,
        number_of_parent_terms,
        number_of_random_terms,
    ).randomise_hpo_terms()
    output_file = os.path.join(
        output_dir,
        Path(phenopacket).stem + "-" + output_file_suffix + Path(phenopacket).suffix,
    )
    phenopacket_rebuilder = PhenopacketRebuilder(phenopacket_contents)
    phenopacket_rebuilder.add_randomised_hpo(randomised_hpo_terms)
    phenopacket_rebuilder.write_phenopacket(Path(output_file))


@click.command()
@click.option(
    "--phenopacket",
    "-P",
    metavar="PATH",
    required=True,
    help="Path to phenopacket to be randomised",
    type=Path,
)
@click.option(
    "--max-real-id",
    "-m",
    metavar="<int>",
    required=True,
    help="Maximum number of real patient HPO ids to retain",
    type=int,
    default=3,
    show_default=True,
)
@click.option(
    "--number-of-parent-terms",
    "-p",
    metavar="<int>",
    required=True,
    help="Number of real patient HPO ids to change to parent terms",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--number-of-random-terms",
    "-r",
    metavar="<int>",
    required=True,
    help="Number of random HPO ids to introduce",
    type=int,
    default=3,
    show_default=True,
)
@click.option(
    "--output-file-suffix",
    "-o",
    metavar="<str>",
    required=True,
    help="Suffix to append to output file",
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="noisy_phenopackets",
    type=Path,
)
def create_noisy_phenopacket(
    phenopacket: Path,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    output_file_suffix: str,
    output_dir: Path,
):
    """Generate a noisy phenopacket from an existing one."""
    noisy_phenopacket(
        phenopacket,
        max_real_id,
        number_of_parent_terms,
        number_of_random_terms,
        output_file_suffix,
        output_dir,
    )


@click.command()
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    required=True,
    help="Path to phenopackets directory",
    type=Path,
)
@click.option(
    "--max-real-id",
    "-m",
    metavar="<int>",
    required=True,
    help="Maximum number of real patient HPO ids to retain",
    type=int,
    default=3,
    show_default=True,
)
@click.option(
    "--number-of-parent-terms",
    "-p",
    metavar="<int>",
    required=True,
    help="Number of real patient HPO ids to change to parent terms",
    type=int,
    default=2,
    show_default=True,
)
@click.option(
    "--number-of-random-terms",
    "-r",
    metavar="<int>",
    required=True,
    help="Number of random HPO ids to introduce",
    type=int,
    default=3,
    show_default=True,
)
@click.option(
    "--output-file-suffix",
    "-o",
    metavar="<str>",
    required=True,
    help="Suffix to append to output file",
)
@click.option(
    "--output-dir",
    "-O",
    metavar="PATH",
    required=True,
    help="Path for creation of output directory",
    default="noisy_phenopackets",
    type=Path,
)
def create_noisy_phenopackets(
    phenopacket_dir: Path,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    output_file_suffix: str,
    output_dir: Path,
):
    """Generate noisy phenopackets from existing ones."""
    phenopackets = files_with_suffix(phenopacket_dir, ".json")
    for phenopacket in phenopackets:
        noisy_phenopacket(
            phenopacket,
            max_real_id,
            number_of_parent_terms,
            number_of_random_terms,
            output_file_suffix,
            output_dir,
        )
