import os
import random
from pathlib import Path

import click
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import OntologyClass, PhenotypicFeature

from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import (
    PhenopacketReader,
    PhenopacketRebuilder,
    PhenopacketWriter,
)

ontology_loaded = False


def load_ontology():
    """Loads ontology once for use in other methods."""
    global ontology_loaded

    if ontology_loaded:
        pass

    ontology_loaded = True
    resource = OntologyResource(slug="hp.obo", local=False)
    ontology = ProntoImplementation(resource)
    return ontology


class RandomisePhenopackets:
    """Randomises the Phenopacket phenotypic profile."""

    def __init__(
        self,
        ontology,
        phenotypic_features: dict,
        number_of_real_id: int,
        number_of_changed_terms: int,
        number_of_random_terms: int,
    ):
        self.ontology = ontology
        self.phenotypic_features = phenotypic_features
        self.number_of_real_id = number_of_real_id
        self.number_of_changed_terms = number_of_changed_terms
        self.number_of_random_terms = number_of_random_terms

    def create_clean_entities(self):
        """Returns a list of HPO terms."""
        entities = list(self.ontology.entities())
        return [x for x in entities if "HP:" in x and "HP:0034334" not in x]

    def max_real_patient_id(self) -> dict:
        """Returns a dictionary of the maximum number of real patient HPO terms."""
        phenotypic_features = list(self.phenotypic_features.items())
        if self.number_of_real_id > len(phenotypic_features):
            if len(phenotypic_features) - 2 > 0:
                self.number_of_real_id = len(phenotypic_features) - 2
            else:
                self.number_of_real_id = 1
        retained_hpo = dict(random.sample(phenotypic_features, self.number_of_real_id))
        return retained_hpo

    def change_to_parent_term(self) -> dict:
        """Returns a dictionary of the HPO terms that have been changed to a parent term."""
        retained_hpo = self.max_real_patient_id()
        remaining_hpo = list(self.phenotypic_features.items() ^ retained_hpo.items())
        if self.number_of_changed_terms > len(remaining_hpo):
            self.number_of_changed_terms = len(remaining_hpo)
        parent_hpo = list(random.sample(remaining_hpo, self.number_of_changed_terms))
        parent = {}
        for p in parent_hpo:
            try:
                parent_id = self.ontology.hierararchical_parents(p[0])[0]
                rels = self.ontology.entity_alias_map(parent_id)
                parent_term = rels[(list(rels.keys())[0])]
                parent_term = "".join(parent_term)
                parent[parent_id] = parent_term
            except IndexError:
                obsolete = self.ontology.entity_metadata_map(p[0])
                updated_term = list(obsolete.values())[0][0]
                parent_id = self.ontology.hierararchical_parents(updated_term)[0]
                rels = self.ontology.entity_alias_map(parent_id)
                parent_term = rels[(list(rels.keys())[0])]
                parent_term = "".join(parent_term)
                parent[parent_id] = parent_term
        return parent

    def random_hpo_terms(self) -> dict:
        """Returns a dictionary of random HPO terms"""
        clean_entities = self.create_clean_entities()
        random_id = list(random.sample(clean_entities, self.number_of_random_terms))
        random_id_dict = {}
        for r in random_id:
            rels = self.ontology.entity_alias_map(r)
            random_term = rels[(list(rels.keys())[0])]
            random_term = "".join(random_term)
            random_id_dict[r] = random_term
        return random_id_dict

    def combine_hpo_terms(self) -> list:
        """Combines real patient HPO terms, parent terms and randomised terms."""
        retained_hpo = self.max_real_patient_id()
        parent = self.change_to_parent_term()
        random_id_dict = self.random_hpo_terms()
        retained_hpo.update(parent)
        retained_hpo.update(random_id_dict)
        hpo_list = []
        for key, value in retained_hpo.items():
            new_hpo = PhenotypicFeature(type=OntologyClass(id=key, label=value))
            hpo_list.append(new_hpo)
        return hpo_list


def noisy_phenopacket(
    phenopacket: Path,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    output_file_suffix: str,
    output_dir: Path,
    ontology,
):
    """Randomises a single phenopacket phenotypic profile, writing to a new .json file"""
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    phenopacket_contents = PhenopacketReader(phenopacket)
    phenotypic_features = phenopacket_contents.remove_excluded_phenotypic_features()
    new_hpo_terms = RandomisePhenopackets(
        ontology,
        phenotypic_features,
        max_real_id,
        number_of_parent_terms,
        number_of_random_terms,
    ).combine_hpo_terms()
    output_file = os.path.join(
        output_dir,
        Path(phenopacket).stem + "-" + output_file_suffix + Path(phenopacket).suffix,
    )
    altered_phenopacket = PhenopacketRebuilder(phenopacket_contents).add_randomised_hpo(
        new_hpo_terms
    )
    altered_phenopacket_message = PhenopacketRebuilder(phenopacket_contents).create_json_message(
        altered_phenopacket
    )
    PhenopacketWriter(altered_phenopacket_message, Path(output_file)).write_file()


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
    ontology = load_ontology()
    noisy_phenopacket(
        phenopacket,
        max_real_id,
        number_of_parent_terms,
        number_of_random_terms,
        output_file_suffix,
        output_dir,
        ontology,
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
    ontology = load_ontology()
    phenopackets = files_with_suffix(phenopacket_dir, ".json")
    for phenopacket in phenopackets:
        noisy_phenopacket(
            phenopacket,
            max_real_id,
            number_of_parent_terms,
            number_of_random_terms,
            output_file_suffix,
            output_dir,
            ontology,
        )
