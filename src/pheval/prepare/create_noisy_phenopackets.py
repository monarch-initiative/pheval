import random
from pathlib import Path

import click
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
    """Loads ontology once for use in other methods."""
    resource = OntologyResource(slug="hp.obo", local=False)
    return ProntoImplementation(resource)


class HpoRandomiser:
    def __init__(self, hpo_ontology):
        self.hpo_ontology = hpo_ontology
        self.phenotypic_abnormalities = set(hpo_ontology.roots(predicates=["HP:0000118"]))

    def retrieve_hpo_term(self, hpo_id: str) -> PhenotypicFeature:
        """Retrieves term for hpo id."""
        rels = self.hpo_ontology.entity_alias_map(hpo_id)
        hpo_term = "".join(rels[(list(rels.keys())[0])])
        return PhenotypicFeature(type=OntologyClass(id=hpo_id, label=hpo_term))

    @staticmethod
    def retain_real_patient_terms(
        phenotypic_features: list[PhenotypicFeature], number_of_real_id
    ) -> list[PhenotypicFeature]:
        """Returns a list of the maximum number of real patient HPO terms."""
        if number_of_real_id > len(phenotypic_features):
            if len(phenotypic_features) - 2 > 0:
                number_of_real_id = len(phenotypic_features) - 2
            else:
                number_of_real_id = 1
        return random.sample(phenotypic_features, number_of_real_id)

    def convert_patient_terms_to_parent(
        self,
        phenotypic_features: list[PhenotypicFeature],
        retained_phenotypic_features: list[PhenotypicFeature],
        number_of_parent_terms: int,
    ) -> list[PhenotypicFeature]:
        """Returns a list of the HPO terms that have been converted to a parent term."""
        remaining_hpo = [
            i for i in phenotypic_features if i not in retained_phenotypic_features
        ]
        number_of_parent_terms = (
            len(remaining_hpo)
            if number_of_parent_terms > len(remaining_hpo)
            else number_of_parent_terms
        )
        hpo_terms_to_be_changed = random.sample(remaining_hpo, number_of_parent_terms)
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

    def create_random_hpo_terms(self, number_of_random_terms: int) -> list[PhenotypicFeature]:
        """Returns a list of random HPO terms"""
        random_ids = list(random.sample(sorted(self.phenotypic_abnormalities), number_of_random_terms))
        return [self.retrieve_hpo_term(random_id) for random_id in random_ids]

    def randomise_hpo_terms(
        self,
        phenotypic_features: list[PhenotypicFeature],
        max_real_id: int,
        number_of_parent_terms: int,
        number_of_random_terms: int,
    ) -> list[PhenotypicFeature]:
        retained_patient_terms = self.retain_real_patient_terms(phenotypic_features, max_real_id)
        return (
            retained_patient_terms
            + self.convert_patient_terms_to_parent(
                phenotypic_features, retained_patient_terms, number_of_parent_terms
            )
            + self.create_random_hpo_terms(number_of_random_terms)
        )


def noisy_phenopacket(
    hpo_randomiser: HpoRandomiser,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    phenopacket: Phenopacket or Family,
) -> Phenopacket or Family:
    # phenopacket_util = PhenopacketUtil(phenopacket)
    # phenotypic_features = phenopacket_util.observed_phenotypic_features()
    phenotypic_features = PhenopacketUtil(
        phenopacket
    ).observed_phenotypic_features()
    random_phenotypes = hpo_randomiser.randomise_hpo_terms(
        phenotypic_features, max_real_id, number_of_parent_terms, number_of_random_terms
    )
    randomised_phenopacket = PhenopacketRebuilder(phenopacket).add_randomised_hpo(random_phenotypes)
    return randomised_phenopacket


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
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    ontology = load_ontology()
    hpo_randomiser = HpoRandomiser(ontology)
    phenopacket_contents = phenopacket_reader(phenopacket)
    created_noisy_phenopacket = noisy_phenopacket(
        hpo_randomiser, max_real_id, number_of_parent_terms, number_of_random_terms, phenopacket_contents
    )
    write_phenopacket(
        created_noisy_phenopacket,
        output_dir.joinpath(
            Path(phenopacket).stem + "-" + output_file_suffix + Path(phenopacket).suffix
        ),
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
    try:
        output_dir.mkdir()
    except FileExistsError:
        pass
    ontology = load_ontology()
    hpo_randomiser = HpoRandomiser(ontology)
    phenopacket_files = files_with_suffix(phenopacket_dir, ".json")
    for phenopacket_file in phenopacket_files:
        phenopacket = phenopacket_reader(phenopacket_file)
        created_noisy_phenopacket = noisy_phenopacket(
            hpo_randomiser, max_real_id, number_of_parent_terms, number_of_random_terms, phenopacket
        )
        write_phenopacket(
            created_noisy_phenopacket,
            output_dir.joinpath(
                Path(phenopacket_file).stem
                + "-"
                + output_file_suffix
                + Path(phenopacket_file).suffix,
            ),
        )
