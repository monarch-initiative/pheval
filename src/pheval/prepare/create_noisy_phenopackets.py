import os
import pathlib
import random
import warnings

import click
from google.protobuf.json_format import MessageToJson
from oaklib.implementations.pronto.pronto_implementation import \
    ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import Family, OntologyClass, Phenopacket, PhenotypicFeature

from pheval.utils.file_utils import DirectoryFiles
from pheval.utils.phenopacket_utils import PhenopacketReader

warnings.filterwarnings("ignore")


class RandomisePhenopackets:
    """Randomises the Phenopacket phenotype."""

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
        entities = list(self.ontology.entities())
        return [x for x in entities if "HP:" in x and "HP:0034334" not in x]

    def max_real_patient_id(self) -> dict:
        phenotypic_features = list(self.phenotypic_features.items())
        if self.number_of_real_id > len(phenotypic_features):
            if len(phenotypic_features) - 2 > 0:
                self.number_of_real_id = len(phenotypic_features) - 2
            else:
                self.number_of_real_id = 1
        retained_hpo = dict(random.sample(phenotypic_features, self.number_of_real_id))
        return retained_hpo

    def change_to_parent_term(self) -> dict:
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


class RebuildPhenopackets:
    """Rebuilds the original phenopacket with the randomised phenotypes."""

    def __init__(self, phenopacket_contents, hpo_list: list, output_file: str):
        self.phenopacket_contents = phenopacket_contents
        self.hpo_list = hpo_list
        self.output_file = output_file

    def add_randomised_hpo(self):
        randomised_ppacket = Phenopacket(
            id=self.phenopacket_contents.pheno.id,
            subject=self.phenopacket_contents.pheno.subject,
            phenotypic_features=self.hpo_list,
            interpretations=self.phenopacket_contents.pheno.interpretations,
            files=self.phenopacket_contents.pheno.files,
            meta_data=self.phenopacket_contents.pheno.meta_data,
        )
        altered_phenopacket = MessageToJson(randomised_ppacket)
        return altered_phenopacket

    def write_altered_phenopacket(self):
        altered_phenopacket = self.add_randomised_hpo()
        if hasattr(self.phenopacket_contents.pheno, "proband"):
            family = Family(
                id=self.phenopacket_contents.pheno.id,
                proband=altered_phenopacket,
                pedigree=self.phenopacket_contents.pheno.pedigree,
                files=self.phenopacket_contents.pheno.files,
                meta_data=self.phenopacket_contents.pheno.meta_data,
            )
            altered_phenopacket = MessageToJson(family)
        with open(self.output_file, "w") as outfile:
            outfile.write(altered_phenopacket)
        outfile.close()


@click.command()
@click.option(
    "--phenopacket-dir",
    "-P",
    metavar="PATH",
    required=True,
    help="Path to phenopackets directory",
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
)
def create_noisy_phenopackets(
    phenopacket_dir: str,
    max_real_id: int,
    number_of_parent_terms: int,
    number_of_random_terms: int,
    output_file_suffix: str,
    output_dir: str,
):
    """Generate noisy phenopackets from existing ones."""
    try:
        os.mkdir(os.path.join(output_dir, ""))
    except FileExistsError:
        pass
    phenopackets = DirectoryFiles(phenopacket_dir, ".json").obtain_files_suffix()
    resource = OntologyResource(slug="hp.obo", local=False)
    ontology = ProntoImplementation(resource)
    for phenopacket in phenopackets:
        phenopacket_full_path = os.path.join(phenopacket_dir, phenopacket)
        phenopacket_contents = PhenopacketReader(phenopacket_full_path)
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
            pathlib.Path(phenopacket).stem
            + "-"
            + output_file_suffix
            + pathlib.Path(phenopacket).suffix,
        )
        RebuildPhenopackets(
            phenopacket_contents, new_hpo_terms, output_file
        ).write_altered_phenopacket()
