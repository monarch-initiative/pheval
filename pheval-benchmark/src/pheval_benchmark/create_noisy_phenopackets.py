import click
from oaklib.resource import OntologyResource
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from phenopackets import Phenopacket, PhenotypicFeature, OntologyClass
from google.protobuf.json_format import Parse, MessageToJson
from pheval_benchmark.assess_prioritisation import directory_files
import json
import random
import pathlib
import os
import warnings

warnings.filterwarnings("ignore")
path_to_obo = os.path.dirname(os.path.realpath(__file__)) + "/resources/obo/hp.obo"
resource = OntologyResource(slug=path_to_obo, local=True)
oi = ProntoImplementation(resource)


def separate_phenotypic_features(phenopacket):
    with open(phenopacket) as phenopacket_:
        ppacket = json.load(phenopacket_)
        ppacket = Parse(json.dumps(ppacket), Phenopacket())
        all_phenotypic_features = ppacket.phenotypic_features
    phenopacket_.close()
    return all_phenotypic_features


def separate_excluded_phenotypes(all_phenotypic_features) -> dict:
    phenotypic_features = {}
    for p in all_phenotypic_features:
        if p.excluded:
            continue
        phenotypic_features[p.type.id] = p.type.label
    return phenotypic_features


def max_real_patient_id(phenotypic_features: dict, number_of_real_id: int) -> dict:
    phenotypic_features = list(phenotypic_features.items())
    if number_of_real_id > len(phenotypic_features):
        if len(phenotypic_features) - 2 > 0:
            number_of_real_id = len(phenotypic_features) - 2
        else:
            number_of_real_id = 1
    retained_hpo = dict(random.sample(phenotypic_features, number_of_real_id))
    return retained_hpo


def change_to_parent_term(phenotypic_features: dict, retained_hpo: dict, number_of_changed_terms: int) -> dict:
    remaining_hpo = list(phenotypic_features.items() ^ retained_hpo.items())
    if number_of_changed_terms > len(remaining_hpo):
        number_of_changed_terms = len(remaining_hpo)
    parent_hpo = list(random.sample(remaining_hpo, number_of_changed_terms))
    parent = {}
    for p in parent_hpo:
        try:
            parent_id = oi.hierararchical_parents(p[0])[0]
            rels = oi.alias_map_by_curie(parent_id)
            parent_term = rels[(list(rels.keys())[0])]
            parent_term = ''.join(parent_term)
            parent[parent_id] = parent_term
        except IndexError:
            pass
    return parent


def random_hpo_terms(number_of_random_terms: int) -> dict:
    all_id = []
    with open("src/pheval_benchmark/resources/obo/hp.obo") as f:
        for line in f:
            if line.startswith("id"):
                line = line.split(": ", 1)[1].rstrip()
                all_id.append(line)
    f.close()
    random_id = list(random.sample(all_id, number_of_random_terms))
    random_id_dict = {}
    for r in random_id:
        rels = oi.alias_map_by_curie(r)
        random_term = rels[(list(rels.keys())[0])]
        random_term = ''.join(random_term)
        random_id_dict[r] = random_term
    return random_id_dict


def combine_hpo_terms(retained_hpo: dict, parent: dict, random_id_dict: dict) -> list:
    retained_hpo.update(parent)
    retained_hpo.update(random_id_dict)
    hpo_list = []
    for key, value in retained_hpo.items():
        new_hpo = PhenotypicFeature(type=OntologyClass(id=key, label=value))
        hpo_list.append(new_hpo)
    return hpo_list


def rebuild_phenopacket(phenopacket, hpo_list: list, output_file: str):
    with open(phenopacket) as phenopacket_:
        ppacket = json.load(phenopacket_)
        ppacket = Parse(json.dumps(ppacket), Phenopacket())
        phenopacket = Phenopacket(id=ppacket.id, subject=ppacket.subject,
                                                   phenotypic_features=hpo_list,
                                                   interpretations=ppacket.interpretations, files=ppacket.files,
                                                   meta_data=ppacket.meta_data)
    phenopacket_.close()
    altered_phenopacket = MessageToJson(phenopacket)
    with open(output_file, 'w') as outfile:
        outfile.write(altered_phenopacket)
    outfile.close()


@click.command()
@click.option("--phenopacket-dir", "-P", metavar='PATH', required=True, help="Path to phenopackets directory")
@click.option("--max-real-id", "-m", metavar='<int>', required=True,
              help="Maximum number of real patient HPO ids to retain", type=int, default=3, show_default=True)
@click.option("--number-of-parent-terms", "-p", metavar='<int>', required=True,
              help="Number of real patient HPO ids to change to parent terms", type=int, default=2, show_default=True)
@click.option("--number-of-random-terms", "-r", metavar='<int>', required=True,
              help="Number of random HPO ids to introduce", type=int, default=3, show_default=True)
@click.option("--output-file-suffix", "-o", metavar='<str>', required=True,
              help="Suffix to append to output file")
@click.option("--output-dir", "-O", metavar="PATH", required=True, help="Path for creation of output directory",
              default="noisy_phenopackets")
def create_noisy_phenopackets(phenopacket_dir: str, max_real_id: int, number_of_parent_terms: int,
                              number_of_random_terms: int, output_file_suffix: str, output_dir: str):
    """ Generate noisy phenopackets from existing ones. """
    try:
        os.mkdir(os.path.join(output_dir, ''))
    except FileExistsError:
        pass
    phenopackets = directory_files(phenopacket_dir)
    for phenopacket in phenopackets:
        phenopacket_full_path = os.path.join(phenopacket_dir, phenopacket)
        all_phenotypic_features = separate_phenotypic_features(phenopacket_full_path)
        phenotypic_features = separate_excluded_phenotypes(all_phenotypic_features)
        retained_id = max_real_patient_id(phenotypic_features, max_real_id)
        parent = change_to_parent_term(phenotypic_features, retained_id, number_of_parent_terms)
        random_hpo = random_hpo_terms(number_of_random_terms)
        new_hpo_terms = combine_hpo_terms(retained_id, parent, random_hpo)
        output_file = os.path.join(output_dir, pathlib.Path(phenopacket).stem + "-" + output_file_suffix + pathlib.Path(
            phenopacket).suffix)
        rebuild_phenopacket(phenopacket_full_path, new_hpo_terms, output_file)
