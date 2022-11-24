import os

import click
import yaml
from google.protobuf.json_format import MessageToJson
from google.protobuf.timestamp_pb2 import Timestamp
from oaklib.implementations.pronto.pronto_implementation import \
    ProntoImplementation
from oaklib.resource import OntologyResource
from phenopackets import (Diagnosis, Family, File, GeneDescriptor,
                          GenomicInterpretation, Individual, Interpretation,
                          MetaData, OntologyClass, Pedigree, Phenopacket,
                          PhenotypicFeature, Resource, VariantInterpretation,
                          VariationDescriptor, VcfRecord)
from pheval_benchmark.utils.utils import DirectoryFiles

path_to_obo = os.path.dirname(os.path.realpath(__file__)).replace(
    "/exomiser", "/resources/obo/hp.obo"
)
resource = OntologyResource(slug=path_to_obo, local=True)
oi = ProntoImplementation(resource)
path_to_geno = os.path.dirname(os.path.realpath(__file__)) + "/geno.owl"
genotype_resource = OntologyResource(slug=path_to_geno, local=True)
go_oi = ProntoImplementation(genotype_resource)


class YamlToFamilyPhenopacketConversion:
    def __init__(self, file, diagnoses):
        with open(file) as yaml_job_file:
            self.job_file = yaml.load(yaml_job_file, Loader=yaml.FullLoader)
        yaml_job_file.close()
        self.output_file = file.replace("yml", "json")
        with open(
            os.path.dirname(os.path.realpath(__file__)) + "/hgnc_complete_set.txt"
        ) as gene_set:
            self.gene_id_symbol = {}
            next(gene_set)
            for line in gene_set:
                l = line.split("\t")
                self.gene_id_symbol[l[1]] = l[0]
            gene_set.close()
        self.diagnoses = diagnoses

    def construct_subject_ped(self):
        ped_file_path = self.job_file["analysis"]["ped"]
        with open(ped_file_path) as ped:
            first_line = ped.readline()
            first_line = first_line.split("\t")
            family_id = first_line[0]
            maternal_id = first_line[3]
            paternal_id = first_line[2]
            persons = []
            ped.seek(0)
            for ped_line in ped:
                p = ped_line.split("\t")
                if int(p[4]) == 1:
                    sex = 2
                if int(p[4]) == 2:
                    sex = 1
                if p[1] == self.job_file["analysis"]["proband"]:
                    person = Pedigree.Person(
                        family_id=family_id,
                        individual_id=p[1],
                        paternal_id=paternal_id,
                        maternal_id=maternal_id,
                        sex=sex,
                        affected_status=int(p[5]),
                    )
                    subject = Individual(
                        id=self.job_file["analysis"]["proband"], sex=sex
                    )
                    persons.append(person)
                else:
                    person = Pedigree.Person(
                        family_id=family_id,
                        individual_id=p[1],
                        sex=sex,
                        affected_status=int(p[5]),
                    )
                    persons.append(person)
        ped.close()
        pedigree = Pedigree(persons=persons)
        return subject, pedigree, family_id

    def create_phenotypic_interpretations(self):
        ids = self.job_file["analysis"]["hpoIds"]
        phenotypic_features = []
        for i in ids:
            try:
                rels = oi.alias_map_by_curie(i)
                term = rels[(list(rels.keys())[0])]
                term = "".join(term)
                hpo = PhenotypicFeature(type=OntologyClass(id=i, label=term))
                phenotypic_features.append(hpo)
            except AttributeError:
                hpo = PhenotypicFeature(type=OntologyClass(id=i))
                phenotypic_features.append(hpo)
        return phenotypic_features

    def create_interpretations(self):
        int_id = self.job_file["analysis"]["proband"] + "-interpretation"
        interpretations = []
        with open(self.diagnoses) as di:  # the diagnosed genes file
            genomic_interpretations = []
            for line in di:
                line = line.strip("\n")
                l = line.split("\t")
                if l[1] == self.job_file["analysis"]["proband"]:
                    if l[14] == "-":
                        interpretation_status = 0
                    if l[14] == "Full":
                        interpretation_status = 3
                    if l[14] == "Partial":
                        interpretation_status = 3
                    try:
                        vcf_record = VcfRecord(
                            genome_assembly=self.job_file["analysis"]["genomeAssembly"],
                            chrom=l[3],
                            pos=int(l[4].replace(",", "")),
                            ref=l[6].split("/")[0],
                            alt=l[6].split("/")[1],
                        )
                        try:
                            gene_context = GeneDescriptor(
                                value_id=self.gene_id_symbol[l[8]], symbol=l[8]
                            )
                        except KeyError:
                            gene_context = GeneDescriptor(symbol=l[8])
                        allelic_state = OntologyClass(
                            id=list(go_oi.basic_search(l[10].lower()))[0],
                            label=l[10].lower(),
                        )
                        variation_descriptor = VariationDescriptor(
                            id=self.job_file["analysis"]["proband"]
                            + ":"
                            + l[3]
                            + ":"
                            + l[4],
                            gene_context=gene_context,
                            vcf_record=vcf_record,
                            allelic_state=allelic_state,
                        )
                        variant_interpretation = VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=variation_descriptor,
                        )
                        genomic_interpretation = GenomicInterpretation(
                            subject_or_biosample_id=self.job_file["analysis"][
                                "proband"
                            ],
                            interpretation_status=interpretation_status,
                            variant_interpretation=variant_interpretation,
                        )
                        genomic_interpretations.append(genomic_interpretation)
                        diagnosis = Diagnosis(
                            genomic_interpretations=genomic_interpretations
                        )
                        interpretation = Interpretation(
                            id=int_id, progress_status="SOLVED", diagnosis=diagnosis
                        )  # All DDD are solved
                    except Exception:
                        vcf_record = VcfRecord(
                            genome_assembly=self.job_file["analysis"]["genomeAssembly"],
                            chrom=l[3],
                            pos=int(l[4].replace(",", "")),
                            ref=".",
                            alt=".",
                        )
                        try:
                            gene_context = GeneDescriptor(
                                value_id=self.gene_id_symbol[l[8]], symbol=l[8]
                            )
                        except KeyError:
                            gene_context = GeneDescriptor(symbol=l[8])
                        allelic_state = OntologyClass(
                            id=list(go_oi.basic_search(l[10].lower()))[0],
                            label=l[10].lower(),
                        )
                        variation_descriptor = VariationDescriptor(
                            id=self.job_file["analysis"]["proband"]
                            + ":"
                            + l[3]
                            + ":"
                            + l[4],
                            gene_context=gene_context,
                            vcf_record=vcf_record,
                            allelic_state=allelic_state,
                        )

                        variant_interpretation = VariantInterpretation(
                            acmg_pathogenicity_classification="NOT_PROVIDED",
                            therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                            variation_descriptor=variation_descriptor,
                        )

                        genomic_interpretation = GenomicInterpretation(
                            subject_or_biosample_id=self.job_file["analysis"][
                                "proband"
                            ],
                            interpretation_status=interpretation_status,
                            variant_interpretation=variant_interpretation,
                        )
                        genomic_interpretations.append(genomic_interpretation)
                        diagnosis = Diagnosis(
                            genomic_interpretations=genomic_interpretations
                        )
                        interpretation = Interpretation(
                            id=int_id, progress_status="SOLVED", diagnosis=diagnosis
                        )
            interpretations.append(interpretation)
        di.close()
        return interpretations

    def create_phenopacket(self):
        subject, ped, family_id = self.construct_subject_ped()
        phenopacket = Phenopacket(
            id=self.job_file["analysis"]["proband"],
            subject=subject,
            phenotypic_features=self.create_phenotypic_interpretations(),
            interpretations=self.create_interpretations(),
        )
        return ped, phenopacket, family_id

    def create_meta_data(self):
        files = [
            File(
                uri=self.job_file["analysis"]["vcf"],
                file_attributes={
                    "fileFormat": "VCF",
                    "genomeAssembly": self.job_file["analysis"]["genomeAssembly"],
                },
            )
        ]
        resources = [
            Resource(
                id="hp",
                name="human phenotype ontology",
                url="http://purl.obolibrary.org/obo/hp.owl",
                version="hp/releases/2019-11-08",
                namespace_prefix="HP",
                iri_prefix="http://purl.obolibrary.org/obo/HP_",
            )
        ]
        timestamp = Timestamp()
        timestamp.GetCurrentTime()
        meta_data = MetaData(
            created=timestamp,
            created_by="pheval-converter",
            resources=resources,
            phenopacket_schema_version="2.0",
        )
        return files, meta_data

    def create_family(self):
        ped, phenopacket, family_id = self.create_phenopacket()
        files, metadata = self.create_meta_data()
        family = Family(
            id=family_id,
            proband=phenopacket,
            pedigree=ped,
            files=files,
            meta_data=metadata,
        )
        return family

    def write_file(self):
        family = self.create_family()
        json = MessageToJson(family)
        with open(self.output_file, "w") as jsfile:
            jsfile.write(json)
        jsfile.close()


@click.command()
@click.option("--file-list", "-f", required=True, help="File to convert")
@click.option("--diagnoses", "-d", required=True, help="Diagnoses file")
def convert_yaml_to_family_phenopacket(file_list, diagnoses):
    files = DirectoryFiles(file_list, ".yml").obtain_files_suffix()
    for file in files:
        YamlToFamilyPhenopacketConversion(file, diagnoses).write_file()
