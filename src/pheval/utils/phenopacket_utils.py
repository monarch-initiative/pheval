import json

# import logging
import os
from collections import defaultdict
from copy import copy
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from google.protobuf.json_format import MessageToJson, Parse
from phenopackets import (
    Family,
    File,
    GenomicInterpretation,
    Interpretation,
    Phenopacket,
    PhenotypicFeature,
)

from pheval.prepare.custom_exceptions import IncorrectFileFormatError


class IncompatibleGenomeAssemblyError(Exception):
    """Exception raised for incompatible genome assembly."""

    def __init__(self, assembly, phenopacket, message="Incompatible Genome Assembly"):
        self.assembly: str = assembly
        self.phenopacket: Path = phenopacket
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} -> {self.assembly} in {self.phenopacket}"


@dataclass
class GenomicVariant:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass
class ProbandCausativeVariant:
    proband_id: str
    assembly: str
    variant: GenomicVariant
    genotype: str
    info: str = None


@dataclass
class ProbandCausativeGene:
    gene_symbol: str
    gene_identifier: str


def read_hgnc_data() -> pd.DataFrame:
    return pd.read_csv(
        os.path.dirname(__file__).replace("utils", "resources/hgnc_complete_set_2022-10-01.txt"),
        delimiter="\t",
        dtype=str,
    )


def create_hgnc_dict() -> defaultdict:
    """Creates reference for updating gene symbols and identifiers."""
    hgnc_df = read_hgnc_data()
    hgnc_data = defaultdict(dict)
    for _index, row in hgnc_df.iterrows():
        previous_names = []
        hgnc_data[row["symbol"]]["ensembl_id"] = row["ensembl_gene_id"]
        hgnc_data[row["symbol"]]["hgnc_id"] = row["hgnc_id"]
        hgnc_data[row["symbol"]]["entrez_id"] = row["entrez_id"]
        hgnc_data[row["symbol"]]["refseq_accession"] = row["refseq_accession"]
        previous = str(row["prev_symbol"]).split("|")
        for p in previous:
            previous_names.append(p.strip('"'))
        hgnc_data[row["symbol"]]["previous_symbol"] = previous_names

    return hgnc_data


def create_gene_identifier_map() -> dict:
    hgnc_df = read_hgnc_data()
    identifier_map = {}
    for _index, row in hgnc_df.iterrows():
        identifier_map[row["ensembl_gene_id"]] = row["symbol"]
        identifier_map[row["hgnc_id"]] = row["symbol"]
        identifier_map[row["entrez_id"]] = row["symbol"]
        identifier_map[row["refseq_accession"]] = row["symbol"]
    return identifier_map


def phenopacket_reader(file: Path):
    """Reads a phenopacket file, returning its contents."""
    file = open(file, "r")
    phenopacket = json.load(file)
    file.close()
    if "proband" in phenopacket:
        return Parse(json.dumps(phenopacket), Family())
    else:
        return Parse(json.dumps(phenopacket), Phenopacket())


class PhenopacketUtil:
    """Retrieves relevant data from a phenopacket."""

    def __init__(self, phenopacket_contents: Phenopacket):
        self.phenopacket_contents = phenopacket_contents

    def sample_id(self) -> str:
        """Retrieve the sample ID from a phenopacket or proband of a family."""
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.subject.id
        else:
            return self.phenopacket_contents.subject.id

    def phenotypic_features(self) -> list[PhenotypicFeature]:
        """Retrieves a list of all HPO terms."""
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.phenotypic_features
        else:
            return self.phenopacket_contents.phenotypic_features

    def observed_phenotypic_features(self) -> list[PhenotypicFeature]:
        """Removes any HPO terms labelled as excluded."""
        phenotypic_features = []
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                continue
            phenotypic_features.append(p)
        return phenotypic_features

    def negated_phenotypic_features(self) -> [PhenotypicFeature]:
        """Retrieve negated phenotypic features."""
        negated_phenotypic_features = []
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                negated_phenotypic_features.append(p)
        return negated_phenotypic_features

    def interpretations(self) -> list[Interpretation]:
        """Returns all interpretations of a phenopacket."""
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.interpretations
        else:
            return self.phenopacket_contents.interpretations

    def causative_variants(self) -> list[ProbandCausativeVariant]:
        """Returns a list of all causative variants listed in a phenopacket."""
        all_variants = []
        interpretation = self.interpretations()
        for i in interpretation:
            for g in i.diagnosis.genomic_interpretations:
                vcf_record = g.variant_interpretation.variation_descriptor.vcf_record
                genotype = g.variant_interpretation.variation_descriptor.allelic_state
                variant_data = ProbandCausativeVariant(
                    self.phenopacket_contents.subject.id,
                    vcf_record.genome_assembly,
                    GenomicVariant(
                        vcf_record.chrom,
                        vcf_record.pos,
                        vcf_record.ref,
                        vcf_record.alt,
                    ),
                    genotype.label,
                    vcf_record.info,
                )
                all_variants.append(variant_data)
        return all_variants

    def files(self) -> list:
        """Returns all files associated with a phenopacket."""
        return self.phenopacket_contents.files

    def vcf_file_data(self, phenopacket_path: Path, vcf_dir: Path) -> File:
        """Retrieves the genome assembly and vcf name from a phenopacket."""
        compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
        vcf_data = [file for file in self.files() if file.file_attributes["fileFormat"] == "vcf"][0]
        if not Path(vcf_data.uri).name.endswith(".vcf") and not Path(vcf_data.uri).name.endswith(
            ".vcf.gz"
        ):
            raise IncorrectFileFormatError(Path(vcf_data.uri), ".vcf or .vcf.gz file")
        if vcf_data.file_attributes["genomeAssembly"] not in compatible_genome_assembly:
            raise IncompatibleGenomeAssemblyError(
                vcf_data.file_attributes["genomeAssembly"], phenopacket_path
            )
        vcf_data.uri = str(vcf_dir.joinpath(Path(vcf_data.uri).name))
        return vcf_data

    @staticmethod
    def _extract_diagnosed_gene(
        genomic_interpretation: GenomicInterpretation,
    ) -> ProbandCausativeGene:
        """Returns the disease causative gene from the variant descriptor field if not empty,
        otherwise, returns from the gene descriptor from a phenopacket."""
        if genomic_interpretation.variant_interpretation.ByteSize() != 0:
            return ProbandCausativeGene(
                genomic_interpretation.variant_interpretation.variation_descriptor.gene_context.symbol,
                genomic_interpretation.variant_interpretation.variation_descriptor.gene_context.value_id,
            )

        else:
            return ProbandCausativeGene(
                gene_symbol=genomic_interpretation.gene.symbol,
                gene_identifier=genomic_interpretation.gene.value_id,
            )

    def diagnosed_genes(self) -> list[ProbandCausativeGene]:
        """Returns a unique list of all causative genes and the corresponding gene identifiers from a phenopacket."""
        pheno_interpretation = self.interpretations()
        genes = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                genes.append(self._extract_diagnosed_gene(g))
                genes = list({gene.gene_symbol: gene for gene in genes}.values())
        return genes

    def diagnosed_variants(self) -> list[GenomicVariant]:
        """Returns a list of all variants from a phenopacket - for use in assess-prioritisation."""
        variants = []
        pheno_interpretation = self.interpretations()
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                variant = GenomicVariant(
                    chrom=g.variant_interpretation.variation_descriptor.vcf_record.chrom,
                    pos=g.variant_interpretation.variation_descriptor.vcf_record.pos,
                    ref=g.variant_interpretation.variation_descriptor.vcf_record.ref,
                    alt=g.variant_interpretation.variation_descriptor.vcf_record.alt,
                )
                variants.append(variant)
        return variants


class PhenopacketRebuilder:
    """Rebuilds a Phenopacket."""

    def __init__(self, phenopacket: Phenopacket or Family):
        self.phenopacket = phenopacket

    def update_interpretations(self, interpretations) -> Phenopacket or Family:
        """Adds the updated interpretations to a phenopacket."""
        phenopacket = copy(self.phenopacket)
        if hasattr(phenopacket, "proband"):
            del phenopacket.proband.interpretations[:]
            phenopacket.proband.interpretations.extend(interpretations)
        else:
            del phenopacket.interpretations[:]
            phenopacket.interpretations.extend(interpretations)
        return phenopacket

    def add_randomised_hpo(self, randomised_hpo) -> Phenopacket or Family:
        """Adds randomised phenotypic profile to phenopacket."""
        phenopacket = copy(self.phenopacket)
        if hasattr(phenopacket, "proband"):
            del phenopacket.proband.phenotypic_features[:]
            phenopacket.proband.phenotypic_features.extend(randomised_hpo)
        else:
            del phenopacket.phenotypic_features[:]
            phenopacket.phenotypic_features.extend(randomised_hpo)
        return phenopacket

    def add_spiked_vcf_path(self, spiked_vcf_file_data: File) -> Phenopacket or Family:
        """Adds spiked vcf path to phenopacket."""
        phenopacket = copy(self.phenopacket)
        phenopacket_files = [
            file for file in phenopacket.files if file.file_attributes["fileFormat"] != "vcf"
        ]
        phenopacket_files.append(spiked_vcf_file_data)
        del phenopacket.files[:]
        phenopacket.files.extend(phenopacket_files)
        return phenopacket


def create_json_message(phenopacket: Phenopacket or Family) -> str:
    """Creates json message for writing to file."""
    return MessageToJson(phenopacket)


def write_phenopacket(phenopacket: Phenopacket or Family, output_file: Path) -> None:
    """Writes a phenopacket."""
    phenopacket_json = create_json_message(phenopacket)
    with open(output_file, "w") as outfile:
        outfile.write(phenopacket_json)
    outfile.close()


class GeneIdentifierUpdater:
    def __init__(self, gene_identifier: str, hgnc_data: dict = None, identifier_map: dict = None):
        self.hgnc_data = hgnc_data
        self.gene_identifier = gene_identifier
        self.identifier_map = identifier_map

    def find_identifier(self, gene_symbol: str) -> str:
        """Finds the specified gene identifier for a gene symbol."""
        if gene_symbol in self.hgnc_data.keys():
            return self.hgnc_data[gene_symbol][self.gene_identifier]
        else:
            for _symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return data[self.gene_identifier]

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier. (e.g.)
        "
        obtain_gene_symbol_from_identifier(query_gene_identifier="HGNC:5")
        "
        """
        return self.identifier_map[query_gene_identifier]

    def _find_alternate_ids(self, gene_symbol: str) -> list[str]:
        """Finds the alternate IDs for a gene symbol."""
        if gene_symbol in self.hgnc_data.keys():
            return [
                self.hgnc_data[gene_symbol]["hgnc_id"],
                "ncbigene:" + self.hgnc_data[gene_symbol]["entrez_id"],
                "ensembl:" + self.hgnc_data[gene_symbol]["ensembl_id"],
                "symbol:" + gene_symbol,
            ]
        else:
            for symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return [
                            data["hgnc_id"],
                            "ncbigene:" + data["entrez_id"],
                            "ensembl:" + data["ensembl_id"],
                            "symbol:" + symbol,
                        ]

    def update_genomic_interpretations_gene_identifier(
        self, interpretations: list[Interpretation]
    ) -> list[Interpretation]:
        """Updates the genomic interpretations of a phenopacket."""
        updated_interpretations = copy(list(interpretations))
        for updated_interpretation in updated_interpretations:
            for g in updated_interpretation.diagnosis.genomic_interpretations:
                g.variant_interpretation.variation_descriptor.gene_context.value_id = (
                    self.find_identifier(
                        g.variant_interpretation.variation_descriptor.gene_context.symbol
                    )
                )
                del g.variant_interpretation.variation_descriptor.gene_context.alternate_ids[:]
                g.variant_interpretation.variation_descriptor.gene_context.alternate_ids.extend(
                    self._find_alternate_ids(
                        g.variant_interpretation.variation_descriptor.gene_context.symbol
                    )
                )
        return updated_interpretations
