import json

# import logging
import os
from collections import defaultdict
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

import pandas as pd
from google.protobuf.json_format import MessageToJson, Parse
from phenopackets import (
    Disease,
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
        """
        Initialise IncompatibleGenomeAssemblyError.

        Attributes:
           assembly (str): Incompatible genome assembly encountered.
           phenopacket (Path): Path to the Phenopacket associated with the error.
           message (str, optional): Custom error message (default is "Incompatible Genome Assembly").
        """
        self.assembly: str = assembly
        self.phenopacket: Path = phenopacket
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} -> {self.assembly} in {self.phenopacket}"


@dataclass
class GenomicVariant:
    """
    Represents a genomic variant.

    Args:
        chrom (str): The chromosome position of the variant recommended to be provided in the following format.
        This includes numerical designations from 1 to 22 representing autosomal chromosomes,
        as well as the sex chromosomes X and Y, and the mitochondrial chromosome MT.
        pos (int): Position of the variant following VCF convention.
        ref (str): Reference allele following VCF convention.
        alt (str): Alternate allele following VCF convention.
    """

    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass
class ProbandCausativeVariant:
    """
    Represents a causative variant associated with a proband

    Args:
        proband_id (str): ID of the proband
        assembly (str): Genome assembly
        variant (GenomicVariant): Genomic variant associated with the proband
        genotype (str): Genotype information for the variant
        info (str, optional): Additional information about the variant (default is an empty string)
    """

    proband_id: str
    assembly: str
    variant: GenomicVariant
    genotype: str
    info: str = ""


@dataclass
class ProbandCausativeGene:
    """
    Represents a causative gene associated with a proband

    Args:
        gene_symbol (str): Symbol representing the gene
        gene_identifier (str): The ENSEMBL gene identifier for the result entry
    Notes:
         While we recommend providing the gene identifier in the ENSEMBL namespace,
         any matching format used in Phenopacket interpretations and result output is acceptable
         for result matching purposes in the analysis.
    """

    gene_symbol: str
    gene_identifier: str


@dataclass(frozen=True, eq=True)
class ProbandDisease:
    """
    Represents a disease associated with a proband

    Args:
        disease_name (str): Name of the disease
        disease_identifier (str): Identifier for the disease result entry in the OMIM namespace

    Notes:
         While we recommend providing the disease identifier in the OMIM namespace,
         any matching format used in Phenopacket interpretations and result output is acceptable
         for result matching purposes in the analysis.
    """

    disease_name: str
    disease_identifier: str


def read_hgnc_data() -> pd.DataFrame:
    """
    Read HGNC data from a file and return it as a Pandas DataFrame.

    Returns:
        pd.DataFrame: DataFrame containing the HGNC data.
    """
    return pd.read_csv(
        os.path.dirname(__file__).replace("utils", "resources/hgnc_complete_set.txt"),
        delimiter="\t",
        dtype=str,
    )


def create_hgnc_dict() -> defaultdict:
    """
    Create a dictionary as a reference for updating gene symbols and identifiers based on HGNC data.


    Returns:
        defaultdict: A dictionary containing gene symbols as keys and their associated gene information.

    Notes:
        The dictionary structure:
        {
            'gene_symbol': {
                'ensembl_id': str,
                'hgnc_id': str,
                'entrez_id': str,
                'refseq_accession': str,
                'previous_symbol': [str, ...]
            },
            ...
        }
    """
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
    """
    Create a mapping of gene identifiers to gene symbols using HGNC data.

    Returns:
        dict: A mapping of gene identifiers to gene symbols.

    Notes:
        The dictionary structure:
        {
            'identifier': 'gene_symbol',
            ...
        }
    """
    hgnc_df = read_hgnc_data()
    identifier_map = {}
    for _index, row in hgnc_df.iterrows():
        identifier_map[row["ensembl_gene_id"]] = row["symbol"]
        identifier_map[row["hgnc_id"]] = row["symbol"]
        identifier_map[row["entrez_id"]] = row["symbol"]
        identifier_map[row["refseq_accession"]] = row["symbol"]
    return identifier_map


def phenopacket_reader(file: Path) -> Union[Phenopacket, Family]:
    """
    Read a Phenopacket file and returns its contents as a Phenopacket or Family object

    Args:
        file (Path): Path to the Phenopacket file

    Returns:
        Union[Phenopacket, Family]: Contents of the Phenopacket file as a Phenopacket or Family object
    """
    file = open(file, "r")
    phenopacket = json.load(file)
    file.close()
    if "proband" in phenopacket:
        return Parse(json.dumps(phenopacket), Family())
    else:
        return Parse(json.dumps(phenopacket), Phenopacket())


class PhenopacketUtil:
    """Class for retrieving data from a Phenopacket or Family object"""

    def __init__(self, phenopacket_contents: Union[Phenopacket, Family]):
        """Initialise PhenopacketUtil

        Args:
            phenopacket_contents (Union[Phenopacket, Family]): Phenopacket or Family object
        """
        self.phenopacket_contents = phenopacket_contents

    def sample_id(self) -> str:
        """
        Retrieve the sample ID from a Phenopacket or proband of a Family

        Returns:
            str: Sample ID
        """
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.subject.id
        else:
            return self.phenopacket_contents.subject.id

    def phenotypic_features(self) -> List[PhenotypicFeature]:
        """
        Retrieve a list of all HPO terms

        Returns:
            List[PhenotypicFeature]: List of HPO terms
        """
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.phenotypic_features
        else:
            return self.phenopacket_contents.phenotypic_features

    def observed_phenotypic_features(self) -> List[PhenotypicFeature]:
        """
        Retrieve a list of all observed HPO terms

        Returns:
            List[PhenotypicFeature]: List of observed HPO terms
        """
        phenotypic_features = []
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                continue
            phenotypic_features.append(p)
        return phenotypic_features

    def negated_phenotypic_features(self) -> List[PhenotypicFeature]:
        """
        Retrieve a list of all negated HPO terms

        Returns:
            List[PhenotypicFeature]: List of negated HPO terms
        """
        negated_phenotypic_features = []
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                negated_phenotypic_features.append(p)
        return negated_phenotypic_features

    def diseases(self) -> List[Disease]:
        """
        Retrieve a list of Diseases associated with the proband

        Returns:
            List[Disease]: List of diseases
        """
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.diseases
        else:
            return self.phenopacket_contents.diseases

    def _diagnosis_from_interpretations(self) -> List[ProbandDisease]:
        """
        Retrieve a list of disease diagnoses associated with the proband from the interpretations object

        Returns:
            List[ProbandDisease]: List of diagnosed diseases
        """
        diagnoses = []
        interpretation = self.interpretations()
        for i in interpretation:
            (
                diagnoses.append(
                    ProbandDisease(
                        disease_name=i.diagnosis.disease.label,
                        disease_identifier=i.diagnosis.disease.id,
                    )
                )
                if i.diagnosis.disease.label != "" and i.diagnosis.disease.id != ""
                else None
            )
        return diagnoses

    def _diagnosis_from_disease(self) -> List[ProbandDisease]:
        """
        Retrieve a list of disease diagnoses associated with the proband from the diseases object

        Returns:
            List[ProbandDisease]: List of diagnosed diseases
        """
        diagnoses = []
        for disease in self.diseases():
            diagnoses.append(
                ProbandDisease(disease_name=disease.term.label, disease_identifier=disease.term.id)
            )
        return diagnoses

    def diagnoses(self) -> List[ProbandDisease]:
        """
        Retrieve a unique list of disease diagnoses associated with the proband from a Phenopacket

        Returns:
            List[ProbandDisease]: List of diagnosed diseases
        """
        return list(set(self._diagnosis_from_interpretations() + self._diagnosis_from_disease()))

    def interpretations(self) -> List[Interpretation]:
        """
        Retrieve a list of interpretations from a Phenopacket

        Returns:
            List[Interpretation]: List of interpretations
        """
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.interpretations
        else:
            return self.phenopacket_contents.interpretations

    def causative_variants(self) -> List[ProbandCausativeVariant]:
        """
        Retrieve a list of causative variants listed in a Phenopacket

        Returns:
            List[ProbandCausativeVariant]: List of proband causative variants
        """
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

    def files(self) -> List[File]:
        """
        Retrieve a list of files associated with a phenopacket

        Returns:
            List[File]: List of files associated with a phenopacket
        """
        return self.phenopacket_contents.files

    def vcf_file_data(self, phenopacket_path: Path, vcf_dir: Path) -> File:
        """
        Retrieve the genome assembly and VCF file name from a phenopacket.

        Args:
            phenopacket_path (Path): The path to the phenopacket file.
            vcf_dir (Path): The directory path where the VCF file is stored.

        Returns:
            File: The VCF file with updated URI pointing to the specified directory.

        Raises:
            IncorrectFileFormatError: If the provided file is not in .vcf or .vcf.gz format.
            IncompatibleGenomeAssemblyError: If the genome assembly of the VCF file is not compatible.

        Note:
            This function searches for a VCF file within the provided list of files, validates its format,
            and checks if the genome assembly is compatible. If the conditions are met, it updates the
            URI of the VCF file to the specified directory and returns the modified file object.
        """
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
        """
        Retrieve the disease causing genes from the variant descriptor field if not empty,
        otherwise, retrieves from the gene descriptor from a phenopacket.
        Args:
            genomic_interpretation (GenomicInterpretation): A genomic interpretation from a Phenopacket
        Returns:
            ProbandCausativeGene: The disease causing gene
        """
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

    def diagnosed_genes(self) -> List[ProbandCausativeGene]:
        """
        Retrieve the disease causing genes from a phenopacket.
        Returns:
            List[ProbandCausativeGene]: List of causative genes
        """
        pheno_interpretation = self.interpretations()
        genes = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                genes.append(self._extract_diagnosed_gene(g))
                genes = list({gene.gene_symbol: gene for gene in genes}.values())
        return genes

    def diagnosed_variants(self) -> List[GenomicVariant]:
        """
        Retrieve a list of all known causative variants from a phenopacket.
        Returns:
            List[GenomicVariant]: List of causative variants
        """
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
    """Class for rebuilding a Phenopacket"""

    def __init__(self, phenopacket: Union[Phenopacket, Family]):
        """Initialise PhenopacketUtil

        Attributes:
            phenopacket (Union[Phenopacket, Family]): Phenopacket or Family object
        """
        self.phenopacket = phenopacket

    def update_interpretations(
        self, interpretations: [Interpretation]
    ) -> Union[Phenopacket, Family]:
        """
        Add the updated interpretations to a Phenopacket or Family.

        Args:
            interpretations (List[Interpretation]): The updated interpretations to be added.

        Returns:
            Union[Phenopacket, Family]: The Phenopacket or Family object with updated interpretations.
        """
        phenopacket = copy(self.phenopacket)
        if hasattr(phenopacket, "proband"):
            del phenopacket.proband.interpretations[:]
            phenopacket.proband.interpretations.extend(interpretations)
        else:
            del phenopacket.interpretations[:]
            phenopacket.interpretations.extend(interpretations)
        return phenopacket

    def add_randomised_hpo(self, randomised_hpo: [PhenotypicFeature]) -> Union[Phenopacket, Family]:
        """
        Add randomised phenotypic profiles to a Phenopacket or Family.

        Args:
            randomised_hpo: The randomised phenotypic profiles to be added.

        Returns:
            Union[Phenopacket, Family] The Phenopacket or Family object with added randomised profiles.
        """
        phenopacket = copy(self.phenopacket)
        if hasattr(phenopacket, "proband"):
            del phenopacket.proband.phenotypic_features[:]
            phenopacket.proband.phenotypic_features.extend(randomised_hpo)
        else:
            del phenopacket.phenotypic_features[:]
            phenopacket.phenotypic_features.extend(randomised_hpo)
        return phenopacket

    def add_spiked_vcf_path(self, spiked_vcf_file_data: File) -> Union[Phenopacket, Family]:
        """
        Add a spiked VCF path to a Phenopacket or Family.

        Args:
        - spiked_vcf_file_data (File): The VCF file data to be added.

        Returns:
        - Phenopacket or Family: The Phenopacket or Family object with the added spiked VCF path.
        """
        phenopacket = copy(self.phenopacket)
        phenopacket_files = [
            file for file in phenopacket.files if file.file_attributes["fileFormat"] != "vcf"
        ]
        phenopacket_files.append(spiked_vcf_file_data)
        del phenopacket.files[:]
        phenopacket.files.extend(phenopacket_files)
        return phenopacket


def create_json_message(phenopacket: Union[Phenopacket, Family]) -> str:
    """
    Create a JSON message for writing to a file.

    Args:
    - phenopacket (Union[Phenopacket, Family]): The Phenopacket or Family object to convert to JSON.

    Returns:
    - str: A JSON-formatted string representation of the Phenopacket or Family object.
    """
    return MessageToJson(phenopacket)


def write_phenopacket(phenopacket: Union[Phenopacket, Family], output_file: Path) -> None:
    """
    Write a Phenopacket or Family object to a file in JSON format.

    Args:
        phenopacket (Phenopacket or Family): The Phenopacket or Family object to be written.
        output_file (Path): The Path object representing the file to write the Phenopacket data.

    Returns:
        None
    """
    phenopacket_json = create_json_message(phenopacket)
    with open(output_file, "w") as outfile:
        outfile.write(phenopacket_json)
    outfile.close()


class GeneIdentifierUpdater:
    """Class for updating gene identifiers within genomic interpretations."""

    def __init__(self, gene_identifier: str, hgnc_data: dict = None, identifier_map: dict = None):
        """
        Initialise the GeneIdentifierUpdater.

        Args:
            gene_identifier (str): The gene identifier to update to.
            hgnc_data (dict): A dictionary containing HGNC data (default: None).
            identifier_map (dict): A dictionary mapping gene identifiers (default: None).
        """

        self.hgnc_data = hgnc_data
        self.gene_identifier = gene_identifier
        self.identifier_map = identifier_map

    def find_identifier(self, gene_symbol: str) -> str:
        """
        Find the specified gene identifier for a gene symbol.

        Args:
            gene_symbol (str): The gene symbol to find the identifier for.

        Returns:
            str: The identified gene identifier.
        """
        if gene_symbol in self.hgnc_data.keys():
            return self.hgnc_data[gene_symbol][self.gene_identifier]
        else:
            for _symbol, data in self.hgnc_data.items():
                for prev_symbol in data["previous_symbol"]:
                    if prev_symbol == gene_symbol:
                        return data[self.gene_identifier]

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier.

        Args:
            query_gene_identifier (str): The gene identifier.

        Returns:
            str: The gene symbol corresponding to the identifier.
        """
        return self.identifier_map[query_gene_identifier]

    def _find_alternate_ids(self, gene_symbol: str) -> List[str]:
        """
        Find the alternate IDs for a gene symbol.

        Args:
            gene_symbol (str): The gene symbol to find alternate IDs for.

        Returns:
            List[str]: List of alternate IDs for the gene symbol.
        """
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
        self, interpretations: List[Interpretation]
    ) -> List[Interpretation]:
        """
        Update the genomic interpretations of a Phenopacket.

        Args:
            interpretations (List[Interpretation]): List of Interpretation objects.

        Returns:
            List[Interpretation]: Updated list of Interpretation objects.
        """
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
