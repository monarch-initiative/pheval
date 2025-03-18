import json
import logging
import os
from copy import copy
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

import polars as pl
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

info_log = logging.getLogger("info")


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


def parse_hgnc_data() -> pl.DataFrame:
    """
    Read HGNC data from a file and return it as a Polars DataFrame.

    Returns:
        pl.DataFrame: DataFrame containing the HGNC data.
    """
    return (
        pl.read_csv(
            os.path.dirname(__file__).replace("utils", "resources/hgnc_complete_set.txt"),
            separator="\t",
            infer_schema=10000000000,
            dtypes={"omim_id": pl.Utf8},
        )
        .select(
            [
                pl.col("hgnc_id").alias("hgnc_id"),
                pl.col("symbol").alias("gene_symbol"),
                pl.col("ensembl_gene_id").alias("ensembl_id"),
                pl.col("entrez_id").alias("entrez_id"),
                pl.col("refseq_accession").alias("refseq_accession"),
                pl.col("prev_symbol").alias("previous_symbol_raw"),
            ]
        )
        .with_columns(
            pl.col("previous_symbol_raw")
            .str.split("|")
            .list.eval(pl.element().str.strip_chars('"'))
            .alias("prev_symbols")
        )
    )


def create_gene_identifier_map() -> pl.DataFrame:
    """
    Create a mapping of gene identifiers to gene symbols using HGNC data.

    Returns:
        pl.DataFrame: A mapping of gene identifiers to gene symbols.
    """
    hgnc_df = parse_hgnc_data()
    return hgnc_df.melt(
        id_vars=["gene_symbol", "prev_symbols"],
        value_vars=["ensembl_id", "hgnc_id", "entrez_id", "refseq_accession"],
        variable_name="identifier_type",
        value_name="identifier",
    ).with_columns(
        pl.col("identifier_type")
        .replace(
            {
                "ensembl_id": "ensembl:",
                "hgnc_id": "",
                "entrez_id": "ncbigene:",
                "refseq_accession": "",
            },
            default="",
        )
        .alias("prefix")
    )


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
                    chrom=str(
                        g.variant_interpretation.variation_descriptor.vcf_record.chrom.replace(
                            "chr", ""
                        )
                    ),
                    pos=int(g.variant_interpretation.variation_descriptor.vcf_record.pos),
                    ref=g.variant_interpretation.variation_descriptor.vcf_record.ref,
                    alt=g.variant_interpretation.variation_descriptor.vcf_record.alt,
                )
                variants.append(variant)
        return variants

    def check_incomplete_variant_record(self) -> bool:
        """
        Check if any variant record in the phenopacket has incomplete information.

        This method iterates through the diagnosed variant records and checks if any of them
        have missing or incomplete information such as empty chromosome, position, reference,
        or alternate allele.

        Returns:
            bool: True if any variant record is incomplete, False otherwise.
        """
        variants = self.diagnosed_variants()
        for variant in variants:
            if (
                variant.chrom == ""
                or variant.pos == 0
                or variant.pos == ""
                or variant.ref == ""
                or variant.alt == ""
            ):
                return True
        return False

    def check_variant_alleles(self) -> bool:
        """
        Check if any variant record in the phenopacket has identical reference and alternate alleles.

        Returns:
            bool: True if the reference and alternate alleles are identical, False otherwise.
        """
        variants = self.diagnosed_variants()
        for variant in variants:
            if variant.ref == variant.alt:
                return True
        return False

    def check_incomplete_gene_record(self) -> bool:
        """
        Check if any gene record in the phenopacket has incomplete information.

        This method iterates through the diagnosed gene records and checks if any of them
        have missing or incomplete information such as gene name, or gene identifier.

        Returns:
            bool: True if any gene record is incomplete, False otherwise.
        """
        genes = self.diagnosed_genes()
        for gene in genes:
            if gene.gene_symbol == "" or gene.gene_identifier == "":
                return True
        return False

    def check_incomplete_disease_record(self) -> bool:
        """
        Check if any disease record in the phenopacket has incomplete information.

        This method iterates through the diagnosed disease records and checks if any of them
        have missing or incomplete information such as empty disease name, or disease identifier.

        Returns:
            bool: True if any disease record is incomplete, False otherwise.
        """
        if len(self.diagnoses()) == 0:
            return True
        return False


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

    def __init__(
        self,
        gene_identifier: str,
        identifier_map: pl.DataFrame = None,
    ):
        """
        Initialise the GeneIdentifierUpdater.

        Args:
            gene_identifier (str): The gene identifier to update to.
            identifier_map (dict): A polars dataframe mapping gene identifiers (default: None).
        """

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
        matches = self.identifier_map.filter(
            (pl.col("gene_symbol") == gene_symbol)
            & (pl.col("identifier_type") == self.gene_identifier)
        )

        if matches.height > 0:
            return matches["identifier"][0]
        prev_symbol_matches = self.identifier_map.filter(
            (pl.col("identifier_type") == self.gene_identifier)
            & (pl.col("prev_symbols").list.contains(gene_symbol))
        )
        if prev_symbol_matches.height > 0:
            return prev_symbol_matches["identifier"][0]
        return None

    def obtain_gene_symbol_from_identifier(self, query_gene_identifier: str) -> str:
        """
        Obtain gene symbol from a gene identifier.

        Args:
            query_gene_identifier (str): The gene identifier.

        Returns:
            str: The gene symbol corresponding to the identifier.
        """
        return self.identifier_map.filter(pl.col("identifier") == query_gene_identifier)[
            "gene_symbol"
        ][0]

    def _find_alternate_ids(self, gene_symbol: str) -> List[str]:
        """
        Find the alternate IDs for a gene symbol.

        Args:
            gene_symbol (str): The gene symbol to find alternate IDs for.

        Returns:
            List[str]: List of alternate IDs for the gene symbol.
        """
        matches = self.identifier_map.filter((pl.col("gene_symbol") == gene_symbol))
        if matches.height > 0:
            return [f"{row['prefix']}{row['identifier']}" for row in matches.rows(named=True)] + [
                f"symbol:{gene_symbol}"
            ]
        prev_symbol_matches = self.identifier_map.filter(
            (pl.col("prev_symbols").list.contains(gene_symbol))
        )
        if prev_symbol_matches.height > 0:
            return [
                f"{row['prefix']}{row['identifier']}"
                for row in prev_symbol_matches.rows(named=True)
            ] + [f"symbol:{gene_symbol}"]
        return None

    def update_genomic_interpretations_gene_identifier(
        self, interpretations: List[Interpretation], phenopacket_path: Path
    ) -> List[Interpretation]:
        """
        Update the genomic interpretations of a Phenopacket.

        Args:
            interpretations (List[Interpretation]): List of Interpretation objects.
            phenopacket_path (Path): The Path to the Phenopacket.

        Returns:
            List[Interpretation]: Updated list of Interpretation objects.
        """
        updated_interpretations = copy(list(interpretations))
        for updated_interpretation in updated_interpretations:
            for g in updated_interpretation.diagnosis.genomic_interpretations:
                updated_gene_identifier = self.find_identifier(
                    g.variant_interpretation.variation_descriptor.gene_context.symbol
                )
                info_log.info(
                    f"Updating gene identifier in {phenopacket_path} from "
                    f"{g.variant_interpretation.variation_descriptor.gene_context.value_id}"
                    f"to {updated_gene_identifier}"
                )
                g.variant_interpretation.variation_descriptor.gene_context.value_id = (
                    updated_gene_identifier
                )
                del g.variant_interpretation.variation_descriptor.gene_context.alternate_ids[:]
                g.variant_interpretation.variation_descriptor.gene_context.alternate_ids.extend(
                    self._find_alternate_ids(
                        g.variant_interpretation.variation_descriptor.gene_context.symbol
                    )
                )
        return updated_interpretations
