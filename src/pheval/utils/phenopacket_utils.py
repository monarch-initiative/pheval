import json
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from google.protobuf.json_format import MessageToJson, Parse
from phenopackets import Family, Phenopacket

from pheval.prepare.custom_exceptions import IncorrectFileFormatError


class IncompatibleGenomeAssemblyError(Exception):
    """Exception raised for incompatible genome assembly."""

    def __init__(self, assembly, phenopacket, message="Incompatible Genome Assembly"):
        self.assembly: str = assembly
        self.phenopacket: str = phenopacket
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} -> {self.assembly} in {self.phenopacket}"


@dataclass
class CausativeVariant:
    phenopacket: Path
    proband_id: str
    assembly: str
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str


class PhenopacketReader:
    """Class for reading Phenopackets and retrieving relevant data."""

    def __init__(self, file: Path):
        self.file_name = file
        file = open(file, "r")
        phenopacket = json.load(file)
        if "proband" in phenopacket:
            self.pheno = Parse(json.dumps(phenopacket), Family())
        else:
            self.pheno = Parse(json.dumps(phenopacket), Phenopacket())
        file.close()

    def phenotypic_features(self):
        """Retrieves a list of all HPO terms."""
        if hasattr(self.pheno, "proband"):
            all_phenotypic_features = self.pheno.proband.phenotypic_features
        else:
            all_phenotypic_features = self.pheno.phenotypic_features
        return all_phenotypic_features

    def remove_excluded_phenotypic_features(self) -> dict:
        """Removes any HPO terms labelled as excluded."""
        phenotypic_features = {}
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                continue
            phenotypic_features[p.type.id] = p.type.label
        return phenotypic_features

    def interpretations(self) -> list:
        """Returns all interpretations of a phenopacket."""
        if hasattr(self.pheno, "proband"):
            pheno_interpretation = self.pheno.proband.interpretations
        else:
            pheno_interpretation = self.pheno.interpretations
        return pheno_interpretation

    def causative_variants(self) -> list[CausativeVariant]:
        """Returns all causative variants listed in a phenopacket - for use in spiked VCF."""
        all_variants = []
        interpretation = self.interpretations()
        for i in interpretation:
            for g in i.diagnosis.genomic_interpretations:
                record = g.variant_interpretation.variation_descriptor.vcf_record
                genotype = g.variant_interpretation.variation_descriptor.allelic_state
                variant = CausativeVariant(
                    self.file_name,
                    self.pheno.subject.id,
                    record.genome_assembly,
                    record.chrom,
                    record.pos,
                    record.ref,
                    record.alt,
                    genotype.label,
                )
                all_variants.append(variant)
        return all_variants

    def files(self) -> list:
        """Returns all files associated with a phenopacket."""
        return self.pheno.files

    def vcf_file_data(self, vcf_dir: Path) -> tuple[str, str]:
        """Retrieves the genome assembly and vcf name from a phenopacket."""
        compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
        ppacket_files = self.files()
        vcf, genome_assembly = "", ""
        for file in ppacket_files:
            vcf_name = ""
            if file.file_attributes["fileFormat"].upper() == "VCF":
                vcf_path = file.uri
                if "/" in vcf_path:
                    vcf_name = vcf_path.rsplit("/")[-1]
                if "/" not in vcf_path:
                    vcf_name = vcf_path
            if not vcf_name.endswith(".vcf") and not vcf_name.endswith(".vcf.gz"):
                raise IncorrectFileFormatError(vcf_name, ".vcf or .vcf.gz file")
            vcf = os.path.join(vcf_dir, vcf_name)
            assembly = file.file_attributes["genomeAssembly"]
            if assembly not in compatible_genome_assembly:
                raise IncompatibleGenomeAssemblyError(assembly, self.file_name)
            genome_assembly = assembly
        return vcf, genome_assembly

    def diagnosed_genes(self) -> list:
        """Returns a unique list of all causative genes from a phenopacket."""
        pheno_interpretation = self.interpretations()
        genes = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                genes.append(g.variant_interpretation.variation_descriptor.gene_context.symbol)
        genes = list(set(genes))
        return genes

    def diagnosed_variants(self) -> list:
        """Returns a list of all causative variants from a phenopacket - for use in assess-prioritisation."""
        variants = []
        pheno_interpretation = self.interpretations()
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                variant = defaultdict(dict)
                variant[
                    "geneSymbol"
                ] = g.variant_interpretation.variation_descriptor.gene_context.symbol
                variant["variant"] = g.variant_interpretation.variation_descriptor.vcf_record
                variants.append(variant)
        return variants


class PhenopacketRebuilder:
    """Rebuilds Phenopacket structure for writing to a file."""

    def __init__(self, phenopacket_contents: PhenopacketReader):
        self.phenopacket_contents = phenopacket_contents

    def add_randomised_hpo(self, hpo_list):
        """Adds randomised phenotypic profile to phenopacket."""
        randomised_ppacket = Phenopacket(
            id=self.phenopacket_contents.pheno.id,
            subject=self.phenopacket_contents.pheno.subject,
            phenotypic_features=hpo_list,
            interpretations=self.phenopacket_contents.pheno.interpretations,
            files=self.phenopacket_contents.pheno.files,
            meta_data=self.phenopacket_contents.pheno.meta_data,
        )
        return randomised_ppacket

    def create_json_message(self, phenopacket):
        """Creates json message for writing to file."""
        if hasattr(self.phenopacket_contents.pheno, "proband"):
            family = Family(
                id=self.phenopacket_contents.pheno.id,
                proband=phenopacket,
                pedigree=self.phenopacket_contents.pheno.pedigree,
                files=self.phenopacket_contents.pheno.files,
                meta_data=self.phenopacket_contents.pheno.meta_data,
            )
            altered_phenopacket = MessageToJson(family)
        else:
            altered_phenopacket = MessageToJson(phenopacket)
        return altered_phenopacket


class PhenopacketWriter:
    """Writes phenopackets back to file."""

    def __init__(self, altered_phenopacket, output_file: Path):
        self.altered_phenopacket = altered_phenopacket
        self.output_file = output_file

    def write_file(self):
        """Writes a json message to a phenopacket .json file."""
        with open(self.output_file, "w") as outfile:
            outfile.write(self.altered_phenopacket)
        outfile.close()
