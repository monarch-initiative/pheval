import json
import os
import typing
from collections import defaultdict, namedtuple
from dataclasses import dataclass
from pathlib import Path

from google.protobuf.json_format import MessageToJson, Parse
from phenopackets import Family, File, Phenopacket, PhenotypicFeature

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


def phenopacket_reader(file: Path):
    """Reads a phenopacket file, returning a Phenopacket or Family object."""
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

    def phenotypic_features(self) -> list[PhenotypicFeature]:
        """Retrieves a list of all HPO terms."""
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.phenotypic_features
        else:
            return self.phenopacket_contents.phenotypic_features

    def remove_excluded_phenotypic_features(self) -> list[PhenotypicFeature]:
        """Removes any HPO terms labelled as excluded."""
        phenotypic_features = []
        all_phenotypic_features = self.phenotypic_features()
        for p in all_phenotypic_features:
            if p.excluded:
                continue
            phenotypic_features.append(p)
        return phenotypic_features

    def interpretations(self) -> list:
        """Returns all interpretations of a phenopacket."""
        if hasattr(self.phenopacket_contents, "proband"):
            return self.phenopacket_contents.proband.interpretations
        else:
            return self.phenopacket_contents.interpretations

    def causative_variants(self, file: Path) -> list[CausativeVariant]:
        """Returns a list of all causative variants listed in a phenopacket."""
        all_variants = []
        interpretation = self.interpretations()
        for i in interpretation:
            for g in i.diagnosis.genomic_interpretations:
                vcf_record = g.variant_interpretation.variation_descriptor.vcf_record
                genotype = g.variant_interpretation.variation_descriptor.allelic_state
                variant_data = CausativeVariant(
                    file,
                    self.phenopacket_contents.subject.id,
                    vcf_record.genome_assembly,
                    vcf_record.chrom,
                    vcf_record.pos,
                    vcf_record.ref,
                    vcf_record.alt,
                    genotype.label,
                )
                all_variants.append(variant_data)
        return all_variants

    def files(self) -> list:
        """Returns all files associated with a phenopacket."""
        return self.phenopacket_contents.files

    def vcf_file_data(self, vcf_dir: Path) -> typing.NamedTuple:
        """Retrieves the genome assembly and vcf name from a phenopacket."""
        VcfData = namedtuple("VcfData", ["vcf_path", "genome_assembly"])
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
            vcf = vcf_dir.joinpath(vcf_name)
            assembly = file.file_attributes["genomeAssembly"]
            if assembly not in compatible_genome_assembly:
                raise IncompatibleGenomeAssemblyError(assembly, self.phenopacket_contents)
            genome_assembly = assembly
        return VcfData(vcf, genome_assembly)

    def diagnosed_genes(self) -> list[str]:
        """Returns a unique list of all causative genes from a phenopacket."""
        pheno_interpretation = self.interpretations()
        genes = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                genes.append(g.variant_interpretation.variation_descriptor.gene_context.symbol)
        genes = list(set(genes))
        return genes

    def diagnosed_variants(self) -> list:
        """Returns a list of all variants from a phenopacket - for use in assess-prioritisation."""
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
    """Rebuilds a Phenopacket."""

    def __init__(self, phenopacket_contents):
        self.phenopacket_contents = phenopacket_contents

    def add_randomised_hpo(self, randomised_hpo):
        """Adds randomised phenotypic profile to phenopacket."""
        self.phenopacket_contents = Phenopacket(
            id=self.phenopacket_contents.id,
            subject=self.phenopacket_contents.subject,
            phenotypic_features=randomised_hpo,
            interpretations=self.phenopacket_contents.interpretations,
            files=self.phenopacket_contents.files,
            meta_data=self.phenopacket_contents.meta_data,
        )

    def add_created_vcf_path(self, vcf_path: Path, genome_assembly: str):
        phenopacket_files = [
            file for file in self.phenopacket_contents.files if file.file_attributes != "VCF"
        ]
        vcf_file_entry = File(
            uri=os.path.abspath(vcf_path),
            file_attributes={"fileFormat": "VCF", "genomeAssembly": genome_assembly},
        )
        phenopacket_files.append(vcf_file_entry)
        self.phenopacket_contents = Phenopacket(
            id=self.phenopacket_contents.id,
            subject=self.phenopacket_contents.subject,
            phenotypic_features=self.phenopacket_contents.phenotypic_features,
            interpretations=self.phenopacket_contents.interpretations,
            files=phenopacket_files,
            meta_data=self.phenopacket_contents.meta_data,
        )

    def create_json_message(self) -> str:
        """Creates json message for writing to file."""
        if hasattr(self.phenopacket_contents, "proband"):
            family = Family(
                id=self.phenopacket_contents.id,
                proband=self.phenopacket_contents,
                pedigree=self.phenopacket_contents.pedigree,
                files=self.phenopacket_contents.files,
                meta_data=self.phenopacket_contents.meta_data,
            )
            return MessageToJson(family)
        else:
            return MessageToJson(self.phenopacket_contents)

    def write_phenopacket(self, output_file: Path):
        phenopacket = self.create_json_message()
        with open(output_file, "w") as outfile:
            outfile.write(phenopacket)
        outfile.close()
