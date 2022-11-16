import json
import os
from collections import defaultdict
from phenopackets import Phenopacket, Family
from google.protobuf.json_format import Parse
from pheval_benchmark.custom_exceptions import IncompatibleGenomeAssemblyError, IncorrectFileFormatError


class PhenopacketReader:
    """ Class for reading Phenopackets and retrieving relevant data. """

    def __init__(self, file: str):
        self.file_name = file
        file = open(file, "r")
        phenopacket = json.load(file)
        if "proband" in phenopacket:
            self.pheno = Parse(json.dumps(phenopacket), Family())
        else:
            self.pheno = Parse(json.dumps(phenopacket), Phenopacket())
        file.close()

    def interpretations(self) -> list:
        if hasattr(self.pheno, 'proband'):
            pheno_interpretation = self.pheno.proband.interpretations
        else:
            pheno_interpretation = self.pheno.interpretations
        return pheno_interpretation
    
    def variants(self) -> list:
        all_variants = []
        interpretation = self.interpretations()
        for i in interpretation:
            for g in i.diagnosis.genomic_interpretations:
                record = g.variant_interpretation.variation_descriptor.vcf_record
                genotype = g.variant_interpretation.variation_descriptor.allelic_state
                variant = ProbandData(self.pheno.subject.id, record.genome_assembly, "chr" + record.chrom, record.pos,
                                      record.ref, record.alt, genotype.label)
                all_variants.append(variant)
        return all_variants
    
    def files(self) -> list:
        return self.pheno.files

    def vcf_file_data(self, vcf_dir):
        compatible_genome_assembly = ["GRCh37", "hg19", "GRCh38", "hg38"]
        ppacket_files = self.files()
        vcf, genome_assembly = "", ""
        for file in ppacket_files:
            vcf_name = ""
            if file.file_attributes["fileFormat"].upper() == "VCF":
                vcf_path = file.uri
                if "/" in vcf_path:
                    vcf_name = vcf_path.rsplit('/')[-1]
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

    @staticmethod
    def diagnosed_genes(pheno_interpretation: list) -> list:
        genes = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                genes.append(g.variant_interpretation.variation_descriptor.gene_context.symbol)
        genes = list(set(genes))
        return genes

    @staticmethod
    def diagnosed_variants(pheno_interpretation: list) -> list:
        variants = []
        for i in pheno_interpretation:
            for g in i.diagnosis.genomic_interpretations:
                variant = defaultdict(dict)
                variant["geneSymbol"] = g.variant_interpretation.variation_descriptor.gene_context.symbol
                variant["variant"] = g.variant_interpretation.variation_descriptor.vcf_record
                variants.append(variant)
        return variants
