#!/usr/bin/python

import os
import json
import csv
import click
from phenopackets import Phenopacket, Family
from google.protobuf.json_format import Parse
from collections import defaultdict
from dataclasses import dataclass, field
import pandas as pd
from statistics import mean


class DirectoryFiles:
    """ Class that retrieves an ordered list of relevant files from a directory"""

    def __init__(self, directory: str, file_suffix: str):
        self.directory = directory
        self.file_suffix = file_suffix

    def obtain_files(self) -> list:
        file_list = []
        directory_ = os.fsencode(os.path.join(self.directory, ''))
        for file in os.listdir(directory_):
            filename = os.fsdecode(file)
            if filename.endswith(self.file_suffix):
                file_list.append(filename)
        file_list.sort()
        return file_list


class PhenopacketReader:
    """ Class for reading Phenopackets and retrieving relevant data. """

    def __init__(self, file: str):
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


@dataclass
class PrioritisationRanks:
    """ Class for keeping track of gene ranks for different runs"""
    gene_index: int = 0
    variant_index: int = 0
    directory: str = ""
    phenopacket: str = ""
    gene: str = ""
    variant: str = ""
    gene_rank_comparison = defaultdict(dict)
    variant_rank_comparison = defaultdict(dict)
    rank: int = 0

    def record_gene_rank(self):
        self.gene_index += 1
        self.gene_rank_comparison[self.gene_index]["Phenopacket"] = self.phenopacket
        self.gene_rank_comparison[self.gene_index]["Gene"] = self.gene
        self.gene_rank_comparison[self.gene_index][self.directory] = int(self.rank)

    def record_variant_rank(self):
        self.variant_index += 1
        self.variant_rank_comparison[self.variant_index]["Phenopacket"] = self.phenopacket
        self.variant_rank_comparison[self.variant_index]["Variant"] = self.variant
        self.variant_rank_comparison[self.variant_index][self.directory] = self.rank

    def generate_output(self, prefix: str):
        gene_comparison = pd.DataFrame.from_dict(self.gene_rank_comparison, orient='index')
        gene_comparison['absolute_rank_difference'] = pd.Series.abs(
            gene_comparison.iloc[:, 2] - gene_comparison.iloc[:, 3])
        gene_comparison.to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")
        # gene_comparison_top5 = gene_comparison[
        #     ((gene_comparison.iloc[:, 2:] > 0) & (gene_comparison.iloc[:, 2:] <= 5)).any(axis=1)]
        # gene_comparison_top5.to_csv(prefix + "-top5_gene_rank_comparison.tsv", sep="\t")
        variant_comparison = pd.DataFrame.from_dict(self.variant_rank_comparison, orient='index')
        variant_comparison['absolute_rank_difference'] = pd.Series.abs(
            variant_comparison.iloc[:, 2] - variant_comparison.iloc[:, 3])
        variant_comparison.to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")
        # variant_comparison_top5 = variant_comparison[
        #     ((variant_comparison.iloc[:, 2:] > 0) & (variant_comparison.iloc[:, 2:] <= 5)).any(axis=1)]
        # variant_comparison_top5.to_csv(prefix + "-top5_variant_rank_comparison.tsv", sep="\t")


@dataclass
class RankStats:
    """Class for keeping track of the rank stats."""
    top: int = 0
    top3: int = 0
    top5: int = 0
    found: int = 0
    total: int = 0
    ranks: list = field(default_factory=list)

    def add_rank(self, rank: int):
        self.ranks.append(1 / rank)
        self.found += 1
        if rank == 1:
            self.top += 1
        if rank != '' and rank < 4:
            self.top3 += 1
        if rank != '' and rank < 6:
            self.top5 += 1

    def percentage_rank(self, value: int) -> float:
        return 100 * value / self.found

    def percentage_top(self) -> float:
        return self.percentage_rank(self.top)

    def percentage_top3(self) -> float:
        return self.percentage_rank(self.top3)

    def percentage_top5(self) -> float:
        return self.percentage_rank(self.top5)

    def percentage_found(self) -> float:
        return 100 * self.found / self.total

    def mean_reciprocal_rank(self) -> float:
        return mean(self.ranks)


class RankStatsWriter:

    def __init__(self, file: str):
        self.file = open(file, 'w')
        self.writer = csv.writer(self.file, delimiter='\t')
        self.writer.writerow(['results_directory_path', 'top', 'top3', 'top5', 'found', 'total', 'mean_reciprocal_rank',
                              'percentage_top', 'percentage_top3', 'percentage_top5', 'percentage_found'])

    def write_row(self, directory, rank_stats: RankStats):
        try:
            self.writer.writerow(
                [directory, rank_stats.top, rank_stats.top3, rank_stats.top5, rank_stats.found, rank_stats.total,
                 rank_stats.mean_reciprocal_rank(), rank_stats.percentage_top(), rank_stats.percentage_top3(),
                 rank_stats.percentage_top5(), rank_stats.percentage_found()])
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


@dataclass
class RankingDict:
    original_results: dict
    identifier: str
    output_results: defaultdict
    ranking_method: str

    def add_gene(self):
        self.output_results[self.identifier]["geneSymbol"] = self.original_results["geneIdentifier"]["geneSymbol"]

    def add_ranking_method_val(self):
        self.output_results[self.identifier][self.ranking_method] = round(self.original_results[self.ranking_method], 4)

    # def add_combined_score(self):
    #     try:
    #         self.output_results[self.identifier]["combinedScore"] = round(self.original_results["combinedScore"], 4)
    #     except KeyError:
    #         self.output_results[self.identifier]["combinedScore"] = "N/A"
    #
    # def add_phenotype_score(self):
    #     try:
    #         self.output_results[self.identifier]["phenotypeScore"] = round(self.original_results["phenotypeScore"], 4)
    #     except KeyError:
    #         self.output_results[self.identifier]["phenotypeScore"] = "N/A"
    #
    # def add_variant_score(self):
    #     try:
    #         self.output_results[self.identifier]["variantScore"] = round(self.original_results["variantScore"], 4)
    #     except KeyError:
    #         self.output_results[self.identifier]["variantScore"] = "N/A"
    #
    # def add_pvalue_score(self):
    #     try:
    #         self.output_results[self.identifier]["pValue"] = round(self.original_results["pValue"], 4)
    #     except KeyError:
    #         self.output_results[self.identifier]["pValue"] = "N/A"

    def add_moi(self):
        self.output_results[self.identifier]["modeOfInheritance"] = self.original_results["modeOfInheritance"]

    def add_contributing_variants(self):
        self.output_results[self.identifier]["contributingVariants"] = self.original_results["contributingVariants"]

    def create_ranking_dict(self) -> dict:
        self.add_gene()
        # self.add_combined_score()
        # self.add_phenotype_score()
        # self.add_variant_score()
        # self.add_pvalue_score()
        self.add_ranking_method_val()
        self.add_moi()
        self.add_contributing_variants()
        return self.output_results


class RankResults:
    """ Class for implementing ranking to results - (Exomiser Specific)"""

    def __init__(self, full_path_to_results_file: str, ranking_method: str):
        self.full_path_to_results_file = full_path_to_results_file
        self.ranking_method = ranking_method

    def rank_results(self):
        exomiser_json_result = defaultdict(dict)
        with open(self.full_path_to_results_file) as jsfile:
            js = json.load(jsfile)
            for result in js:
                gene = result["geneScores"]
                for g in gene:
                    if self.ranking_method in g:
                        if "contributingVariants" in g:
                            identifier = g["geneIdentifier"]["geneSymbol"] + "_" + g["modeOfInheritance"]
                            exomiser_json_result_dict = RankingDict(g, identifier, exomiser_json_result,
                                                                    self.ranking_method).create_ranking_dict()
        jsfile.close()
        if self.ranking_method == "pValue":
            exomiser_results = sorted(exomiser_json_result_dict.items(), key=lambda x: x[1][self.ranking_method],
                                      reverse=False)
        else:
            exomiser_results = sorted(exomiser_json_result_dict.items(), key=lambda x: x[1][self.ranking_method],
                                      reverse=True)
        rank, count, previous, result = 0, 0, None, {}
        for key, info in exomiser_results:
            count += 1
            if info[self.ranking_method] != previous:
                rank += count
                previous = info[self.ranking_method]
                count = 0
            result[key] = rank
        ranks = dict(exomiser_results)
        for key, value in result.items():
            ranks[key]["rank"] = value
        return ranks


class AssessmentOfPrioritisation:
    """ Class for assessing gene and variant prioritisation. """

    def __init__(self, ranks: dict, threshold: float, ranking_method: str, prioritisation_ranks: PrioritisationRanks,
                 genes: list, variants: list):
        self.ranks = ranks
        self.threshold = threshold
        self.ranking_method = ranking_method
        self.prioritisation_ranks = prioritisation_ranks
        self.genes = genes
        self.variants = variants

    def assess_gene(self, rank_stats: RankStats):
        for g in self.genes:
            rank_stats.total += 1
            for key, info in self.ranks.items():
                if g == info["geneSymbol"] and float(self.threshold) != 0.0:
                    if self.ranking_method != "pValue" and float(self.threshold) < float(info[self.ranking_method]):
                        rank = info["rank"]
                        self.prioritisation_ranks.rank, self.prioritisation_ranks.gene = rank, g
                        rank_stats.add_rank(rank)
                        break
                    if self.ranking_method == "pValue" and float(self.threshold) > float(info[self.ranking_method]):
                        rank = info["rank"]
                        self.prioritisation_ranks.rank, self.prioritisation_ranks.gene = rank, g
                        rank_stats.add_rank(rank)
                        break
                if g == info["geneSymbol"] and float(self.threshold) == 0.0:
                    rank = info["rank"]
                    self.prioritisation_ranks.rank, self.prioritisation_ranks.gene = rank, g
                    rank_stats.add_rank(rank)
                    break
                self.prioritisation_ranks.rank = 0
            self.prioritisation_ranks.record_gene_rank()

    def assess_variant(self, rank_stats: RankStats):
        for variant in self.variants:
            rank_stats.total += 1
            key, info = variant.items()
            gene = key[1]
            variant = info[1]
            for key1, info1 in self.ranks.items():
                if gene == info1["geneSymbol"]:
                    cont = info1["contributingVariants"]
                    for cv in cont:
                        if variant.chrom == cv["contigName"] and int(cv["start"]) == variant.pos and cv["ref"] == \
                                variant.ref and cv["alt"] == variant.alt and float(self.threshold) != 0.0:
                            if self.ranking_method != "pValue" and float(self.threshold) < float(
                                    info1[self.ranking_method]):
                                rank = info1["rank"]
                                self.prioritisation_ranks.rank, self.prioritisation_ranks.variant = rank, "_".join(
                                    [variant.chrom, str(variant.pos), variant.ref, variant.alt])
                                rank_stats.add_rank(rank)
                                break
                            if self.ranking_method == "pValue" and float(self.threshold) > float(
                                    info1[self.ranking_method]):
                                rank = info1["rank"]
                                self.prioritisation_ranks.rank, self.prioritisation_ranks.variant = rank, "_".join(
                                    [variant.chrom, str(variant.pos), variant.ref, variant.alt])
                                rank_stats.add_rank(rank)
                                break
                        if variant.chrom == cv["contigName"] and int(cv["start"]) == variant.pos and cv["ref"] == \
                                variant.ref and cv["alt"] == variant.alt and float(self.threshold) == 0.0:
                            rank = info1["rank"]

                            self.prioritisation_ranks.rank, self.prioritisation_ranks.variant = rank, "_".join(
                                [variant.chrom, str(variant.pos), variant.ref, variant.alt])
                            rank_stats.add_rank(rank)
                            break
                    break
                self.prioritisation_ranks.rank = 0
            self.prioritisation_ranks.record_variant_rank()


@click.command()
@click.option("--directory-list", "-d",
              required=True, metavar='FILE',
              help="A .txt file containing the full path of results directories,"
                   " with each new directory contained on a new line.")
@click.option("--phenopacket-dir", "-p", required=True, metavar='PATH',
              help="Full path to directory containing phenopackets.")
@click.option("--output-prefix", "-o", metavar='<str>', required=True, help=" Output file prefix. ")
@click.option("--ranking-method", "-r", metavar='<str>', default="combinedScore", show_default=True,
              help="Ranking method for gene prioritisation. Options include: "
                   "combinedScore, phenotypeScore, variantScore and pValue")
@click.option("--threshold", "-t", metavar='<float>', default=float(0.0), required=False, help="Score threshold.")
def assess_prioritisation(directory_list, phenopacket_dir, ranking_method, output_prefix, threshold):
    """ Assesses percentage of top, top3 and top5 genes and variants found Exomiser results
    from confirmed cases in phenopackets. """
    gene_stats_writer = RankStatsWriter(output_prefix + "-gene_summary.tsv")
    variants_stats_writer = RankStatsWriter(output_prefix + "-variant_summary.tsv")
    directories = open(directory_list).read().splitlines()
    prioritisation_ranks = PrioritisationRanks()
    for directory in directories:
        prioritisation_ranks.directory = directory
        gene_rank_stats, variant_rank_stats = RankStats(), RankStats()
        exomiser_json_results = DirectoryFiles(directory, ".json").obtain_files()
        for exomiser_result in exomiser_json_results:
            prioritisation_ranks.phenopacket = exomiser_result
            phenopacket = os.path.join(phenopacket_dir, exomiser_result.replace("-exomiser.json", ".json"))
            interpretations = PhenopacketReader(phenopacket).interpretations()
            genes = PhenopacketReader.diagnosed_genes(interpretations)
            variants = PhenopacketReader.diagnosed_variants(interpretations)
            exomiser_full_path = os.path.join(directory, exomiser_result)
            ranking_dict = RankResults(exomiser_full_path, ranking_method).rank_results()
            assess = AssessmentOfPrioritisation(ranking_dict, threshold, ranking_method, prioritisation_ranks, genes,
                                                variants)
            assess.assess_gene(gene_rank_stats)
            assess.assess_variant(variant_rank_stats)
        gene_stats_writer.write_row(directory, gene_rank_stats)
        variants_stats_writer.write_row(directory, variant_rank_stats)
    prioritisation_ranks.generate_output(output_prefix)
    gene_stats_writer.close()
    variants_stats_writer.close()
