#!/usr/bin/python

import os
import json
import csv
import click
from collections import defaultdict
from dataclasses import dataclass, field
import pandas as pd
from statistics import mean


class PhenopacketReader:

    def __init__(self, file: str):
        file = open(file, "r")
        self.pheno = json.load(file)
        file.close()

    def interpretations(self):
        if "proband" in self.pheno:
            pheno_interpretation = self.pheno["proband"]["interpretations"]
        else:
            pheno_interpretation = self.pheno["interpretations"]
        return pheno_interpretation

    @staticmethod
    def diagnosed_genes(pheno_interpretation):
        genes = []
        genomic_interpretations = []
        for i in pheno_interpretation:
            genomic_interpretations = i["diagnosis"]["genomicInterpretations"]
        for g in genomic_interpretations:
            gene = g["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"]
            genes.append(gene)
        genes = list(set(genes))
        return genes

    @staticmethod
    def diagnosed_variants(pheno_interpretation):
        variants = []
        genomic_interpretations = []
        for i in pheno_interpretation:
            genomic_interpretations = i["diagnosis"]["genomicInterpretations"]
        for g in genomic_interpretations:
            gene = defaultdict(dict)
            gene_symbol = g["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"]
            variant = g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]
            gene["geneSymbol"] = gene_symbol
            gene["variant"] = variant
            variants.append(gene)
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

    def add_gene_rank(self):
        self.gene_index += 1
        self.gene_rank_comparison[self.gene_index]["Phenopacket"] = self.phenopacket
        self.gene_rank_comparison[self.gene_index]["Gene"] = self.gene
        self.gene_rank_comparison[self.gene_index][self.directory] = int(self.rank)

    def add_variant_rank(self):
        self.variant_index += 1
        self.variant_rank_comparison[self.variant_index]["Phenopacket"] = self.phenopacket
        self.variant_rank_comparison[self.variant_index]["Variant"] = self.variant
        self.variant_rank_comparison[self.variant_index][self.directory] = self.rank

    def output_gene_comparison(self, prefix: str):
        gene_comparison = pd.DataFrame.from_dict(self.gene_rank_comparison, orient='index')
        gene_comparison.index.name = "index"
        gene_comparison.to_csv(prefix + "-gene_rank_comparison.tsv", sep="\t")

    def output_variant_comparison(self, prefix: str):
        variant_comparison = pd.DataFrame.from_dict(self.variant_rank_comparison, orient='index')
        variant_comparison.index.name = "index"
        variant_comparison.to_csv(prefix + "-variant_rank_comparison.tsv", sep="\t")

    def output_top5_gene_comparison(self, prefix: str):
        gene_comparison_top5 = pd.DataFrame.from_dict(self.gene_rank_comparison, orient='index')
        gene_comparison_top5 = gene_comparison_top5[
            ((gene_comparison_top5.iloc[:, 2:] > 0) & (gene_comparison_top5.iloc[:, 2:] <= 5)).any(axis=1)]
        gene_comparison_top5 = gene_comparison_top5.reset_index(drop=True)
        gene_comparison_top5.to_csv(prefix + "-top5_gene_rank_comparison.tsv", sep="\t")

    def output_top5_variant_comparison(self, prefix: str):
        variant_comparison_top5 = pd.DataFrame.from_dict(self.variant_rank_comparison, orient='index')
        variant_comparison_top5 = variant_comparison_top5[
            ((variant_comparison_top5.iloc[:, 2:] > 0) & (variant_comparison_top5.iloc[:, 2:] <= 5)).any(axis=1)]
        variant_comparison_top5 = variant_comparison_top5.reset_index(drop=True)
        variant_comparison_top5.to_csv(prefix + "-top5_variant_rank_comparison.tsv", sep="\t")


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

    def mean_reciprocal_rank(self):
        return mean(self.ranks)


class RankStatsWriter:

    def __init__(self, file: str):
        self.file = open(file, 'w')
        self.writer = csv.writer(self.file, delimiter='\t')
        self.writer.writerow(
            ['results_directory_path', 'top', 'percentage_top', 'top3', 'percentage_top3', 'top5', 'percentage_top5',
             'found', 'percentage_found', 'total', "mean_reciprocal_rank"])

    def write_row(self, directory, rank_stats: RankStats):
        try:
            self.writer.writerow(
                [directory, rank_stats.top, rank_stats.percentage_top(), rank_stats.top3, rank_stats.percentage_top3(),
                 rank_stats.top5, rank_stats.percentage_top5(), rank_stats.found, rank_stats.percentage_found(),
                 rank_stats.total, rank_stats.mean_reciprocal_rank()])
        except IOError:
            print("Error writing ", self.file)

    def close(self):
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


def directory_files(directory: str) -> list:
    """ Iterates over all files in a directory and creates a list of all .json results files,
    assumes that the results file end with -exomiser-results.json """
    exomiser_json_result_list = []
    directory_ = os.fsencode(os.path.join(directory, ''))
    for file in os.listdir(directory_):
        filename = os.fsdecode(file)
        if filename.endswith(".json"):
            exomiser_json_result_list.append(filename)
    exomiser_json_result_list.sort()
    return exomiser_json_result_list


def ranking(exomiser_full_path: str, ranking_method: str) -> dict:
    """ Takes in a json file and return`s a dictionary with a specified ranking."""
    exomiser_json_result_dict = defaultdict(dict)
    with open(exomiser_full_path) as jsfile:
        js = json.load(jsfile)
        for result in js:
            gene = result["geneScores"]
            for g in gene:
                if ranking_method in g:
                    if "contributingVariants" in g:
                        identifier = g["geneIdentifier"]["geneSymbol"] + "_" + g["modeOfInheritance"]
                        exomiser_json_result_dict[identifier]["geneSymbol"] = g["geneIdentifier"]["geneSymbol"]
                        exomiser_json_result_dict[identifier][ranking_method] = round(g[ranking_method], 4)
                        exomiser_json_result_dict[identifier]["modeOfInheritance"] = g["modeOfInheritance"]
                        exomiser_json_result_dict[identifier]["contributingVariants"] = g["contributingVariants"]
    if ranking_method == "pValue":
        exomiser_results = sorted(exomiser_json_result_dict.items(), key=lambda x: x[1][ranking_method], reverse=False)
    else:
        exomiser_results = sorted(exomiser_json_result_dict.items(), key=lambda x: x[1][ranking_method], reverse=True)
    rank, count, previous, result = 0, 0, None, {}
    for key, info in exomiser_results:
        count += 1
        if info[ranking_method] != previous:
            rank += count
            previous = info[ranking_method]
            count = 0
        result[key] = rank
    ranks = dict(exomiser_results)
    for key, value in result.items():
        ranks[key]["rank"] = value
    return ranks


def assess_gene(ranks: dict, genes: list, rank_stats: RankStats, threshold, ranking_method,
                gene_ranks: PrioritisationRanks):
    """ Assigns a simple ranking to the gene when found within the json results file.
     Iterates through the json array in order of output from Exomiser."""
    for g in genes:
        rank_stats.total += 1
        for key, info in ranks.items():
            if g == info["geneSymbol"] and float(threshold) != 0.0:
                if ranking_method != "pValue" and float(threshold) < float(info[ranking_method]):
                    rank = info["rank"]
                    gene_ranks.rank, gene_ranks.gene = rank, g
                    rank_stats.add_rank(rank)
                    break
                if ranking_method == "pValue" and float(threshold) > float(info[ranking_method]):
                    rank = info["rank"]
                    gene_ranks.rank, gene_ranks.gene = rank, g
                    rank_stats.add_rank(rank)
                    break
            if g == info["geneSymbol"] and float(threshold) == 0.0:
                rank = info["rank"]
                gene_ranks.rank, gene_ranks.gene = rank, g
                rank_stats.add_rank(rank)
                break
            gene_ranks.rank = 0
        gene_ranks.add_gene_rank()


def assess_variant(ranks: dict, variants: list, rank_stats: RankStats, threshold, ranking_method,
                   variant_ranks: PrioritisationRanks):
    for variant in variants:
        rank_stats.total += 1
        key, info = variant.items()
        gene = key[1]
        variant = info[1]
        for key1, info1 in ranks.items():
            if gene == info1["geneSymbol"]:
                cont = info1["contributingVariants"]
                for cv in cont:
                    if variant["chrom"] == cv["contigName"] and int(cv["start"]) == int(variant["pos"]) \
                            and cv["ref"] == variant["ref"] and cv["alt"] == variant["alt"] and float(threshold) != 0.0:
                        if ranking_method != "pValue" and float(threshold) < float(info1[ranking_method]):
                            rank = info1["rank"]
                            variant_ranks.rank, variant_ranks.variant = rank, variant["chrom"] + "_" + variant[
                                "pos"] + "_" + variant["ref"] + '_' + variant["alt"]
                            rank_stats.add_rank(rank)
                            break
                        if ranking_method == "pValue" and float(threshold) > float(info1[ranking_method]):
                            rank = info1["rank"]
                            variant_ranks.rank, variant_ranks.variant = rank, variant["chrom"] + "_" + variant[
                                "pos"] + "_" + variant["ref"] + '_' + variant["alt"]
                            rank_stats.add_rank(rank)
                            break
                    if variant["chrom"] == cv["contigName"] and int(cv["start"]) == int(variant["pos"]) \
                            and cv["ref"] == variant["ref"] and cv["alt"] == variant["alt"] and float(threshold) == 0.0:
                        rank = info1["rank"]
                        variant_ranks.rank, variant_ranks.variant = rank, variant["chrom"] + "_" + variant[
                            "pos"] + "_" + variant["ref"] + '_' + variant["alt"]
                        rank_stats.add_rank(rank)
                        break
                break
            variant_ranks.rank = 0
        variant_ranks.add_variant_rank()


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
    gene_stats_writer = RankStatsWriter(output_prefix + "-gene_prioritisation.tsv")
    variants_stats_writer = RankStatsWriter(output_prefix + "-variant_prioritisation.tsv")
    directories = open(directory_list).read().splitlines()
    for directory in directories:
        gene_ranks = PrioritisationRanks()
        variant_ranks = PrioritisationRanks()
        gene_ranks.directory, variant_ranks.directory = directory, directory
        gene_rank_stats = RankStats()
        variant_rank_stats = RankStats()
        exomiser_json_results = directory_files(directory)
        for exomiser_result in exomiser_json_results:
            variant_ranks.phenopacket, gene_ranks.phenopacket = exomiser_result, exomiser_result
            phenopacket = os.path.join(phenopacket_dir, exomiser_result.replace("-exomiser.json", ".json"))
            interpretations = PhenopacketReader(phenopacket).interpretations()
            genes = PhenopacketReader.diagnosed_genes(interpretations)
            variants = PhenopacketReader.diagnosed_variants(interpretations)
            exomiser_full_path = os.path.join(directory, exomiser_result)
            ranking_dict = ranking(exomiser_full_path, ranking_method)
            assess_gene(ranking_dict, genes, gene_rank_stats, threshold, ranking_method, gene_ranks)
            assess_variant(ranking_dict, variants, variant_rank_stats, threshold, ranking_method, variant_ranks)
        gene_stats_writer.write_row(directory, gene_rank_stats)
        variants_stats_writer.write_row(directory, variant_rank_stats)
        gene_ranks.output_gene_comparison(output_prefix)
        variant_ranks.output_variant_comparison(output_prefix)
        gene_ranks.output_top5_gene_comparison(output_prefix)
        variant_ranks.output_top5_variant_comparison(output_prefix)
    gene_stats_writer.close()
    variants_stats_writer.close()
