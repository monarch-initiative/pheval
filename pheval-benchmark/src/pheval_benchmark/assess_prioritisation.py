#!/usr/bin/python

import os
import json
import csv
import click
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class RankStats:
    """Class for keeping track of the rank stats."""
    top: int = 0
    top3: int = 0
    top5: int = 0
    found: int = 0
    total: int = 0

    def add_rank(self, rank: int):
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


class RankStatsWriter:

    def __init__(self, file: str):
        self.file = open(file, 'w')
        self.writer = csv.writer(self.file, delimiter='\t')
        self.writer.writerow(['results_directory_path', 'top', 'top3', 'top5', 'percentage_found', 'found', 'total'])

    def write_row(self, directory, rank_stats: RankStats):
        try:
            self.writer.writerow([directory, rank_stats.percentage_top(), rank_stats.percentage_top3(),
                                  rank_stats.percentage_top5(), rank_stats.percentage_found(), rank_stats.found,
                                  rank_stats.total])
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
    json_result_list = []
    directory_ = os.fsencode(directory)
    for file in os.listdir(directory_):
        filename = os.fsdecode(file)
        if filename.endswith(".json"):
            json_result_list.append(filename)
    return json_result_list


def diagnosed_genes(ppacket: str) -> list:
    """ Retrieves a list of unique diagnosed genes contained within the interpretations block of a phenopacket. """
    genomic_interpretations = []
    with open(ppacket) as pfile:
        genes = []
        pheno = json.load(pfile)
        if "interpretations" in pheno["proband"]:
            interpretations = pheno["proband"]["interpretations"]
        else:
            interpretations = pheno["interpretations"]
        for i in interpretations:
            genomic_interpretations = i["diagnosis"]["genomicInterpretations"]
        for g in genomic_interpretations:
            gene = g["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"]
            genes.append(gene)
        genes = list(set(genes))
    pfile.close()
    return genes


def diagnosed_variants(ppacket: str) -> list:
    """ Returns a list of dictionaries containing both gene and variant data. """
    with open(ppacket) as pfile:
        variants = []
        pheno = json.load(pfile)
        if "interpretations" in pheno["proband"]:
            interpretations = pheno["proband"]["interpretations"]
        else:
            interpretations = pheno["interpretations"]
        for i in interpretations:
            genomic_interpretations = i["diagnosis"]["genomicInterpretations"]
        for g in genomic_interpretations:
            gene = defaultdict(dict)
            gene_symbol = g["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"]
            variant = g["variantInterpretation"]["variationDescriptor"]["vcfRecord"]
            gene["geneSymbol"] = gene_symbol
            gene["variant"] = variant
            variants.append(gene)
    pfile.close()
    return variants


def ranking(exomiser_full_path: str, ranking_method: str) -> dict:
    """ Takes in a json file and return`s a dictionary with the default ranking output by
    Exomiser's .genes.tsv file. """
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


def assess_variant(ranks: dict, variants: list, rank_stats: RankStats):
    for variant in variants:
        rank_stats.total += 1
        key, info = variant.items()
        gene = key[1]
        var = info[1]
        for key1, info1 in ranks.items():
            if gene == info1["geneSymbol"]:
                cont = info1["contributingVariants"]
                for cv in cont:
                    if var["chrom"] == cv["contigName"] and int(cv["start"]) == int(var["pos"]) \
                            and cv["ref"] == var["ref"] and cv["alt"] == var["alt"]:
                        rank = info1["rank"]
                        rank_stats.add_rank(rank)
                        break
                break


def assess_gene(ranks: dict, genes: list, rank_stats: RankStats):
    """ Assigns a simple ranking to the gene when found within the json results file.
     Iterates through the json array in order of output from Exomiser."""
    for g in genes:
        rank_stats.total += 1
        for key, info in ranks.items():
            if g == info["geneSymbol"]:
                rank = info["rank"]
                rank_stats.add_rank(rank)
                break


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
def assess_prioritisation(directory_list, phenopacket_dir, ranking_method, output_prefix):
    """ A CLI that given a text file containing a list of result directories output from Exomiser
    returns the percentage of top, top3 and top5 genes found in confirmed cases from phenopackets."""
    gene_stats_writer = RankStatsWriter(output_prefix+"-gene-prioritisation.txt")
    variants_stats_writer = RankStatsWriter(output_prefix+"-variant-prioritisation.txt")
    directories = open(directory_list).read().splitlines()
    for directory in directories:
        gene_rank_stats = RankStats()
        variant_rank_stats = RankStats()
        exomiser_json_results = directory_files(directory)
        for exomiser_result in exomiser_json_results:
            phenopacket = phenopacket_dir + exomiser_result.replace("-exomiser-results.json", ".json")
            genes = diagnosed_genes(phenopacket)
            variants = diagnosed_variants(phenopacket)
            exomiser_full_path = directory + exomiser_result
            ranking_dict = ranking(exomiser_full_path, ranking_method)
            assess_gene(ranking_dict, genes, gene_rank_stats)
            assess_variant(ranking_dict, variants, variant_rank_stats)
        gene_stats_writer.write_row(directory, gene_rank_stats)
        variants_stats_writer.write_row(directory, variant_rank_stats)

    gene_stats_writer.close()
    variants_stats_writer.close()


if __name__ == '__main__':
    assess_prioritisation()
