#!/usr/bin/python

import os, json, csv, click
from collections import defaultdict

directory_list = click.option("--directory_list", "-d",
                              required=True,metavar='FILE',
                              help="A .txt file containing the full path of results directories,"
                                   " with each new directory contained on a new line.")
path_to_ppackets = click.option("--path_to_ppackets", "-ppackets" ,required=True, metavar='PATH',
                                    help="Full path to directory containing phenopackets.")
outputfile = click.option("--outputfile", "-o", metavar='<str>',required=True, help="Name of output file. ")
ranking_method = click.option("--ranking_method", "-r", metavar='<str>',default="combinedScore",show_default=True,
                              help="Ranking method for gene prioritisation. Options include: "
                                                                                             "combinedScore, phenotypeScore, variantScore and pValue")

def directory_files(directory): # THIS IS KEPT THE SAME WHEN RANKING CHANGES
    """ Iterates over all files in a directory and creates a list of all .json results files. """
    json_result_list=[]
    directory_ = os.fsencode(directory)
    for file in os.listdir(directory_):
        filename = os.fsdecode(file)
        if filename.endswith(".json"):
            json_result_list.append(filename)
    return json_result_list


def unique_gene_list(path_to_ppackets, file): # THIS IS KEPT THE SAME WHEN RANKING CHANGES
    """ Retrieves a list of unique diagnosed genes contained within the interpretations block of a phenopacket. """
    ppacket = path_to_ppackets + file
    genomic_interpretations = []
    with open(ppacket) as pfile:
        genes=[]
        pheno = json.load(pfile)
        try:
            interpretations = pheno["proband"]["interpretations"] # for family
        except KeyError:
            interpretations = pheno["interpretations"]
        for i in interpretations:
            genomic_interpretations = i["diagnosis"]["genomicInterpretations"] # for family
        for g in genomic_interpretations:
            gene = g["variantInterpretation"]["variationDescriptor"]["geneContext"]["symbol"] # for family
            genes.append(gene)
        genes = list(set(genes))
    pfile.close()
    return(genes)

def ranking(json_file, directory, ranking_method): # this same ranking I guess but for variantScore, phenotypeScore, pValue too
    """ Takes in a json file and returns a dictionary with the default ranking output by Exomiser's .genes.tsv file. """
    jsondict = defaultdict(dict)
    json_file = directory + json_file
    with open(json_file) as jsfile:
        js = json.load(jsfile)
        for result in js:
            gene = result["geneScores"]
            for g in gene:
                if ranking_method in g:
                    if "contributingVariants" in g:
                        identifier = g["geneIdentifier"]["geneSymbol"] + "_" + g["modeOfInheritance"]
                        jsondict[identifier]["geneSymbol"] = g["geneIdentifier"]["geneSymbol"]
                        jsondict[identifier][ranking_method] = round(g[ranking_method], 4)
                        jsondict[identifier]["modeOfInheritance"] = g["modeOfInheritance"]
                        jsondict[identifier]["contributingVariants"] = g["contributingVariants"]
    if ranking_method == "pValue":
        res = sorted(jsondict.items(), key=lambda x: x[1][ranking_method], reverse=False)
    else:
        res = sorted(jsondict.items(), key = lambda x: x[1][ranking_method], reverse=True)
    rank, count, previous, result = 0, 0, None, {}
    for key, info in res:
        count += 1
        if info[ranking_method] != previous:
            rank += count
            previous = info[ranking_method]
            count = 0
        result[key] = rank
    res = dict(res)
    for key, value in result.items():
        res[key]["rank"] = value
    return res

def assign_ranking(ranking,genes,found,total,top,top3,top5): # THIS WILL CHANGE FOR DIFFERENT RANKING COULD POTENTIALLY JUST CREATE A WHOLE DIFFERENT FUNCTION OR IF STATEMENTS - unsure until I figure it out
    """ Assigns a simple ranking to the gene when found within the json results file.
     Iterates through the json array in order of output from Exomiser."""
    rank = 100000000
    for g in genes:
        found += 1
        for key, info in ranking.items():
            if g == info["geneSymbol"]:
                rank = info["rank"]
                break
        if rank == 1:
            top += 1
        if rank != '' and rank < 4:
            top3 += 1
        if rank != '' and rank < 6:
            top5 += 1
        total += 1
    return top, top3, top5, found, total

@click.command()
@directory_list
@path_to_ppackets
@ranking_method
@outputfile
def generate_output(directory_list, path_to_ppackets, ranking_method, outputfile): # THIS WILL ALSO CHANGE FOR DIFFERENT RANKING
    """ A CLI that given a text file containing a list of result directories output from Exomiser
    returns the percentage of top, top3 and top5 genes found in confirmed cases from phenopackets."""
    with open(outputfile , 'w') as out:
        tsv_writer = csv.writer(out, delimiter='\t')
        tsv_writer.writerow(['results_directory_path', 'top', 'top3', 'top5', 'found', 'total'])
    out.close()
    directory_list = open(directory_list).read().splitlines()
    for directory in directory_list:
        top=0
        top3=0
        top5=0
        total=0
        found=0
        json_files = directory_files(directory)
        for file in json_files:
            genes2 = unique_gene_list(path_to_ppackets, file)
            default_ranking_dict = ranking(file, directory,ranking_method)
            top, top3, top5, found, total = assign_ranking(default_ranking_dict, genes2, found,total,top,top3,top5)
        top_p = 100*top/found
        top3_p = 100*top3/found
        top5_p = 100*top5/found
        with open(outputfile , 'a') as out:
            tsv_writer = csv.writer(out, delimiter='\t')
            tsv_writer.writerow([directory, top_p, top3_p, top5_p, found, total])
        out.close()



if __name__ == '__main__':
    generate_output()



