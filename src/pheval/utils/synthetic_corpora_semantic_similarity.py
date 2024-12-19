# Imports
import os
from os import path
from pathlib import Path
import numpy as np
import pandas as pd
import logging
import random

from pydantic import BaseModel
from typing import List, Any, Optional

# Pheval phenopacket utils
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil

# Ontology tools
from oaklib import OntologyResource, get_adapter
from oaklib.implementations import ProntoImplementation

# Parallel process lib
from multiprocessing import Pool

# Pheval base tools imports
from utils import ensure_file_exists, ensure_dir_exists


### Utility functions
def ensure_file_exists(*files: str):
    """Ensures the existence of files passed as parameter
    Raises:
        FileNotFoundError: If any file passed as a parameter doesn't exist a FileNotFound Exception will be raised
    """
    for file in files:
        if not path.isfile(file):
            raise FileNotFoundError(f"File {file} not found")


def ensure_dir_exists(dirs: List[str]):
    """Ensures the existence of directories passed as parameter
    Raises:
        DirectoryNotFoundError: If any file passed as a parameter doesn't exist a FileNotFound Exception will be raised
    """
    for dirname in dirs:
        if not path.isdir(dirname):
            raise FileNotFoundError(f"Directory {dirname} not found")


### Semantic similariy code
def avg_best_matches(termset_sim_obj):
    scores = {"jaccard_similarity":[], "ancestor_information_content":[], "phenodigm_score":[]}
    for subj_id, best_match_obj in termset_sim_obj.subject_best_matches.items():
        for score_type in scores:
            sim_score = getattr(best_match_obj.similarity, score_type)
            scores[score_type].append(sim_score)
    
    # Avg scores
    for k in scores:
        scores[k] = np.average(scores[k])
    
    return scores


def load_ontology_and_compare_terms(input_args):
    
    # Expand input arguments and create output datastructure
    onto_file, termset1, termset2, filename = input_args
    output_data = {"filename":[],
                   "jaccard_similarity":[],
                   "ancestor_information_content":[],
                   "phenodigm_score":[]}

    resource = ProntoImplementation(OntologyResource(slug=str(onto_file), local=True))
    dv = resource.termset_pairwise_similarity(termset1, termset2)
    tscores = avg_best_matches(dv)
            
    output_data["filename"].append(filename)
    for sim_metric in output_data:
        if sim_metric == "filename":
            continue
        output_data[sim_metric].append(tscores[sim_metric])
    
    return output_data


### Main function to bring everything together
def oaklib_phenopackets_comapre(hpo_obo_filepath: str,
                                comparison_dirs: List[str], 
                                output_dir: str, 
                                output_prefix: str = "semantic_compare",
                                num_proc: int = 1):
    """
    Expects a list of directories as input where each directory is a set of phenopackets. Each set of phenopackets should belong
    to the same corpora (i.e. share the same set of file names). The underlying assumption here is that the pheontype terms are 
    different or "scrambled" for all directories except the first in the list.
    
    The first input directory within the argument will be compared against all others. For example...
     - If input directories are A, B, C, D... Then A-->B, A-->C, A-->D comparisons will be made
     - If input directories are B, C, D... Then B-->C, B-->D comparisons will be made

    """
    
    # Filepath and input arg checks
    ensure_file_exists(hpo_obo_filepath)
    ensure_dir_exists(comparison_dirs)
    ensure_dir_exists([output_dir])
    if len(comparison_dirs) < 2:
        raise ValueError("- Number of input directories must be greater than 1... Found {}".format(comparison_dirs))
    
    # Ensure finenames match between base directory, and all others (for simetric comparisons)
    base_fnames = set([fname for fname in os.listdir(comparison_dirs[0]) if fname.endswith(".json")])
    for cdir in comparison_dirs[1:]:
        comp_fnames = set([fname for fname in os.listdir(cdir) if fname.endswith(".json")])
        if comp_fnames != base_fnames:
            raise ValueError("- Input file set {} does not match file set found in {}".format(comparison_dirs[0], cdir))
    
    print("- File path checks succesfull...")
    print("- Loading phenotype information from input phenopackets...")
    tot_phenopackets = len(base_fnames)
    tot_dirs = len(comparison_dirs)
    
    # Extract phenotype terms as list for all sets of phenopackets
    corpora_phenotype_terms = {}
    for i, cname in enumerate(comparison_dirs):
        corpora_phenotype_terms.update({i:{fname:[] for fname in base_fnames}})
        for fname in base_fnames:
            
            # Open phenopacket file and pull relevant information (we want observed phenotypes)
            fpath = os.path.join(cname, fname)
            phenopacket_util = PhenopacketUtil(phenopacket_reader(fpath))
            observed_phenotypes = phenopacket_util.observed_phenotypic_features()
            phenotype_ids = [observed_phenotype.type.id for observed_phenotype in observed_phenotypes]
            
            # Fill in our corpora-->sample-->phenotype datastructure
            corpora_phenotype_terms[i][fname] = phenotype_ids

    print("- Phenotype terms extracted from {} input phenopackets...".format(format(tot_phenopackets, ',')))
    print("- Computing {} sets of semantic similarity scores across {} directories...".format(format(tot_phenopackets,','), tot_dirs))
    
    # Hard code our similarity metrics here
    sim_keys = ["jaccard_similarity", "ancestor_information_content", "phenodigm_score"]

    # HPO termset similarity comparisons between base phenopacket terms and new
    for i in range(1, len(comparison_dirs)):
        print("- Comparing directorie indices {}-->{}...".format(0, i))
        
        # Grab directory names for output filename formatting
        dname_base = comparison_dirs[0].split('/')[-1]
        dname_comp = comparison_dirs[i].split('/')[-1]
        outfilename = os.path.join(output_dir, "{}_{}_VS_{}.tsv".format(output_prefix, dname_base, dname_comp))
        
        # Create input arguments for parallel "map" function
        # Should be [[ontology_filpath, phenotype_list1, phenotype_list2, common_filename], []...]
        parallel_args = []
        for fname in base_fnames:
            terms1, terms2 = corpora_phenotype_terms[0][fname], corpora_phenotype_terms[i][fname]
            pargs = [hpo_obo_filepath, terms1, terms2, fname]
            parallel_args.append(pargs)
        
        # Perform comparisons in parallel
        pool = Pool(processes=num_proc)
        res = pool.map_async(load_ontology_and_compare_terms, parallel_args).get()
        print("- Semantic similarity comparisons completed...")
        print("- Combining and writing results to {}...".format(outfilename))

        # Combine and write results
        df_keys = list(res[0].keys())
        combined_df = {k:[] for k in df_keys}
        for r in res:
            for k in df_keys:
                combined_df[k] += r[k]
        
        pd.DataFrame(combined_df).to_csv(outfilename, sep='\t', index=False)
        print("- Data successfully written. {}/{} comparisons made...".format(i, len(comparison_dirs)-1))


def plot_semantic_similarity_results(res_files: List[str]):
    
    # make sure data exists before plotting
    for fpath in res_files:
        ensure_file_exists(fpath)
    
    # plot individual comparisons
    all_data = {}
    sim_keys = {}
    for fpath in res_files:
        res_df = pd.read_csv(fpath, sep='\t')

        all_data.update({fpath.split('/')[-1]:{}})
        fig = plt.figure().set_size_inches(12, 4)
        col = 0

        for sim_metric in df.columns:
            if sim_metric == "filename":
                continue
            
            plot_vals = np.asarray(df[sim_metric].astype(float))
            all_data[fpath].update({sim_metric:plot_vals})
            sim_keys.update({sim_metric:''})
            tot_comps = len(plot_vals)
            avg = round(np.average(plot_vals), 4)
            
            ax = plt.subplot2grid((1, 3), (0, col))
            ax.hist(plot_vals, edgecolor='black')
            ax.set_xlabel("{} score".format(sim_metric))
            ax.set_ylabel("Number of comparisons")
            ax.set_title("{} similarity scores{}(N={}, Avg={})".format(sim_metric, '\n', format(tot_comps, ','), avg))
            col += 1
    
        plt.tight_layout()
        plt.show()
    

    # Combine comparison to single plot
    if len(res_files) >= 2:
        fig = plt.figure().set_size_inches(12, 4)
        col = 0
        
        for sim_metric in sim_keys:
            ax = plt.subplot2grid((1, 3), (0, col))
            col += 1
            for fpath in all_data:
                fname = fpath.split('/')[-1]
                plot_data = all_data[fpath][sim_metric]
                ax.hist(plot_data, edgecolors='black', density=True, lable=fname, alpha=.5)
            
            ax.set_xlabel("{} score".format(sim_metric))
            ax.set_ylabel("Number of comparisons")
            ax.set_title("{} similarity scores comparisons)".format(sim_metric))
        
        plt.tight_layout()
        plt.show()
