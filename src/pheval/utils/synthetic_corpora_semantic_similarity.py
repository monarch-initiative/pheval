# Imports
import os
from os import path
from pathlib import Path
import numpy as np
import pandas as pd
import logging
import random
from typing import List, Any, Optional, Union

# Data vis
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Pheval phenopacket utils, Ontology tools, parallel processing
from pheval.utils.phenopacket_utils import phenopacket_reader
from pheval.utils.phenopacket_utils import PhenopacketUtil
from semsimian import Semsimian
import multiprocessing as mp

# Pheval imports
###from pheval.utils.file_utils import ensure_file_exists, ensure_dir_exists


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


def divide_workload(data_list, num_proc: int=1) -> list:
    """
    Meant to divide up the elements in data_list into num_proc equal portions
    by iteratively adding each element to a basket never repeating the same basket until all baskets have an equal amount
    If num_proc == 1 then the original input list will be returned nested in a top layer list i.e. [data_list]
    """

    # Deal with our edge case at the very begginning which then is used as input into the second potential edge case
    ndata_elements = len(data_list)
    if ndata_elements < num_proc:
        num_proc = ndata_elements

    # Edge case
    if num_proc <= 1:
        return [data_list]
    else:
        baskets = [[] for i in range(0, num_proc)]
        index_count = 0
        for d in data_list:
            baskets[index_count].append(d)
            if index_count == (num_proc-1):
                index_count = 0
            else:
                index_count += 1

        #print("- Workload divided into {} portions with each portion recieving {} elements respectively...".format(num_proc, [format(len(b), ',') for b in baskets]))
        return baskets


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


def semsim_termset_compare(hp_db_path, input_args):
    
    # Load in data
    semsim_obj = Semsimian(spo=None, resource_path=hp_db_path)

    # Create output datastructure
    sim_metrics = ["jaccard_similarity", "ancestor_information_content", "phenodigm_score"]
    output_data = {"filename":[],
                   "jaccard_similarity":[],
                   "ancestor_information_content":[],
                   "phenodigm_score":[]}
    
    # For each set of input arguments, compute jaccar, ic, and phenodigm sim scores and add to our output dict
    for i, argset in enumerate(input_args):
        termset1, termset2, filename = argset
        output_data["filename"].append(filename)
        for sim in sim_metrics:
            score = semsim_obj.termset_comparison(set(termset1), set(termset2), score_metric=sim)
            output_data[sim].append(score)

        if i % 100 == 0:
            print("- processed {}/{}".format(i+1, format(len(input_args), ',')))

    return output_data


### Main function to bring everything together
def semsim_phenopackets_comapre(hp_db_path: str,
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
    ensure_file_exists(hp_db_path)
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
    print('wooo', len(comparison_dirs))
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
            pargs = [terms1, terms2, fname]
            parallel_args.append(pargs)
        
        # Divide input arguments into batchs for async processing, and instantiate num_procs worth of Semsimian objects to compute with
        div_parallel_args = divide_workload(parallel_args, num_proc=num_proc)
        div_semsim_objs = [hp_db_path for i in range(0, num_proc)]

        # Setup parallel processing overhead, kick off jobs via asynchronous processing, and retrieve results
        output = mp.Queue()
        pool = mp.Pool(processes=num_proc)
        results = [pool.apply_async(semsim_termset_compare, args=(sem, inargs)) for sem, inargs in zip(div_semsim_objs,
                                                                                                       div_parallel_args)]
        # .get() method needs to be called
        output = [p.get() for p in results]
        print("- Semantic similarity comparisons completed...")
        print("- Combining and writing results to {}...".format(outfilename))

        # Merge results from previous step into single data structure
        # Each element in the output list is a dictionary {sival:[pval,pval,...], }
        merged_data = {}
        for p in output:
            for k,v in p.items():
                if k not in merged_data:
                    merged_data.update({k:[]})
                merged_data[k] += v

        # Convert and write from dataframe
        pd.DataFrame(merged_data).to_csv(outfilename, sep='\t', index=False)
        print("- Data successfully written. {}/{} comparisons made...".format(i, len(comparison_dirs)-1))



def plot_semantic_similarity_results(res_files: List[str], 
                                     comp_labels: List[str] = [], 
                                     save_fig: Optional[Union[str, bool]] = False):
    
    ## make sure data exists before plotting
    #for fpath in res_files:
    #    ensure(file_exists)
    
    # plot individual comparisons
    
    all_data = {}
    sim_keys = {}
    out_figs = []
    for fpath in res_files:
        res_df = pd.read_csv(fpath, sep='\t')
        all_data.update({fpath:{}})
        fig = plt.figure().set_size_inches(12, 4)
        col = 0
        fname = fpath.split('/')[-1].split(".tsv")[0]

        for sim_metric in res_df.columns:
            if sim_metric == "filename":
                continue
            
            plot_vals = np.asarray(res_df[sim_metric].astype(float))
            all_data[fpath].update({sim_metric:plot_vals})
            sim_keys.update({sim_metric:''})
            tot_comps = len(plot_vals)
            avg = round(np.average(plot_vals), 4)
            
            ax = plt.subplot2grid((1, 3), (0, col))
            ax.hist(plot_vals, edgecolor='black')
            ax.set_xlabel("{} score".format(sim_metric))
            ax.set_ylabel("Number of comparisons")
            ax.set_title("{} (Avg={})".format(sim_metric, avg))
            col += 1
        
        plt.suptitle("{} Similarity score results (N={})".format(fname, format(len(res_df), ',')))
        plt.tight_layout()
    
    # Combine comparison to single plot
    if len(res_files) >= 2:
        fig = plt.figure().set_size_inches(12, 4)
        col = 0
        
        for sim_metric in sim_keys:
            ax = plt.subplot2grid((1, 3), (0, col))
            vplot_data =[]
            vlabs = []
            col += 1
            dir_ind = 0
            
            for fpath in all_data:
                fname = fpath.split('/')[-1]
                plot_data = all_data[fpath][sim_metric]
                vplot_data.append(plot_data)
                vlabs.append(fname)
                dir_ind += 1

            
            ax.violinplot(vplot_data) #, labels=vlabs)
            ax.set_ylabel("{} score".format(sim_metric))
            ax.set_xticks([i for i in range(1, len(vplot_data)+1)])
            if len(comp_labels) > 0:
                vlabs = comp_labels
                
            ax.set_xticklabels(vlabs, rotation=0)##, size=labsize)
            ax.set_title("{} similarity scores{}(N={})".format(sim_metric, "\n", format(len(vplot_data[0]), ',')))
        
        plt.tight_layout()
    
    # Save the plots to a PDF file
    if save_fig != False:
        pdf = PdfPages(save_fig)
        fig_nums = plt.get_fignums()
        figs = [plt.figure(n) for n in fig_nums]
        for fig in figs:   
            fig.savefig(pdf, format='pdf') 

        pdf.close()
