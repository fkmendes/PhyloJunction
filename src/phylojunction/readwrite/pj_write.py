import os
import typing as ty
import numpy as np
import pandas as pd # type: ignore

# pj imports
from phylojunction.data.tree import AnnotatedTree
import phylojunction.pgm.pgm as pgm

# for debugging
from tabulate import tabulate # type: ignore

def write_text_output(outfile_handle: ty.IO, content_string_list: ty.List[str]) -> None:
    content_string = "\n".join(content_string_list)
    outfile_handle.write(content_string)


def write_data_df(outfile_handle: ty.IO, data_df: pd.DataFrame, format="csv") -> None:
    """Write a pandas DataFrame to the specified output file stream using the specified format

    Args:
        outfile_handle (file): Output file object to write to
        data_df (pandas.DataFrame): A data frame containing random variable values to print to file
        format (str): Extension for output file

    Returns:
        None
    """
    
    if format == "csv":
        data_df.to_csv(outfile_handle, index=False)
    elif format == "tsv":
        data_df.to_csv(outfile_handle, sep="\t", index=False)


def prep_data_df(node_pgm_list: ty.List[pgm.NodePGM]) -> ty.Tuple[ty.List[str], ty.List[pd.DataFrame]]:
    """Return two pandas DataFrame's, with scalar and tree random variables
    
    Arguments
        node_pgm_list (pgm.NodePGM):

    Returns
        (tuple): Tuple with two lists as elements, one with file-suffix strings, another with pandas.DataFrame's   
    """
    
    data_df_names_list: ty.List[str] = [] # return
    data_df_list: ty.List[pd.DataFrame] = [] # return

    # creating object to hold scalar r.v. output
    scalar_data_df = pd.DataFrame() 

    for node_pgm in node_pgm_list:
        rv_name = node_pgm.node_pgm_name
        node_val = node_pgm.value # list
        n_sim = len(node_pgm)
        n_repl = node_pgm.repl_size
        
        # if stochastic node was sampled and at least one value was indeed saved in list
        if node_pgm.is_sampled and node_val:
            # trees go in different file
            # TODO: move the tree part of the code here to a separate module
            if isinstance(node_val[0], AnnotatedTree):
                # creating object to hold tree node output
                tree_data_df = pd.DataFrame() 
                tree_data_df["simulation"] = np.array([np.repeat(i+1, n_repl) for i in range(n_sim)]).flatten()
                tree_data_df["replicate"] = np.array([[i+1 for i in range(n_repl)] for j in range(n_sim)]).flatten()
                tree_data_df["tree"] = [
                    val.tree.as_string(schema="newick", suppress_annotations=True, suppress_internal_taxon_labels=True).strip("\"").strip() \
                        for val in node_val
                ]
                _nrow_tree_data_df = len(tree_data_df.index)
                
                tree_data_summary_df = pd.DataFrame()
                tree_data_summary_df["simulation"] = tree_data_df["simulation"]
                tree_data_summary_df["replicate"] = tree_data_df["replicate"]
                dummy_list1 = [ 0 for i in range(_nrow_tree_data_df) ]
                dummy_list2 = [ 0.0 for i in range(_nrow_tree_data_df) ]
                tree_data_summary_df["n_total"] = dummy_list1
                tree_data_summary_df["n_extant"] = dummy_list1
                tree_data_summary_df["n_extinct"] = dummy_list1
                tree_data_summary_df["n_sa"] = dummy_list1
                tree_data_summary_df["origin_age"] = dummy_list2
                tree_data_summary_df["root_age"] = dummy_list2
                for idx, val in enumerate(node_val):              
                    tree_data_summary_df.at[idx, "n_total"] = len(val.tree.leaf_nodes()) - val.n_sa # conservative: should match n_extant + n_extinct!
                    tree_data_summary_df.at[idx, "n_extant"] = val.n_extant_terminal_nodes
                    tree_data_summary_df.at[idx, "n_extinct"] = val.n_extinct_terminal_nodes
                    tree_data_summary_df.at[idx, "n_sa"] = val.n_sa
                    tree_data_summary_df.at[idx, "origin_age"] = "{:,.4f}".format(val.origin_age)
                    tree_data_summary_df.at[idx, "root_age"] = "{:,.4f}".format(val.root_age)

                # adding to return
                data_df_list.append(tree_data_df)
                data_df_names_list.append("_trs.tsv")
                data_df_list.append(tree_data_summary_df)
                data_df_names_list.append("_trs_summaries.tsv")

            # can only be scalar at the moment, if not tree
            else:
                scalar_data_df["simulation"] = [i for i in range(n_sim)]
                scalar_data_df[rv_name] = node_val
                
    # we wait until all nodes are seen before adding scalar_data_df to the top of the return list
    data_df_list = [ scalar_data_df ] + data_df_list
    data_df_names_list = ["_rvs.csv"] + data_df_names_list

    # return scalar_data_df, tree_data_df
    return data_df_names_list, data_df_list


def dump_pgm_data(dir_string: str, pgm_obj: pgm.ProbabilisticGraphicalModel, prefix: str="") -> None:
    """Write stochastic-node sampled values in specified directory 
    
    Args:
        dir_string (str): Where to save the files to be written
        pgm_obj (pgm.ProbabilisticGraphicalModel): Probabilistic graphical model instance whose sampled values we are extracting and writing to file
        prefix (str): String to preceed file names

    Returns:
        None
    """

    sorted_node_pgm_list: ty.List[pgm.NodePGM] = pgm_obj.get_sorted_node_pgm_list()

    # populating data df that will be dumped and their file names
    data_df_names_list: ty.List[str]
    data_df_list: ty.List[pd.DataFrame]
    data_df_names_list, data_df_list = prep_data_df(sorted_node_pgm_list) # still prefixless

    # sort out file path
    if not dir_string.endswith("/"):
        dir_string += "/"
    
    output_file_path = dir_string    
    if prefix:
        output_file_path += prefix

    # full file paths
    data_df_full_fp_list = [output_file_path + fn for fn in data_df_names_list]

    # debugging
    # print("data_df_full_fp_list:")
    # print("    " + "\n    ".join(n for n in data_df_full_fp_list))
    # print("_trs_summaries.tsv")
    # print(tabulate(data_df_list[2], headers="keys", tablefmt="pretty"))
    
    # write!
    for idx, full_fp in enumerate(data_df_full_fp_list):
        with open(full_fp, "w") as data_out:
            f: str = ""
            
            if full_fp.endswith("csv"):
                f = "csv"
            elif full_fp.endswith("tsv"):
                f = "tsv"
            
            write_data_df(data_out, data_df_list[idx], format=f)


def get_write_inference_rev_scripts(all_sims_model_spec_list: ty.List[str], all_sims_mcmc_logging_spec_list: ty.List[str], dir_list: ty.List[str], prefix: str="", write2file: bool=False) -> ty.List[str]:
    """Get and/or write full inference .Rev scripts

    Args:
        all_sims_model_spec_list (str): List of strings specifying just the model part of a .Rev script, one element per simulation
        all_sims_mcmc_logging_spec_list (str): List of strings specifying just the MCMC and logging part of a .Rev script, one element per simulation
        dir_list (str): List of three string specifying directories (inference root, scripts, results)
        prefix (str): String prefix to place before the name of files being written
        write2file (bool): If 'True', function writes to file. Defaults to 'False'.

    Returns:
        list of str(s): A list of full .Rev script string specifications, one per simulation
    """
    
    def _iterate_specs_over_sims(n: int):
        """
        Yields:
            (str): String specifying rev script for one simulation
        """
        # i-th sim
        for i in range(n):
            _ith_model_spec = all_sims_model_spec_list[i]
            _ith_mcmc_spec = all_sims_mcmc_logging_spec_list[i]
            _ith_string = _ith_model_spec + "\n" + _ith_mcmc_spec

            yield _ith_string

    _n_sim = len(all_sims_model_spec_list)
    _inference_root_dir, _scripts_dir, _results_dir = dir_list

    _prefix = str()
    if _prefix: _prefix += "_"
    else: _prefix = ""         

    # return
    _all_sims_spec_strs_list = [spec_str for spec_str in _iterate_specs_over_sims(_n_sim)]

    # in addition to returning rev script strings (one per simulation),
    # we also write to file
    if write2file:
        for dir_name in dir_list:
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)

        for i in range(_n_sim):
            _ith_file_name = prefix + "sim_" + str(i) + ".rev"
            _ith_file_path = _scripts_dir + _ith_file_name            
            
            with open(_ith_file_path, "w") as rev_out:
                rev_out.write(_all_sims_spec_strs_list[i])
        
    return _all_sims_spec_strs_list