import os
import typing as ty
import numpy as np
import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import statistics as stat
import pickle

# pj imports
from phylojunction.data.tree import AnnotatedTree
import phylojunction.pgm.pgm as pgm

# for debugging
from tabulate import tabulate # type: ignore

def write_text_output(outfile_handle: ty.IO, content_string_list: ty.List[str]) -> None:
    content_string = "\n".join(content_string_list)
    outfile_handle.write(content_string)


def write_fig_to_file(outfile_path: str, fig: plt.Figure) -> None:
    fig.savefig(outfile_path)


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


# TODO: move prep_data_df to separate file, pj_parse 
def prep_data_df(node_pgm_list: ty.List[pgm.NodePGM]) -> ty.Tuple[ty.List[str], ty.List[pd.DataFrame]]:
    """Return two pandas DataFrame's, with scalar and tree random variables
    
    Arguments
        node_pgm_list (pgm.NodePGM):

    Returns
        (tuple): Tuple with two lists as elements, one with file-suffix strings, another with pandas.DataFrame's   
    """
    
    data_df_names_list: ty.List[str] = [] # return
    data_df_list: ty.List[pd.DataFrame] = [] # return

    # creating object to hold r.v. values and summaries
    scalar_data_df = pd.DataFrame()
    scalar_repl_df = pd.DataFrame()
    scalar_repl_summary_df = pd.DataFrame()
    tree_data_df = pd.DataFrame()
    tree_data_summary_df = pd.DataFrame()
    tree_repl_summary_df = pd.DataFrame()

    there_is_scalar_repl_df = False
    there_is_tree_repl_summary_df = False

    # main loop: nodes!
    for node_pgm in node_pgm_list:
        rv_name = node_pgm.node_pgm_name
        node_val = node_pgm.value # list
        n_sim = len(node_pgm)
        n_repl = node_pgm.repl_size
        
        # if stochastic node was sampled and at least one value was indeed saved in list
        if node_pgm.is_sampled and node_val:

            # trees go in separate files
            if isinstance(node_val[0], AnnotatedTree):
                ######################################
                # DataFrame holding trees themselves #
                ######################################
                tree_data_df["sample"] = np.array([np.repeat(i+1, n_repl) for i in range(n_sim)]).flatten()
                tree_data_df["replicate"] = np.array([[i+1 for i in range(n_repl)] for j in range(n_sim)]).flatten()
                tree_data_df["tree"] = [
                    val.tree.as_string(schema="newick", suppress_annotations=True, suppress_internal_taxon_labels=True).strip("\"").strip() \
                        for val in node_val
                ]
                
                ####################################
                # DataFrame holding tree summaries #
                ####################################
                # initialization and some population
                nrow_tree_data_df = len(tree_data_df.index)
                int_dummy_list = np.array([0 for i in range(nrow_tree_data_df)])
                float_dummy_list = np.array([0.0 for i in range(nrow_tree_data_df)])
                tree_data_summary_df["sample"] = tree_data_df["sample"]
                tree_data_summary_df["replicate"] = tree_data_df["replicate"]
                tree_data_summary_df["origin_age"] = float_dummy_list
                tree_data_summary_df["root_age"] = float_dummy_list
                tree_data_summary_df["n_total"] = int_dummy_list
                tree_data_summary_df["n_extant"] = int_dummy_list
                tree_data_summary_df["n_extinct"] = int_dummy_list
                tree_data_summary_df["n_sa"] = int_dummy_list
                
                # deprecated
                # tree_data_summary_df = pd.DataFrame({
                #     "simulation": tree_data_df["simulation"],
                #     "replicate": tree_data_df["replicate"],
                #     "origin_age": float_dummy_list,
                #     "root_age": float_dummy_list,
                #     "n_total": int_dummy_list,
                #     "n_extant": int_dummy_list,
                #     "n_extinct": int_dummy_list,
                #     "n_sa": int_dummy_list,
                # })

                # if there is more than one state, we paste another subdataframe to the main one
                total_n_states = node_val[0].state_count
                if total_n_states > 1:
                    n_in_state_subdf_dict: ty.Dict[str, np.array] = dict()

                    for ith_state in range(total_n_states):
                        n_in_state_subdf_dict["n_" + str(ith_state)] = int_dummy_list
                    
                    n_in_state_subdf = pd.DataFrame(n_in_state_subdf_dict)
                    tree_data_summary_df = pd.concat((tree_data_summary_df, n_in_state_subdf), axis=1)

                # populating tree_data_summary_df
                for idx, val in enumerate(node_val):
                    tree_data_summary_df.at[idx, "root_age"] = "{:,.4f}".format(val.root_age)
                    if val.origin_age:
                        tree_data_summary_df.at[idx, "origin_age"] = "{:,.4f}".format(val.origin_age)     
                    tree_data_summary_df.at[idx, "n_total"] = len(val.tree.leaf_nodes()) - val.n_sa # conservative: should match n_extant + n_extinct!
                    tree_data_summary_df.at[idx, "n_extant"] = val.n_extant_terminal_nodes
                    tree_data_summary_df.at[idx, "n_extinct"] = val.n_extinct_terminal_nodes
                    tree_data_summary_df.at[idx, "n_sa"] = val.n_sa
                    # if we have more than one state, we will have counts for leaves in the different states
                    # NOTE: at the moment, all leaves are counted (in the future, maybe just the sampled ones)
                    if total_n_states > 1:
                        for ith_state in range(val.state_count):
                            tree_data_summary_df.at[idx, "n_" + str(ith_state)] = val.state_count_dict[ith_state]

                # replicates (tree plating)!
                if n_repl > 1:
                    ##############################################
                    # DataFrame holding tree replicate summaries #
                    ##############################################
                    there_is_tree_repl_summary_df = True
                    # initialization and some population
                    float_dummy_list = np.array([0.0 for i in range(2 * n_sim)])
                    tree_repl_summary_df["sample"] = np.array([i+1 for x in range(n_sim) for i in [x, x]]) # 1, 1, 2, 2, 3, 3...
                    tree_repl_summary_df["summary"] = np.array([summary for x in range(n_sim) for summary in ["average", "std. dev."]]) # average, std. dev., average, std. dev....
                    tree_repl_summary_df["origin_age"] = float_dummy_list
                    tree_repl_summary_df["root_age"] = float_dummy_list
                    tree_repl_summary_df["n_total"] = float_dummy_list
                    tree_repl_summary_df["n_extant"] = float_dummy_list
                    tree_repl_summary_df["n_extinct"] = float_dummy_list
                    tree_repl_summary_df["n_sa"] = float_dummy_list

                    repls_all_stats_dict: ty.Dict[str, float] = dict() # all replicates
                    j = 0
                    for i in range(n_sim):
                        start_idx = i * n_repl
                        end_idx = start_idx + n_repl
                        ith_sim_repls = node_val[start_idx:end_idx]

                        repls_all_stats_dict = dict() # all replicates
                        for repl_obj in ith_sim_repls:
                            obj_stat_dict = repl_obj.get_stats_dict()

                            for st, stat_v in obj_stat_dict.items():
                                float_stat_v: float = 0.0
                                
                                try: 
                                    float_stat_v = float(stat_v)
                                except:
                                    # TODO: throw cannot convert to float Exception
                                    pass
                                
                                try:
                                    repls_all_stats_dict[st].append(float_stat_v)
                                
                                except:
                                    repls_all_stats_dict[st] = [ float_stat_v ]

                        if "Origin age" in repls_all_stats_dict:
                            tree_repl_summary_df.at[j, "origin_age"] = stat.mean(repls_all_stats_dict["Origin age"])
                            tree_repl_summary_df.at[j+1, "origin_age"] = stat.stdev(repls_all_stats_dict["Origin age"])
                        tree_repl_summary_df.at[j, "root_age"] = stat.mean(repls_all_stats_dict["Root age"])
                        tree_repl_summary_df.at[j+1, "root_age"] = stat.stdev(repls_all_stats_dict["Root age"])
                        tree_repl_summary_df.at[j, "n_total"] = stat.mean(repls_all_stats_dict["Total taxon count"])
                        tree_repl_summary_df.at[j+1, "n_total"] = stat.stdev(repls_all_stats_dict["Total taxon count"])
                        tree_repl_summary_df.at[j, "n_extant"] = stat.mean(repls_all_stats_dict["Extant taxon count"])
                        tree_repl_summary_df.at[j+1, "n_extant"] = stat.stdev(repls_all_stats_dict["Extant taxon count"])
                        tree_repl_summary_df.at[j, "n_extinct"] = stat.mean(repls_all_stats_dict["Extinct taxon count"])
                        tree_repl_summary_df.at[j+1, "n_extinct"] = stat.stdev(repls_all_stats_dict["Extinct taxon count"])
                        tree_repl_summary_df.at[j, "n_sa"] = stat.mean(repls_all_stats_dict["Direct ancestor count"])
                        tree_repl_summary_df.at[j+1, "n_sa"] = stat.stdev(repls_all_stats_dict["Direct ancestor count"])
                        j += 2
                
                # debugging
                # print(tabulate(tree_repl_summary_df, tablefmt="plain", showindex=False).lstrip())

                ############################
                # Finally adding to return #
                ############################
                data_df_list.append(tree_data_df)
                data_df_names_list.append("trs.tsv")
                data_df_list.append(tree_data_summary_df)
                data_df_names_list.append("trs_summaries.tsv")
                if there_is_tree_repl_summary_df:
                    data_df_list.append(tree_repl_summary_df)
                    data_df_names_list.append("trs_replicate_summary.tsv")

            # can only be scalar at the moment, if not tree
            else:
                scalar_data_df["simulation"] = [i for i in range(n_sim)]

                # single replicate per sample
                if n_repl == 1:
                    scalar_data_df[rv_name] = node_val

                # replicates (plating)!
                else:
                    there_is_scalar_repl_df = True
                    #############################################
                    # DataFrame holding scalar replicate values #
                    #############################################
                    # initialization and some population
                    scalar_repl_df["sample"] = np.array([np.repeat(i+1, n_repl) for i in range(n_sim)]).flatten()
                    scalar_repl_df["replicate"] = np.array([[i+1 for i in range(n_repl)] for j in range(n_sim)]).flatten()
                    
                    ################################################
                    # DataFrame holding scalar replicate summaries #
                    ################################################
                    # initialization and some population
                    scalar_repl_summary_df["sample"] = np.array([i+1 for x in range(n_sim) for i in [x, x]]) # 1, 1, 2, 2, 3, 3...
                    scalar_repl_summary_df["summary"] = np.array([summary for x in range(n_sim) for summary in ["average", "std. dev."]]) # average, std. dev., average, std. dev....
                    
                    # below we populate both replicate-value and replicate-summary dataframes
                    j = 0
                    repls_values_list: ty.List[float] = []
                    repls_avg_sd_list: ty.List[float] = []
                    for i in range(n_sim):
                        start_idx = i * n_repl
                        end_idx = start_idx + n_repl
                        ith_sim_repls = node_val[start_idx:end_idx]
                        repls_values_list += ith_sim_repls
                        repls_avg_sd_list.append(stat.mean(ith_sim_repls))
                        repls_avg_sd_list.append(stat.stdev(ith_sim_repls))
                    
                    scalar_repl_df[rv_name] = repls_values_list # populating replicate-value dataframe
                    scalar_repl_summary_df[rv_name] = repls_avg_sd_list # populating replicate-summary dataframe

    # we wait until all nodes are seen before adding scalar_data_df and potentially scalar_repl_summary_df to the top of the return list
    scalar_dfs_list: ty.List[pd.DataFrame] = [ scalar_data_df ]
    scalar_dfs_names_list: ty.List[pd.str] = [ "rvs.csv" ]
    if there_is_scalar_repl_df:
        scalar_dfs_list += [ scalar_repl_df, scalar_repl_summary_df ]
        scalar_dfs_names_list += [ "rvs_replicates.csv", "rvs_replicate_summary.csv" ]

    data_df_list = scalar_dfs_list + data_df_list
    data_df_names_list = scalar_dfs_names_list + data_df_names_list

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
        output_file_path += prefix + "_"

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


def dump_serialized_pgm(dir_string: str, pgm_obj: pgm.ProbabilisticGraphicalModel, prefix: str="") -> None:
    """Write serialized PGM in specified directory 
    
    Args:
        dir_string (str): Where to save the serialized model to be written
        pgm_obj (pgm.ProbabilisticGraphicalModel): Probabilistic graphical model instance to be serialized and saved
        prefix (str): String to preceed file name

    Returns:
        None
    """
    
    if not dir_string.endswith("/"): dir_string += "/"
    if prefix: prefix += "_"

    with open(dir_string + prefix + "pgm.pickle", "wb") as picklefile:
        pickle.dump(pgm_obj, picklefile)
        

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

    n_sim = len(all_sims_model_spec_list)
    _inference_root_dir, _scripts_dir, _results_dir = dir_list

    _prefix = str()
    if _prefix: _prefix += "_"
    else: _prefix = ""         

    # return
    _all_sims_spec_strs_list = [spec_str for spec_str in _iterate_specs_over_sims(n_sim)]

    # in addition to returning rev script strings (one per simulation),
    # we also write to file
    if write2file:
        for dir_name in dir_list:
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)

        for i in range(n_sim):
            _ith_file_name = prefix + "sim_" + str(i) + ".rev"
            _ith_file_path = _scripts_dir + _ith_file_name            
            
            with open(_ith_file_path, "w") as rev_out:
                rev_out.write(_all_sims_spec_strs_list[i])
        
    return _all_sims_spec_strs_list