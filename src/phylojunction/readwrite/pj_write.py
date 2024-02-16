import os
import typing as ty
import numpy as np
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import statistics as stat
import pickle

# pj imports
from phylojunction.data.tree import AnnotatedTree
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmdp

# for debugging
from tabulate import tabulate  # type: ignore


def write_str_list(
        outfile_handle: ty.IO,
        content_string_list: ty.List[str]) -> None:
    content_string = "\n".join(content_string_list)
    outfile_handle.write(content_string)


def write_fig_to_file(
        outfile_path: str,
        fig_obj: plt.Figure) -> None:
    fig_obj.savefig(outfile_path, bbox_inches="tight", dpi=300)


def write_data_df(
        outfile_handle: ty.IO,
        data_df: pd.DataFrame,
        format="csv") -> None:
    """Write a pandas DataFrame to output file stream

    Args:
        outfile_handle (file): Output file object to write to
        data_df (pandas.DataFrame): A data frame containing random
            variable values to print to file
        format (str): Extension for output file

    Returns:
        None
    """

    if format == "csv":
        data_df.to_csv(outfile_handle, index=False)
    elif format == "tsv":
        data_df.to_csv(outfile_handle, sep="\t", index=False)


def initialize_scalar_dataframe(
        sample_size: int,
        n_repl: int = 1,
        summaries_avg_over_repl: bool = False) -> pd.DataFrame:
    """_summary_

    Args:
        sample_size (int): Number of samples (i.e., independent full
            model simulations)
        n_repl (int, optional):  How many times the scalar variable
            is replicated. Defaults to 1.
        summaries_avg_over_repl (bool, optional): Dataframe will hold
            statistics (avg., st. dev.) summarized over replicates.
            Defaults to False.

    Returns:
        pd.DataFrame:  DataFrame with certain columns holding 0 or 0.0 values
    """

    return_df = pd.DataFrame()

    # summarizing (avg and stdev, which is why the 2*) over replicates
    if summaries_avg_over_repl:
        return_df["sample"] = \
            np.array([(i + 1) for i in range(sample_size * 2)]).flatten()
        return_df["summary"] = \
            np.array([summary for x in range(sample_size) for
                      summary in ["average", "std. dev."]])

    # scalar value for all replicates and all samples
    else:
        return_df["sample"] = \
            np.array([np.repeat(i + 1, n_repl)
                      for i in range(sample_size)]).flatten()
        return_df["replicate"] = \
            np.array([[i + 1 for i in range(n_repl)]
                      for j in range(sample_size)]).flatten()

    return return_df


def initialize_tree_dataframe(
        sample_size: int,
        n_repl: int = 1,
        summaries: bool = False,
        summaries_avg_over_repl: bool = False) -> pd.DataFrame:
    """Initialize pandas DataFrame to hold tree information

    Args:
        sample_size (int): Number of samples (i.e., independent full
            model simulations)
        n_repl (int, optional): How many times the tree is replicated.
            Defaults to 1.
        summaries (bool, optional): Dataframe will hold individual-tree
            statistics. Defaults to False.
        summaries_avg_over_repl (bool, optional): Dataframe will hold
            statistics summarized over replicates. Defaults to False.

    Returns:
        pd.DataFrame: DataFrame with certain columns holding 0 or
            0.0 values
    """

    return_df = pd.DataFrame()

    # tree summaries
    if summaries:
        # summaries for all replicated trees for all samples
        if not summaries_avg_over_repl:
            float_dummy_list = np.array(
                [0.0 for i in range(n_repl * sample_size)])
            int_dummy_list = np.array(
                [0 for i in range(n_repl * sample_size)])
            return_df["sample"] = \
                np.array([np.repeat(i + 1, n_repl)
                          for i in range(sample_size)]).flatten()
            return_df["replicate"] = \
                np.array([[i + 1 for i in range(n_repl)]
                          for j in range(sample_size)]).flatten()
            return_df["origin_age"] = float_dummy_list
            return_df["root_age"] = float_dummy_list
            return_df["n_total"] = int_dummy_list
            return_df["n_extant"] = int_dummy_list
            return_df["n_extinct"] = int_dummy_list
            return_df["n_sa"] = int_dummy_list

        # summarizing stats over replicates
        else:
            float_dummy_list = np.array(
                [0.0 for i in range(2 * sample_size)])
            return_df["sample"] = \
                np.array([i + 1 for x in range(sample_size)
                          for i in [x, x]])  # 1, 1, 2, 2, 3, 3...
            # average, std. dev., average, std. dev....
            return_df["summary"] = \
                np.array([summary for x in range(sample_size)
                          for summary in ["average", "std. dev."]])
            return_df["origin_age"] = float_dummy_list
            return_df["root_age"] = float_dummy_list
            return_df["n_total"] = float_dummy_list
            return_df["n_extant"] = float_dummy_list
            return_df["n_extinct"] = float_dummy_list
            return_df["n_sa"] = float_dummy_list

    # tree newick strings, for all replicates and all samples
    else:
        return_df["sample"] = \
            np.array([np.repeat(i + 1, n_repl)
                      for i in range(sample_size)]).flatten()
        return_df["replicate"] = \
            np.array([[i + 1 for i in range(n_repl)]
                      for j in range(sample_size)]).flatten()

    return return_df


def prep_data_df(
        dag_obj: pgm.DirectedAcyclicGraph,
        write_nex_states: bool = False) -> \
    ty.Tuple[ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]],
             ty.List[ty.Dict[str, pd.DataFrame]]]:
    """Return two pandas DataFrame's, with scalar and tree random variables.

    Args:
        dag_obj (DirectedAcyclicGraph): DAG object holding the
            simulated data to be tabulated.
        write_nex_states (bool): Whether to write .nex file with states

    Returns:
        (tuple): Tuple with two lists as elements, one with file-suffix
                 strings, another with pandas.DataFrame's
    """

    sample_size = dag_obj.sample_size
    node_dag_list = dag_obj.get_sorted_node_dag_list()

    # tree nodes (one dataframe per different tree node)

    # { "[tr_node_name1].tsv": pd.DataFrame with newick strings,
    # for all replicates and samples, "[tr_node_name2].tsv": ... }
    tree_value_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # { "[tr_node_name1].tsv": pd.DataFrame with newick strings
    # (reconstructed trees), for all replicates and samples,
    # "[tr_node_name2].tsv": ... }
    tree_rec_value_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # same as above, but with branches annotated with states
    tree_ann_value_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # same as above, but with branches annotated with states
    tree_rec_ann_value_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # { "[tr_node_name1]_summaries.tsv": pd.DataFrame with summary stats,
    # for all replicates and samples }
    tree_summary_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # { "[tr_node_name1]_replicate_summary.tsv": pd.DataFrame summary stats
    # averaged over replicates, for all samples }
    tree_repl_summary_df_dict: ty.Dict[str, pd.DataFrame] = dict()

    # { "[tr_node_name1]_sampleidx_repl_idx": states_str, ... }
    tree_living_nd_states_str_dict: ty.Dict[str, str] = dict()

    # { "[tr_node_name1]_sampleidx_repl_idx": nexus_str, ... }
    tree_living_nd_states_str_nexus_dict: ty.Dict[str, str] = dict()

    # { "[tr_node_name1]_sampleidx_repl_idx_anc_states": states_str, ... }
    tree_internal_nd_states_str_dict: ty.Dict[str, str] = dict()

    # creating object to hold r.v. values and summaries
    scalar_constant_df: pd.DataFrame = pd.DataFrame()

    # { # of replicates: pd.DataFrame summary stats averaged over
    # replicates, for all samples ... }
    scalar_value_df_dict: ty.Dict[int, pd.DataFrame] = dict()

    scalar_repl_summary_df: pd.DataFrame = pd.DataFrame()

    # NOTE: we have a separate dataframe for each tree whose
    # replicates are being summarized, but a single one for scalars
    # because a user might have trees with very different summary
    # stats; scalars always have only average and std. dev.

    # main loop: nodes!
    for node_dag in node_dag_list:
        rv_name = node_dag.node_name
        node_val = node_dag.value  # list
        n_repl = node_dag.repl_size

        # if stochastic node is constant and at least one value was indeed
        # saved in list
        if isinstance(node_dag, pgm.StochasticNodeDAG) and node_val:
            ################################
            # Fixed-value stochastic nodes #
            ################################

            # scalar constants (no support for replication via 2D-lists yet)
            if isinstance(node_val[0], (str, int, float, np.float64)) \
                    and not node_dag.is_sampled:

                if scalar_constant_df.empty:
                    if sample_size == 0:
                        sample_size = len(node_val)

                    scalar_constant_df = \
                        initialize_scalar_dataframe(
                            sample_size,
                            n_repl=1,
                            summaries_avg_over_repl=False)

                if len(node_val) > 1:
                    scalar_constant_df[rv_name] = node_val

                # vectorize single value up to sample size
                else:
                    scalar_constant_df[rv_name] = \
                        [node_val[0] for i in range(sample_size)]

            ############################
            # Sampled stochastic nodes #
            ############################

            ###########
            # Scalars #
            ###########
            if isinstance(node_val[0], (str, int, float, np.float64)) \
                    and node_dag.is_sampled:

                ######################################################
                # DataFrame holding scalar values (replicates if     #
                # there are any), one per total number of replicates #
                ######################################################
                try:
                    scalar_value_df_dict[n_repl].empty

                except (KeyError, AttributeError, NameError) as e:
                    scalar_value_df_dict[n_repl] = \
                        initialize_scalar_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries_avg_over_repl=False)

                ################################################
                # DataFrame holding scalar replicate summaries #
                ################################################
                if n_repl > 1 and scalar_repl_summary_df.empty:
                    scalar_repl_summary_df = \
                        initialize_scalar_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries_avg_over_repl=True)

                # below we populate both replicate-value and
                # replicate-summary dataframes
                j = 0
                repls_values_list: ty.List[float] = []
                repls_avg_sd_list: ty.List[float] = []
                for i in range(sample_size):
                    start_idx = i * n_repl
                    end_idx = start_idx + n_repl
                    ith_sim_repls = node_val[start_idx:end_idx]
                    repls_values_list += ith_sim_repls

                    if n_repl > 1:
                        repls_avg_sd_list.append(stat.mean(ith_sim_repls))
                        sd = stat.stdev(ith_sim_repls)

                        if sd < 1E-10:
                            repls_avg_sd_list.append(0.0)

                        else:
                            repls_avg_sd_list.append(sd)

                # populating replicate-value dataframe and storing it
                # for its number of replicates
                scalar_value_df_dict[n_repl][rv_name] = repls_values_list

                if n_repl > 1:
                    # populating replicate-summary dataframe
                    scalar_repl_summary_df[rv_name] = repls_avg_sd_list

            #########
            # Trees #
            #########
            elif isinstance(node_val[0], AnnotatedTree):

                ######################################
                # DataFrame holding trees themselves #
                # (one per tree node)                #
                ######################################
                try:
                    tree_value_df_dict[rv_name].empty

                except (KeyError, AttributeError, NameError) as e:
                    tree_value_df_dict[rv_name] = \
                        initialize_tree_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries=False,
                            summaries_avg_over_repl=False)

                    tree_ann_value_df_dict[rv_name] = \
                        initialize_tree_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries=False,
                            summaries_avg_over_repl=False)

                    tree_rec_value_df_dict[rv_name] = \
                        initialize_tree_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries=False,
                            summaries_avg_over_repl=False)

                    tree_rec_ann_value_df_dict[rv_name] = \
                        initialize_tree_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries=False,
                            summaries_avg_over_repl=False)

                # complete trees #
                tree_value_df_dict[rv_name][rv_name] = [
                    ith_val.tree.as_string(
                        schema="newick",
                        suppress_annotations=True,
                        suppress_internal_taxon_labels=True,
                        suppress_internal_node_labels=False).strip("\"")
                    .strip() for ith_val in node_val
                ]

                tree_ann_value_df_dict[rv_name][rv_name] = [
                    ith_val.tree.as_string(
                        schema="newick",
                        suppress_annotations=False,
                        suppress_internal_taxon_labels=True,
                        suppress_internal_node_labels=True).strip("\"")
                    .strip() for ith_val in node_val
                ]

                # print("complete tree = " + node_val[0].tree.as_string(
                #         schema="newick",
                #         suppress_annotations=True,
                #         suppress_internal_taxon_labels=True,
                #         suppress_internal_node_labels=True).strip("\""))

                # reconstructed trees #
                # suppress_rooting=False guarantees that
                # [&R] is not added to tree Newick string as
                # a result of re-rooting the reconstructed tree
                rec_tree_list = []
                rec_tree_ann_list = []
                for ith_val in node_val:
                    ith_rec_tree = ith_val.extract_reconstructed_tree()

                    rec_tree_list.append(
                        ith_rec_tree.as_string(
                            schema="newick",
                            suppress_annotations=True,
                            suppress_internal_taxon_labels=True,
                            suppress_rooting=True,
                            suppress_internal_node_labels=False)
                        .strip("\"").strip()
                    )

                    rec_tree_ann_list.append(
                        ith_rec_tree.as_string(
                            schema="newick",
                            suppress_annotations=False,
                            suppress_internal_taxon_labels=True,
                            suppress_rooting=True,
                            suppress_internal_node_labels=False)
                        .strip("\"").strip()
                    )

                tree_rec_value_df_dict[rv_name][rv_name] = rec_tree_list
                tree_rec_ann_value_df_dict[rv_name][rv_name] = \
                    rec_tree_ann_list

                ####################################
                # DataFrame holding tree summaries #
                # (one per tree node)              #
                ####################################
                try:
                    tree_summary_df_dict[rv_name].empty

                except (KeyError, AttributeError, NameError) as e:
                    tree_summary_df_dict[rv_name] = \
                        initialize_tree_dataframe(
                            sample_size,
                            n_repl=n_repl,
                            summaries=True,
                            summaries_avg_over_repl=False)

                total_n_states = node_val[0].state_count

                # populating tree_summary_df_dict,
                # iterating over samples and replicates
                idx = 0
                for replicate_tree in node_val:
                    tree_summary_df_dict[rv_name].at[idx, "root_age"] = \
                        float("{:,.4f}".format(replicate_tree.root_age))

                    if replicate_tree.origin_age:
                        tree_summary_df_dict[rv_name].at[idx, "origin_age"] = \
                            float("{:,.4f}".format(replicate_tree.origin_age))

                    # conservative: should match n_extant + n_extinct!
                    tree_summary_df_dict[rv_name].at[idx, "n_total"] = \
                        len(replicate_tree.tree.leaf_nodes()) \
                        - replicate_tree.n_sa_nodes
                    tree_summary_df_dict[rv_name].at[idx, "n_extant"] = \
                        replicate_tree.n_extant_terminal_nodes
                    tree_summary_df_dict[rv_name].at[idx, "n_extinct"] = \
                        replicate_tree.n_extinct_terminal_nodes
                    tree_summary_df_dict[rv_name].at[idx, "n_sa"] = \
                        replicate_tree.n_sa_nodes

                    # if we have more than one state, we will have counts
                    # for leaves in the different states
                    # NOTE: at the moment, all leaves are counted (in the
                    # future, maybe just the sampled ones)
                    if total_n_states > 1:
                        for ith_state in range(replicate_tree.state_count):
                            tree_summary_df_dict[rv_name].at[idx, "n_" + str(ith_state)] = \
                                replicate_tree.state_count_dict[ith_state]

                    idx += 1

                ##############################################
                # DataFrame holding tree replicate summaries #
                # (one per tree node)                        #
                ##############################################
                if n_repl > 1:
                    try:
                        tree_repl_summary_df_dict[rv_name].empty

                    except (KeyError, AttributeError, NameError) as e:
                        tree_repl_summary_df_dict[rv_name] = \
                            initialize_tree_dataframe(
                                sample_size,
                                n_repl=n_repl,
                                summaries=True,
                                summaries_avg_over_repl=True)

                    # all replicates
                    repls_all_stats_dict: ty.Dict[str, float] = dict()
                    j = 0
                    for i in range(sample_size):
                        start_idx = i * n_repl
                        end_idx = start_idx + n_repl
                        ith_sim_repls = node_val[start_idx:end_idx]

                        repls_all_stats_dict = dict()  # all replicates
                        for repl_obj in ith_sim_repls:
                            obj_stat_dict = repl_obj.get_stats_dict()

                            for st, stat_v in obj_stat_dict.items():
                                float_stat_v: float = 0.0

                                try:
                                    float_stat_v = float(stat_v)

                                except ValueError:
                                    # TODO: throw cannot convert to
                                    # float Exception
                                    pass

                                try:
                                    repls_all_stats_dict[st].append(
                                        float_stat_v)

                                except KeyError:
                                    repls_all_stats_dict[st] = [float_stat_v]

                        if "Origin age" in repls_all_stats_dict:
                            tree_repl_summary_df_dict[rv_name].at[j, "origin_age"] = \
                                "{:,.5f}".format(stat.mean(repls_all_stats_dict["Origin age"]))
                            sd = stat.stdev(repls_all_stats_dict["Origin age"])

                            if sd < 1E-10:
                                tree_repl_summary_df_dict[rv_name].at[j + 1, "origin_age"] = 0.0

                            else:
                                tree_repl_summary_df_dict[rv_name].at[j + 1, "origin_age"] = \
                                    "{:,.5f}".format(sd)

                        tree_repl_summary_df_dict[rv_name].at[j, "root_age"] = \
                            "{:,.5f}".format(stat.mean(repls_all_stats_dict["Root age"]))
                        sd = stat.stdev(repls_all_stats_dict["Root age"])

                        if sd < 1E-10:
                            tree_repl_summary_df_dict[rv_name].at[j + 1, "root_age"] = 0.0

                        else:
                            tree_repl_summary_df_dict[rv_name].at[j + 1, "root_age"] = \
                                "{:,.5f}".format(sd)

                        tree_repl_summary_df_dict[rv_name].at[j, "n_total"] = \
                            stat.mean(repls_all_stats_dict["Total taxon count"])
                        tree_repl_summary_df_dict[rv_name].at[j + 1, "n_total"] = \
                            "{:,.5f}".format(stat.stdev(repls_all_stats_dict["Total taxon count"]))
                        tree_repl_summary_df_dict[rv_name].at[j, "n_extant"] = \
                            stat.mean(repls_all_stats_dict["Extant taxon count"])
                        tree_repl_summary_df_dict[rv_name].at[j + 1, "n_extant"] = \
                            "{:,.5f}".format(stat.stdev(repls_all_stats_dict["Extant taxon count"]))
                        tree_repl_summary_df_dict[rv_name].at[j, "n_extinct"] = \
                            stat.mean(repls_all_stats_dict["Extinct taxon count"])
                        tree_repl_summary_df_dict[rv_name].at[j + 1, "n_extinct"] = \
                            "{:,.5f}".format(stat.stdev(repls_all_stats_dict["Extinct taxon count"]))
                        tree_repl_summary_df_dict[rv_name].at[j, "n_sa"] = \
                            stat.mean(repls_all_stats_dict["Direct ancestor count"])
                        tree_repl_summary_df_dict[rv_name].at[j + 1, "n_sa"] = \
                            "{:,.5f}".format(stat.stdev(repls_all_stats_dict["Direct ancestor count"]))

                        j += 2

                #########################################
                # Doing dictionaries with tsv and nexus #
                # strings for states                    #
                #########################################
                sample_idx = 0
                repl_idx = 0
                idx = 0
                for replicate_tree in node_val:
                    while (sample_idx < sample_size):
                        while (repl_idx < n_repl):
                            living_nd_states_tsv_fp = \
                                rv_name + "_sample" + str(sample_idx + 1) \
                                + "_repl" + str(repl_idx + 1) + ".tsv"
                            internal_nd_states_tsv_fp = \
                                rv_name + "_sample" + str(sample_idx + 1) \
                                + "_repl" + str(repl_idx + 1) + "_anc_states.tsv"
                            living_nd_states_str, int_node_states_str = \
                                node_val[idx].get_taxon_states_str()
                            tree_living_nd_states_str_dict[living_nd_states_tsv_fp] = \
                                living_nd_states_str
                            tree_internal_nd_states_str_dict[internal_nd_states_tsv_fp] = \
                                int_node_states_str

                            if write_nex_states:
                                states_nex_fp = rv_name + "_sample" \
                                    + str(sample_idx + 1) + "_repl" \
                                    + str(repl_idx + 1) + ".nex"
                                nexus_str = node_val[idx].get_taxon_states_str(nexus=True)
                                tree_living_nd_states_str_nexus_dict[
                                    states_nex_fp] = nexus_str

                            repl_idx += 1
                            idx += 1

                        repl_idx = 0
                        sample_idx += 1

                    # TODO
                    # while (sample_idx < n_samples):
                    #     while (repl_idx < n_repl):
                    #         nex_fp = node_name + "_sample" + str(sample_idx) + "_repl" + str(repl_idx) + ".nex"
                    #         nexus_str = node_val[idx].get_taxon_states_nexus_str()
                    #         tree_nexus_str_states_dict[nex_fp] = nexus_str
                    #         repl_idx += 1
                    #     sample_idx += 1

    # debugging
    # print("\n\nPrinting dataframe for scalars with fixed value:\n")
    # print(tabulate(scalar_constant_df, headers=scalar_constant_df.head(), tablefmt="pretty", showindex=False).lstrip())

    # for n_repl, scalar_df in scalar_value_df_dict.items():
    #     print("\n\nPrinting dataframe with scalars with " + str(n_repl) + " replicate(s):\n")
    #     print(tabulate(scalar_df, scalar_df.head(), tablefmt="pretty", showindex=False).lstrip())

    # print("\n\nPrinting dataframe with scalar replicate summaries:\n")
    # print(tabulate(scalar_repl_summary_df, scalar_repl_summary_df.head(), tablefmt="pretty", showindex=False).lstrip())

    # print("\n\nPrinting dataframe with complete tree values, with replicates if those exist:\n")
    # for tr_node_name, tr_df in tree_value_df_dict.items():
    #     print(tabulate(tr_df, tr_df.head(), tablefmt="plain", showindex=False).lstrip() + "\n")

    # print("\n\nPrinting dataframe with reconstructed tree values, with replicates if those exist:\n")
    # for tr_node_name, tr_df in tree_rec_value_df_dict.items():
    #     print(tabulate(tr_df, tr_df.head(), tablefmt="plain", showindex=False).lstrip() + "\n")

    # print("\n\nPrinting dataframe with tree summaries, with replicates if those exist:\n")
    # for tr_node_name, tr_df in tree_summary_df_dict.items():
    #     print(tabulate(tr_df, tr_df.head(), tablefmt="pretty", showindex=False).lstrip() + "\n")

    # print("\n\nPrinting dataframe with tree summaries summarized over replicates:\n")
    # for tr_node_name, tr_df in tree_repl_summary_df_dict.items():
    #     print(tr_node_name)
    #     print(tabulate(tr_df, tr_df.head(), tablefmt="pretty", showindex=False).lstrip())

    # print("\n\nPrinting tsv files with states from living nodes:\n")
    # for tsv_fp, states_str in tree_living_nd_states_str_dict.items():
    #     print(tsv_fp)
    #     print(states_str)

    # print("\n\nPrinting nexus files with states from all trees:\n")
    # for nex_fp, nexus_str in tree_nexus_str_states_dict.items():
    #     print(nex_fp)
    #     print(nexus_str)

    # print("\n\nPrinting tsv files with states from internal nodes in reconstructed tree:\n")
    # for tsv_fp, states_str in tree_internal_nd_states_str_dict.items():
    #     print(tsv_fp)
    #     print(states_str)

    ############################
    # Finally adding to return #
    ############################
    scalar_output_stash: \
        ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]] = []
    tree_output_stash: \
        ty.List[ty.Union[ty.Dict[str, pd.DataFrame], ty.Dict[str, str]]] = []

    scalar_output_stash.extend([scalar_constant_df,
                                scalar_value_df_dict,
                                scalar_repl_summary_df])

    tree_output_stash.extend([tree_value_df_dict,
                              tree_ann_value_df_dict,
                              tree_rec_value_df_dict,
                              tree_rec_ann_value_df_dict,
                              tree_summary_df_dict,
                              tree_repl_summary_df_dict,
                              tree_living_nd_states_str_dict,
                              tree_living_nd_states_str_nexus_dict,
                              tree_internal_nd_states_str_dict])

    return scalar_output_stash, tree_output_stash


def prep_data_filepaths_dfs(
    scalar_output_stash:
        ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]],
    tree_output_stash:
        ty.List[ty.Union[ty.Dict[str, pd.DataFrame], ty.Dict[str, str]]]) -> \
            ty.Tuple[ty.List[str], ty.List[ty.Union[pd.DataFrame, str]]]:
    """Prepare list of file paths and list of pandas DataFrames.
    
    Args:
        scalar_output_stash (list): List of either pandas dataframes,
            or dictionaries with number of replicates as keys, and
            pandas dataframes as values. These contain scalar simulated
            data.
        tree_output_stash (ty.List[ty.Dict[str, pd.DataFrame]]):
            List of dictionaries with tree node names as keys, and
            pandas' dataframes as values.
            These contains tree simulated data.

    Returns:
        (tuple): List of filepath strings and and list of pandas
            dataframes to be written to disk.
    """

    output_fp_list: ty.List[str] = []
    output_df_str_list: ty.List[pd.DataFrame] = []

    # scalar variables
    if not scalar_output_stash[0].empty:
        output_fp_list.append("scalar_constants.csv")
        output_df_str_list.append(scalar_output_stash[0])

    # one .csv file per plating replication
    if scalar_output_stash[1]:
        for n_repl, df in scalar_output_stash[1].items():
            output_fp_list.append("scalar_rvs_repl" + str(n_repl) + ".csv")
            output_df_str_list.append(df)

    if not scalar_output_stash[2].empty:
        output_fp_list.append("scalar_rvs_stats_summary.csv")
        output_df_str_list.append(scalar_output_stash[2])

    # tree variables
    if tree_output_stash[0]:
        for rv_name, df in tree_output_stash[0].items():
            output_fp_list.append(rv_name + "_complete.tsv")
            output_df_str_list.append(df)

    if tree_output_stash[1]:
        for rv_name, df in tree_output_stash[1].items():
            output_fp_list.append(rv_name + "_annotated_complete.tsv")
            output_df_str_list.append(df)

    if tree_output_stash[2]:
        for rv_name, df in tree_output_stash[2].items():
            output_fp_list.append(rv_name + "_reconstructed.tsv")
            output_df_str_list.append(df)

    if tree_output_stash[3]:
        for rv_name, df in tree_output_stash[3].items():
            output_fp_list.append(rv_name + "_annotated_reconstructed.tsv")
            output_df_str_list.append(df)

    if tree_output_stash[4]:
        for rv_name, df in tree_output_stash[4].items():
            output_fp_list.append(rv_name + "_stats.csv")
            output_df_str_list.append(df)

    if tree_output_stash[5]:
        for rv_name, df in tree_output_stash[5].items():
            output_fp_list.append(rv_name + "_stats_summary.csv")
            output_df_str_list.append(df)

    if tree_output_stash[6]:
        for state_tsv_fp, state_string in tree_output_stash[6].items():
            output_fp_list.append(state_tsv_fp)
            output_df_str_list.append(state_string)

    # contains states of living nodes
    # in nexus format and will be empty if user did not
    # as for it
    if tree_output_stash[7]:
        for state_nex_fp, state_nex_string in tree_output_stash[7].items():
            output_fp_list.append(state_nex_fp)
            output_df_str_list.append(state_nex_string)

    if tree_output_stash[8]:
        for state_tsv_fp, state_string in tree_output_stash[8].items():
            output_fp_list.append(state_tsv_fp)
            output_df_str_list.append(state_string)

    return output_fp_list, output_df_str_list


def dump_pgm_data(dir_string: str,
                  dag_obj: pgm.DirectedAcyclicGraph,
                  prefix: str = "",
                  write_nex_states: bool = False) -> None:
    """Write stochastic-node sampled values in specified directory

    Args:
        dir_string (str): Where to save the files to be written
        dag_obj (DirectedAcyclicGraph): DAG object whose sampled
            values we are extracting and writing to file.
        prefix (str): String to preceed file names

    Returns:
        None
    """

    sorted_node_dag_list: \
        ty.List[pgm.NodeDAG] = dag_obj.get_sorted_node_dag_list()

    # populating data stashes that will be dumped and their file names
    scalar_output_stash: \
        ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]]
    tree_output_stash: ty.List[ty.Dict[str, pd.DataFrame]]
    output_fp_list: ty.List[str]
    output_df_str_list: ty.List[pd.DataFrame]

    # still prefixless #
    scalar_output_stash, tree_output_stash = \
        prep_data_df(dag_obj, write_nex_states)

    output_fp_list, output_df_str_list = \
        prep_data_filepaths_dfs(scalar_output_stash, tree_output_stash)

    # sort out file path #
    if not dir_string.endswith("/"):
        dir_string += "/"

    output_file_path = dir_string

    if prefix:
        output_file_path += prefix + "_"

    # full file paths #
    data_df_full_fp_list = \
        [output_file_path + fn for fn in output_fp_list]

    # debugging
    # print("\n\n" + "\n".join(data_df_full_fp_list))
    # for df_or_str in output_df_str_list:
    #     if isinstance(df_or_str, pd.DataFrame):
    #         print(tabulate(df, df.head(), tablefmt="plain",
    #               showindex=False).lstrip())

    # write! #
    for idx, full_fp in enumerate(data_df_full_fp_list):
        f: str = ""
        to_print = output_df_str_list[idx]

        if full_fp.endswith("csv"):
            f = "csv"

        elif full_fp.endswith("tsv"):
            f = "tsv"

        elif full_fp.endswith("nex"):
            f = "nex"

        if isinstance(to_print, pd.DataFrame):
            with open(full_fp, "w") as data_out:
                write_data_df(data_out, to_print, format=f)

        elif isinstance(to_print, str):
            # ignore nexus file #
            if f == "nex" and not write_nex_states:
                continue

            with open(full_fp, "w") as data_out:
                write_str_list(data_out, [to_print])


def dump_serialized_pgm(
        file_name: str,
        dag_obj: pgm.DirectedAcyclicGraph,
        cmd_log_list: ty.List[str],
        prefix: str = "",
        to_folder: bool = False) -> None:
    """Write serialized DAG in specified directory.

    Args:
        file_name (str): Serialized file name
        dag_obj (DirectedAcyclicGraph): DAG object to be serialized
            and saved.
        prefix (str): String to preceed file name.
    """

    if file_name:
        if prefix:
            prefix += "_"

        # PySimpleGUI
        if to_folder:
            if not file_name.endswith("/"):
                file_name += "/"

            with open(file_name + prefix + ".pickle", "wb") as picklefile:
                pickle.dump((dag_obj, cmd_log_list), picklefile)

        # PySide6
        else:
            head, tail = os.path.split(file_name)
            if not head.endswith("/"):
                head += "/"

            with open(head + prefix + tail + ".pickle", "wb") as picklefile:
                pickle.dump((dag_obj, cmd_log_list), picklefile)


def get_write_inference_rev_scripts(
        all_sims_model_spec_list: ty.List[str],
        all_sims_mcmc_logging_spec_list: ty.List[str],
        dir_list: ty.List[str],
        prefix: str = "",
        write2file: bool = False) -> ty.List[str]:
    """Get and/or write full inference .Rev scripts

    Args:
        all_sims_model_spec_list (str): List of strings specifying
            just the model part of a .Rev script, one element per simulation
        all_sims_mcmc_logging_spec_list (str): List of strings
            specifying just the MCMC and logging part of a .Rev script,
            one element per simulation
        dir_list (str): List of three string specifying directories
            (inference root, scripts, results)
        prefix (str): String prefix to place before the name of files being
            written
        write2file (bool): If 'True', function writes to file.
            Defaults to 'False'.

    Returns:
        list of str(s): A list of full .Rev script string specifications,
            one per simulation
    """

    def iterate_specs_over_sims(n: int):
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
    inference_root_dir, scripts_dir, results_dir = dir_list

    prefix = str()
    if prefix: prefix += "_"
    else: prefix = ""

    # return
    all_sims_spec_strs_list = \
        [spec_str for spec_str in iterate_specs_over_sims(n_sim)]

    # in addition to returning rev script strings (one per simulation),
    # we also write to file
    if write2file:
        for dir_name in dir_list:
            if not os.path.isdir(dir_name):
                os.mkdir(dir_name)

        for i in range(n_sim):
            _ith_file_name = prefix + "sim_" + str(i) + ".rev"
            _ith_file_path = scripts_dir + _ith_file_name

            with open(_ith_file_path, "w") as rev_out:
                rev_out.write(all_sims_spec_strs_list[i])

    return all_sims_spec_strs_list


if __name__ == "__main__":
    # testing below

    # initializing model
    # model_fp = "examples/multiple_scalar_tree_plated.pj"
    model_fp = "examples/geosse.pj"
    dag_obj = cmdp.script2dag(model_fp, in_pj_file=True)

    # either: to see what's inside dataframes (uncomment debugging)
    prep_data_df(dag_obj)

    # or: to see file names and dataframes after organized (uncomment debugging)
    # dump_pgm_data("./", dag_obj, prefix="test")
