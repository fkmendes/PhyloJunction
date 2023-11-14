import sys
sys.path.extend(["../inference", "../../phylojunction"])
import typing as ty
import numpy as np

# pj imports
import phylojunction.inference.revbayes.rb_dn_parametric as rbpar
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec

def get_mcmc_logging_spec_list(a_node_pgm_name: str, moves_str: str, n_sim: int, mcmc_chain_length: int , prefix: str, results_dir: str) -> ty.List[str]:
    """Generate list of strings where each element (one per simulation) will configure RevBayes' MCMC (moves) and logging inside a .Rev script

    Args:
        a_node_pgm_name (str): One of the probabilistic graphical model nodes, required to set the model object
        moves_str (str): All moves to be carried out during MCMC, as a string containing new line characters
        n_sim (int): Number of simulations
        prefix (str): Prefix to preceed MCMC result .log file names
        results_dir (str): String specifying directory where MCMC result .log files shall be put by RevBayes

    Returns:
        (str): List of strings (one per simulation), each will configure RevBayes' MCMC (moves) and logging inside a .Rev script
    """
    mcmc_logging_spec_str: str = ""
    mcmc_logging_spec_str += "mymodel = model(" + a_node_pgm_name + ")\n\n"
    mcmc_logging_spec_str += "monitors[1] = mnScreen()\n"

    mcmc_logging_spec_list: ty.List[str] = list()
    for i in range(n_sim):
        mcmc_logging_spec_list.append(
            moves_str + "\n" + 
            mcmc_logging_spec_str + 
            "monitors[2] = mnModel(\"" + results_dir + prefix + "run" + str(i+1) + ".log\")\n\n" + \
            "mymcmc = mcmc(mymodel, moves, monitors)\n" + \
            "mymcmc.run(" + str(mcmc_chain_length) + ")\n\n" + \
            "quit()"
        )

    return mcmc_logging_spec_list


def pgm_obj_to_rev_inference_spec(pgm_obj: pgm.ProbabilisticGraphicalModel, inference_root_dir: str, mcmc_chain_length: int=1000, prefix: str="") -> ty.Tuple[ty.List[str], ty.List[str], ty.List[str]]:
    """Prepare string to put in .Rev script, from model specified by user in PhyloJunction

    Args:
        pgm_obj (ProbabilisticGraphicalModel): PGM object created through script or commands typed by user
        mcmc_chain_length (int): Number of MCMC iterations to carry out with RevBayes
        inference_root_dir (str): Directory where scripts and results should be put, and from where RevBayes should be called
        prefix (str): Prefix will be appended to directory and file names. Defaults to "".

    Returns:
        tuple: Tuple containing three lists, with strings to be printed to files (2) and directory names (1)
    """
    
    all_nodes_all_sims_spec_list: ty.List[ty.List[str]] = []
    all_nodes_moves_str: str = str()
    sorted_node_pgm_list: ty.List[pgm.NodePGM] = pgm_obj.get_sorted_node_pgm_list()
    node_name: str = str()
    n_sim = 0

    ########################
    # Going over PGM nodes #
    ########################
    node_count = 1
    for node_pgm in sorted_node_pgm_list:
        node_name = node_pgm.node_name
        is_clamped = node_pgm.is_clamped
        node_inference_spec_str: str
        node_inference_spec_list: ty.List[str] = [] # will contain all sims for this node 

        if isinstance(node_pgm, pgm.StochasticNodePGM):
            ################
            # Sampled node #
            ################
            if node_pgm.is_sampled:
                node_operator_weight = node_pgm.operator_weight
                dn_obj = node_pgm.sampling_dn

                if dn_obj.DN_NAME == "DnSSE": continue # TODO
                
                # getting distribution spec
                n_sim, n_repl, rev_str_list = rbpar.get_rev_str_from_dn_parametric_obj(dn_obj)

                # all simulations for this PGM node
                for ith_sim in range(n_sim):
                    ######################
                    # Observed data prep #
                    ######################
                    if is_clamped:
                        node_inference_spec_str = "truth_" + node_name + " <- "
                        start = ith_sim * n_repl
                        end = start + n_repl
                        
                        if type(node_pgm.value) in (list, np.ndarray):
                            if type(node_pgm.value[0]) in (float, int, str):
                                if n_repl > 1:
                                    node_inference_spec_str += "[" + ", ".join(str(v) for v in node_pgm.value[start:end]) + "]\n\n"
                                    # with replicates, we need a for loop
                                    node_inference_spec_str += "for (i in 1:" + str(n_repl) + ") {\n    " +  node_name + "[i] ~ " + rev_str_list[ith_sim] + "\n"
                                    node_inference_spec_str += "    " + node_name + "[i].clamp(" + "truth_" + node_name + "[i])"
                                    node_inference_spec_str += "\n}"
                                
                                # no replicates, single value
                                else:
                                    node_inference_spec_str += str(node_pgm.value[0]) + "\n\n"
                                    node_inference_spec_str += node_name + ".clamp(" + "truth_" + node_name + ")"

                            else:
                                print("parsing rb script, found list of objects in clamped node... ignoring for now")
                            # print(node_pgm.node_name + ": adding rev spec for sim = " + str(ith_sim))
                    
                        # with replicates, we need a for loop
                        # if n_repl > 1:
                        #     node_inference_spec_str += "for (i in 1:" + str(n_repl) + ") {\n    " +  node_name + "[i] ~ " + rev_str_list[ith_sim] + "\n"
                        #     node_inference_spec_str += "    " + node_name + "[i].clamp(" + "truth_" + node_name + "[i])"
                        #     node_inference_spec_str += "\n}"

                        # no replicates, single value
                        # else:
                        #     node_inference_spec_str += node_name + ".clamp(" + "truth_" + node_name + ")"

                        node_inference_spec_list.append(node_inference_spec_str)
                    
                    # not clamped, single distributed-as statement
                    else:
                        node_inference_spec_list.append(node_name + " ~ " + rev_str_list[ith_sim])

                # print("after all sims processed, I have")
                # print(node_inference_spec_str)                    
            
                # TODO
                # elif isinstance(node_pgm, DeterministicNodePGM):
                #     pj2rev_deterministic(node_pgm) # will determine which class inside, and convert it to rev

            ############################
            # Clamped stochastic node, #
            # e.g.,                    #
            # x <- 1.0,                #
            # x <- [1, 2]              #
            ############################
            else:
                n_vals = len(node_pgm.value)
                n_sim = pgm_obj.sample_size

                # in case user enters clamped node first
                # and immediately wants to see the rev script
                if n_sim == 0:
                    n_sim = n_vals

                # all simulations for this PGM node
                for ith_sim in range(n_sim):

                    #############
                    # Important #
                    #############
                    # in PJ, we assume that if a node has the same number of (clamped) assigned values
                    # as the the sample size (n_samples) of stochastic nodes (or if there are no stochastic
                    # nodes added to the PGM at all yet, see above code), then we take each of its
                    # values and split them among the sample-sized-rev scripts, e.g., if we do:
                    #
                    # clamped_node <- [1, 2] # with no other node in the PGM
                    # 
                    # Then we will have:
                    # clamped_node <- 1 # first .Rev script
                    # clamped_node <- 2 # second .Rev script
                    if n_vals == n_sim:
                        node_inference_spec_str = node_name + " <- " + node_pgm.value[ith_sim]

                    # if clamped nodes have a single value, this value goes into each and every
                    # rev scripts
                    elif n_vals == 1:
                        node_inference_spec_str = node_name + " <- " + node_pgm.value[0]
                    
                    # if clamped nodes have more values than sampled stochastic nodes,
                    # all these values go into each and every rev script
                    else:
                        raise ec.NodeInferenceDimensionalityError(node_name)
                    
                    # TODO: fix this
                    # node_inference_spec_str = node_name + " <- "

                    # # node_pgm.value is always a list
                    # if type(node_pgm.value[0]) in (float, int, str):
                    #     if len(node_pgm.value) == 1:
                    #         node_inference_spec_str += str(node_pgm.value[0])
                    #     else:
                    #         node_inference_spec_str += "[" + ", ".join(str(v) for v in node_pgm.value) + "]"

                    node_inference_spec_list.append(node_inference_spec_str)
                    

            # one item per pgm_node, each item with n_sim items
            all_nodes_all_sims_spec_list.append(node_inference_spec_list) # 1D: nodes, 2D: sims

            ##############################
            # MCMC move/operator weights #
            ############################## 
            if not is_clamped:
                all_nodes_moves_str += "\nmoves[" + str(node_count) + "] = mvSlide(" + node_name + ", delta=1, weight=" + str(node_operator_weight * 5) + ")\n"
                node_count += 1

    # print("before reorganizing, all_nodes_all_sims_spec_list is")
    # print(all_nodes_all_sims_spec_list)

    # finished looking at all nodes and all simulations, time to reorganize
    all_sims_model_spec_list = ["" for i in range(n_sim)] # 1D: simulations, 2D: nodes
    line_sep: str = "" # spacing between node specification lines
    
    try:
        # if there is more than one node, we will add empty lines between their specificatoin
        if len(all_nodes_all_sims_spec_list[0]) > 1:
            line_sep = "\n\n"
    # will be an empty list if deterministic for now
    except: pass
    
    # converting into, 1D: sims, 2D: nodes
    for j, jth_node_pgm_list_all_sims in enumerate(all_nodes_all_sims_spec_list):
        for i, ith_sim_this_node_str in enumerate(jth_node_pgm_list_all_sims):
            # if there are > 1 nodes and we are not looking at the last node
            if len(all_nodes_all_sims_spec_list) > 1 and j < (len(all_nodes_all_sims_spec_list)-1): 
                all_sims_model_spec_list[i] += ith_sim_this_node_str + line_sep
            
            # if it's the last one, we avoid end lines
            else:
                all_sims_model_spec_list[i] += ith_sim_this_node_str # to avoid endlines after last node

    ##############################
    # Sorting out inference dirs #
    ##############################
    script_subdir_name: str
    res_subdir_name: str
    if prefix:
        script_subdir_name = prefix + "_rev_scripts/"
        res_subdir_name = prefix + "_inference_results/"
        prefix += "_"

    else:
        script_subdir_name = "rev_scripts/"
        res_subdir_name = "inference_results/"
        prefix = ""

    if not inference_root_dir.endswith("/"):
        inference_root_dir += "/"

    scripts_dir = inference_root_dir + script_subdir_name
    results_dir = inference_root_dir + res_subdir_name 

    # the idea is for RB to be called from inference_root_dir
    all_sims_mcmc_logging_spec_list = get_mcmc_logging_spec_list(node_name, all_nodes_moves_str, n_sim, mcmc_chain_length, prefix, res_subdir_name)

    dir_list = [inference_root_dir, scripts_dir, results_dir]
    # print("\n".join(all_sims_model_spec_list[0]))

    # print(all_nodes_moves_str)

    # print(all_sims_mcmc_logging_spec_list[0])

    return all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list