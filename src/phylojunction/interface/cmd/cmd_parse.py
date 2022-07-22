import re
import io
import typing as ty
import matplotlib.pyplot as plt # type: ignore

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
# import user_interface.phylojunction_inference as pjinf
# import user_interface.phylojunction_io as pjio
import phylojunction.interface.cmd.cmd_parse_utils as cmdu
import phylojunction.interface.grammar.dn_grammar as dngrammar
import phylojunction.interface.grammar.det_fn_grammar as detgrammar
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

def script2pgm(script_file_path_or_model_spec: str, in_pj_file: bool=True) -> pgm.ProbabilisticGraphicalModel:
    """Go through lines in file object or specification string and populate probabilistic graphical model

    Args:
        script_file_path_or_model_spec (str): Path to script .pj file or full string specifying model directly
    """

    def _execute_spec_lines(all_lines_list: ty.List[str], pgm_obj: pgm.ProbabilisticGraphicalModel) -> None:
        for line in all_lines_list:
            # clear padding whitespaces
            line = line.lstrip().rstrip()

            _ = cmdline2pgm(pgm_obj, line)

    pgm_obj = pgm.ProbabilisticGraphicalModel()

    all_lines_list: ty.List[str] = []

    if in_pj_file:
        with open(script_file_path_or_model_spec, "r") as infile:
            all_lines_list = infile.readlines()

    else:
        all_lines_list = script_file_path_or_model_spec.split("\n")

    # side-effect: populates pgm_obj
    _execute_spec_lines(all_lines_list, pgm_obj)
        
    # debugging
    # seeing trees in script strings in main()
    # for node_pgm_name, node_pgm in pgm_obj.node_name_val_dict.items():    
    #     if node_pgm_name == "trs":
    #         fig = plt.figure() # note that pjgui uses matplotlib.figure.Figure (which is part of Matplotlib's OOP class library)
    #                                  # here, we instead use pyplot's figure, which is the Matlab-like state-machine API
    #         ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    #         ax.patch.set_alpha(0.0)
    #         ax.xaxis.set_ticks([])
    #         ax.yaxis.set_ticks([])
    #         ax.spines['left'].set_visible(False)
    #         ax.spines['bottom'].set_visible(False)
    #         ax.spines['right'].set_visible(False)
    #         ax.spines['top'].set_visible(False)
    #
    #         print(node_pgm.value[0].tree.as_string(schema="newick"))
    #
    #         pjtr.plot_ann_tree(
    #             node_pgm.value[0],
    #             ax,
    #             use_age=False,
    #             start_at_origin=True,
    #             sa_along_branches=True
    #         )
            
    #         plt.show()

    # debugging
    # as if we had clicked "See" in the inference tab
    # all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list = pjinf.pgm_obj_to_rev_inference_spec(pgm_obj, "./inference_test/", mcmc_chain_length=1000, prefix="test")
    # for i, ith_sim_model_spec in enumerate(all_sims_model_spec_list):
    #     print(ith_sim_model_spec)
    #     print(all_sims_mcmc_logging_spec_list[i])

    # pjio.output_inference_rev_scripts(all_sims_model_spec_list, all_sims_mcmc_logging_spec_list, dir_list, prefix="test")
    
    return pgm_obj

def cmdline2pgm(pgm_obj: pgm.ProbabilisticGraphicalModel, cmd_line: str) -> str:
    """Update ProbabilisticGraphicalModel object using user-entered command line in PJ's language

    Args:
        pgm_obj (ProbabilisticGraphicalModel): Object that will hold all stochastic nodes (random variables) created by user via commands
        cmd_line (str): Command line string provided by user through static script or via GUI

    Returns:
        str: Command string if there was nothing wrong with it
    """

    # skip lines without useful commands
    if cmd_line.startswith("#") or cmd_line.startswith("\n") or not cmd_line: return ""

    ##################################
    # Checking command line is valid #
    ##################################
    if re.match(cmdu.cannot_start_with_this_regex, cmd_line):
        raise ec.ScriptSyntaxError(cmd_line, "A special character or digit was detected as the first character in this line, but this is not allowed. Exiting...")

    if re.search(cmdu.cannot_end_with_this_regex, cmd_line):
        raise ec.ScriptSyntaxError(cmd_line, "Cannot end command with special characters, apart from \']\' and \')\'. Exiting...")

    assignment_call_count = len(re.findall(cmdu.assign_regex, cmd_line))
    sampled_as_call_count = len(re.findall(cmdu.sampled_as_regex, cmd_line))
    deterministic_call_count = len(re.findall(cmdu.deterministic_regex, cmd_line))

    if (assignment_call_count + sampled_as_call_count + deterministic_call_count) > 1:
        raise ec.ScriptSyntaxError(cmd_line, "Only _one_ of the following can be done:\n    (1) [variable] <- [value]\n    (2) [variable] ~ [sampling distribution]" + \
            "\n    (3) [variable] := [deterministic function]\nExiting...")

    elif (assignment_call_count + sampled_as_call_count + deterministic_call_count) == 0:
        print("command line that caused the bombing = " + cmd_line)
        raise ec.ScriptSyntaxError(cmd_line, "One of the following must be done:\n    (1) [variable] <- [value]\n    (2) [variable] ~ [sampling distribution]" + \
            "\n    (3) [variable] := [deterministic function]\nExiting...")

    ##########################
    # Executing command line #
    ##########################
    
    # (1) variable assignment
    elif assignment_call_count == 1:
        try:
            rv_name, operator_str, rv_spec = re.split(cmdu.assign_regex, cmd_line)
            parse_variable_assignment(pgm_obj, rv_name, rv_spec, cmd_line)

        except:
            raise ec.ScriptSyntaxError(cmd_line, "Something went wrong during variable assignment. Variable name and value could not be tokenized. Exiting...")


    # (2) sampling distribution assignment
    elif sampled_as_call_count == 1:
        try:
            rv_name, operator_str, rv_dn_spec = re.split(cmdu.sampled_as_regex, cmd_line)
            parse_samp_dn_assignment(pgm_obj, rv_name, rv_dn_spec, cmd_line)

        except:
            raise ec.ScriptSyntaxError(cmd_line, "Something went wrong during sampling distribution assignment. Variable name and sampling distribution specification could not be tokenized. Exiting...")


    # (3) deterministic function assignment
    elif deterministic_call_count == 1:
        try:
            det_nd_name, operator_str, det_fn_spec = re.split(cmdu.deterministic_regex, cmd_line)
            parse_deterministic_function_assignment(pgm_obj, det_nd_name, det_fn_spec, cmd_line)
        except:
            raise ec.ScriptSyntaxError(cmd_line, "Something went wrong during deterministic function assignment. Variable name and deterministic function specification could not be tokenized. Exiting...")

    # cmd_line is valid, returning (used by GUI, when displaying command in History window)!
    return cmd_line


#####################
# Parsing functions #
#####################

# function called by (1. variable assignment)
def parse_variable_assignment(pgm_obj: pgm.ProbabilisticGraphicalModel, stoch_node_name: str, stoch_node_spec: str, cmd_line: str) -> None:
    """Create StochasticNodePGM instance from command string with "<-" operator, then add it to ProbabilisticGraphiclModel instance

    This node will not be sampled and is in practice a constant node.

    Args:
        pgm_obj (ProbabilisticGraphicalModel): Object that will hold all nodes created by user via commands
        stoch_node_name (str): Name of stochastic node being created
        stoch_node_spec (str): Stochastic node specification string (whatever is right of "<-" operator in PJ command)
        cmd_line (str): Command line string provided by user through static script or via GUI

    Returns:
        None
    """

    def create_add_stoch_node_pgm(a_stoch_node_name, sample_size, a_val_obj_list):
        stoch_node = pgm.StochasticNodePGM(a_stoch_node_name, sample_size=sample_size, value=a_val_obj_list, clamped=True)
        pgm_obj.add_node(stoch_node)

    #############
    # IMPORTANT #
    #############
    values_list = list() # arguments of dn/det parameters will come out as lists

    # if argument is a list
    if re.match(cmdu.vector_value_regex, stoch_node_spec):
        values_list = cmdu.parse_val_vector(stoch_node_spec)
    
    # if scalar variable
    else:
        if len(re.findall("\.", stoch_node_spec)) > 1:
            raise ec.ScriptSyntaxError(cmd_line, "Something went wrong during variable assignment. It looks like there were two dots \".\" when only one is allowed. Exiting...")
        values_list.append(stoch_node_spec)

    val_obj_list = cmdu.val_or_obj(pgm_obj, values_list)

    # if rv_spec is scalar
    # if re.match(int_or_float_regex, rv_spec) and len(re.findall(int_or_float_regex, rv_spec)) > 1:
    #     raise ScriptSyntaxError(rv_spec, "Something went wrong during variable assignment. If not a vector, the assigned value must consist of a single integer or float. Exiting...")

    # # ok!
    # else:
    #     scalar_values = rv_spec

    # # if rv_spec is vector
    # if re.match(vector_value_regex, rv_spec):
    #     scalar_values = parse_val_vector(rv_spec)

    # if not re.match(int_or_float_regex, rv_spec) and not re.match(vector_value_regex, rv_spec):
    #     raise ScriptSyntaxError(rv_spec, "Something went wrong during variable assignment. The value assigned to a variable must be a numeric scalar or a vector (not a character). Exiting...")

    n_samples = len(val_obj_list)
    create_add_stoch_node_pgm(stoch_node_name, n_samples, val_obj_list)


# function called by (2. sampling distribution assignment)
def parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_dn_spec, cmd_line):
    """Create StochasticNodePGM instance from command string with "~" operator, then add it to ProbabilisticGraphiclModel instance

    This node will be sampled from a distribution and is a random variable.

    Args:
        pgm_obj (ProbabilisticGraphicalModel): Object that will hold all nodes created by user via commands
        stoch_node_name (str): Name of stochastic node being created
        stoch_node_dn_spec (str): Specification string for distribution from which stochastic node will be sampled (whatever is right of "~" operator in PJ command)
        cmd_line (str): Command line string provided by user through static script or via GUI

    Returns:
        None
    """

    def create_add_rv_pgm(a_stoch_node_name: str, sample_size: int, replicate_size: int, a_dn_obj: pgm.DistributionPGM, parent_pgm_nodes: ty.List[pgm.NodePGM], clamped: bool):
        # set dn inside rv, then call .sample
        stoch_node_pgm = pgm.StochasticNodePGM(a_stoch_node_name, sample_size=sample_size, replicate_size=replicate_size, sampled_from=a_dn_obj, parent_nodes=parent_pgm_nodes, clamped=clamped)
        pgm_obj.add_node(stoch_node_pgm)

    if len(re.search(cmdu.sampling_dn_spec_regex, stoch_node_dn_spec).groups()) != 2:
        raise ec.ScriptSyntaxError(stoch_node_dn_spec, "Something went wrong during sampling distribution specification. Exiting...")

    # ok!
    else:
        dn_name, dn_spec = re.search(cmdu.sampling_dn_spec_regex, stoch_node_dn_spec).groups()

        if not dn_name in dngrammar.PJDnGrammar.dn_grammar_dict:
            raise ec.ScriptSyntaxError(stoch_node_dn_spec, "Something went wrong during sampling distribution specification. Distribution name not recognized. Exiting...")

        # parses, e.g., "par1=arg1, par2=arg2" into { par1:arg1, par2:arg2 }
        spec_dict, parent_pgm_nodes = cmdu.parse_spec(pgm_obj, dn_spec, cmd_line)
        dn_obj = dngrammar.PJDnGrammar.create_dn_obj(dn_name, spec_dict)

        # deal with clamping
        clamped = False
        try:
            clamp = spec_dict["clamp"][0]

            if clamp in ("\"true\"", "\"T\"", "\"True\""):
                clamped = True
            if clamp in ("\"false\"", "\"F\"", "\"False\""):
                clamped = False

        except: pass

        sample_size = 1
        try:
            sample_size = int(spec_dict["n"][0])
        except:
            pass

        # deal with replicates in a single simulation
        repl_size = 1
        try:
            repl_size = int(spec_dict["nr"][0])
        except: pass # user did not provide replicate size

        create_add_rv_pgm(stoch_node_name, sample_size, repl_size, dn_obj, parent_pgm_nodes, clamped)


# function called by (3. deterministic function assignment)
def parse_deterministic_function_assignment(pgm_obj: pgm.ProbabilisticGraphicalModel, det_nd_name: str, det_fn_spec: str, cmd_line: str) -> None:
    """Create DeterministicNodePGM instance from command string with ":=" operator, then add it to ProbabilisticGraphiclModel instance

    This node is not sampled (not a random variable) and is deterministically initialized from a deterministic function.

    Args:
        pgm_obj (ProbabilisticGraphicalModel): Object that will hold all nodes created by user via commands
        det_nd_name (str): Name of deterministic node being created
        det_fn_spec (str): Specification string for deterministic function that will specify the deterministic node (whatever is right of ":=" operator in PJ command)
        cmd_line (str): Command line string provided by user through static script or via GUI

    Returns:
        None
    """

    def create_add_det_nd_pgm(det_nd_name: str, det_obj: ty.Any, parent_pgm_nodes: ty.List[pgm.NodePGM]):
        det_nd_pgm = pgm.DeterministicNodePGM(det_nd_name, value=det_obj, parent_nodes=parent_pgm_nodes)
        pgm_obj.add_node(det_nd_pgm)
        # deterministic node is of class DeterministicNodePGM, which
        # derives NodePGM -- we do not need to initialize NodePGM
        # as when a new StochasticNodePGM is created (see above in
        # create_add_rv_pgm())
        #
        # this check is just to make sure we are adding a class
        # deriving from NodePGM
        # if isinstance(det_obj, NodePGM):
        #     det_obj.node_pgm_name = det_nd_name
        #     det_obj.parent_nodes = parent_pgm_nodes
        #     det_nd_pgm = det_obj
        #     pgm_obj.add_node(det_nd_pgm)

    sampling_dn_matches = re.search(cmdu.sampling_dn_spec_regex, det_fn_spec)
    if not sampling_dn_matches is None and len(sampling_dn_matches.groups()) != 2:
            raise ec.ScriptSyntaxError(det_fn_spec, "Something went wrong during deterministic function specification. Exiting...")

    # ok!
    else:
        det_fn_matches = re.search(cmdu.sampling_dn_spec_regex, det_fn_spec)        
        if not det_fn_matches is None:
            det_fn_name, dn_spec = det_fn_matches.groups()
        
            if not det_fn_name in detgrammar.PJDetFnGrammar.det_fn_grammar_dict:
                raise ec.ScriptSyntaxError(det_fn_spec, "Something went wrong during sampling distribution specification. Distribution name not recognized. Exiting...")

            # parses, e.g., "par1=arg1, par2=arg2" into { par1:arg1, par2:arg2 }
            spec_dict, parent_pgm_nodes = cmdu.parse_spec(pgm_obj, dn_spec, cmd_line)

            det_fn_obj = detgrammar.PJDetFnGrammar.create_det_fn_obj(det_fn_name, spec_dict)

            create_add_det_nd_pgm(det_nd_name, det_fn_obj, parent_pgm_nodes)


if __name__ == "__main__":

    script_str1 = "a <- 1.0"
    script_str2 = "a <- [1, 2, 3]"
    script_str3 = "a <- 1\nb <- 2\nc <- [a, b, 3]"
    script_str4 = "a ~ lognormal(n=10, mean=10.0, sd=1.0)"
    script_str5 = "l0 <- 1.0\na := sse_rate(name=\"lambda\", value=l0, event=\"b_speciation\", states=[2, 0, 1])"
    script_str6 = "l0rate := sse_rate(name=\"lambda\", value=1.0, event=\"w_speciation\")"
    script_str6_2 = "l0rate := sse_rate(name=\"lambda\", value=[1.0, 1.1], event=\"w_speciation\")"
    script_str7 = script_str6 + "\nmeh := sse_wrap(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)"
    script_str8 = script_str7 + "\ntr ~ discrete_sse(n=2, meh=meh, start_state=[0,0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str9 = "n <- 1000\nrv ~ unif(n=n, min=0.0, max=1.0)"
    script_str10 = \
        "l0rate := sse_rate(name=\"lambda0\", value=1.0, states=[0,0,0], event=\"w_speciation\")\n" + \
        "l1rate := sse_rate(name=\"lambda1\", value=1.0, states=[1,1,1], event=\"w_speciation\")\n" + \
        "m0rate := sse_rate(name=\"mu0\", value=0.8, states=[0,0,0], event=\"extinction\")\n" + \
        "m1rate := sse_rate(name=\"mu1\", value=0.8, states=[1,1,1], event=\"extinction\")\n" + \
        "q01rate := sse_rate(name=\"q01\", value=0.6, states=[0,1], event=\"transition\")\n" + \
        "q10rate := sse_rate(name=\"q10\", value=0.6, states=[1,0], event=\"transition\")\n" + \
        "meh := sse_wrap(flat_rate_mat=[l0rate, l1rate, m0rate, m1rate, q01rate, q10rate], n_states=2, n_epochs=1)\n" + \
        "tr ~ discrete_sse(n=1, meh=meh, start_state=[0], stop=\"age\", stop_value=2.5, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str11 = "m <- [1.0, 1.1]\ns <- [1.0, 0.9]\na ~ lognormal(n=2, mean=m, sd=s, log_space=\"false\")"
    script_str12 = "a ~ lognormal(n=3, mean=[0.0, 0.01, 0.011], sd=[1.0], log_space=\"true\")"
    script_str13 = "a ~ lognormal(mean=0.0, sd=1.0)"
    script_str14 = "a ~ lognormal(n=3, mean=[0.0, 0.01, 0.011], sd=1.0)"
    script_str15 = "a ~ normal(n=3, mean=[0.0, 0.01, 0.011], sd=[1.0])"
    script_str16 = "a ~ exponential(n=3, rate=1.0)"
    script_str17 = "a ~ exponential(rate=1.0, rate_parameterization=\"true\")"
    # TODO: must make str18 be allowed together with str20, by multiplying the single value drawn from uniform
    script_str18 = "l0 ~ unif(n=1, min=0.8, max=0.8)\nl0rate := sse_rate(name=\"lambda0\", value=l0, states=[0,0,0], event=\"w_speciation\")"
    script_str19 = "l0 ~ unif(n=3, min=0.8, max=0.8)\nl0rate := sse_rate(name=\"lambda0\", value=l0, states=[0,0,0], event=\"w_speciation\")"
    script_str20 = script_str18 + "\nmeh := sse_wrap(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)" + \
        "\ntr ~ discrete_sse(n=3, meh=meh, start_state=[0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str21 = script_str19 + "\nmeh := sse_wrap(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)" + \
                    "\ntr ~ discrete_sse(n=3, meh=meh, start_state=[0,0,0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str22 = "rv1 ~ lognormal(n=2, mean=0.0, sd=0.25)\nrv2 ~ lognormal(n=2, mean=-1.5, sd=0.25)\nrv3 ~ normal(n=2, nr=10, mean=rv1, sd=rv2, clamp=\"true\")"
    script_str23 = "l0ratet1 := sse_rate(name=\"lambdat1\", value=1.0, event=\"w_speciation\")\n" + \
        "l0ratet2 := sse_rate(name=\"lambdat2\", value=0.25, event=\"w_speciation\")\n" + \
        "l0ratet3 := sse_rate(name=\"lambdat3\", value=3.0, event=\"w_speciation\")\n" + \
        "l0ratet4 := sse_rate(name=\"lambdat4\", value=0.4, event=\"w_speciation\")\n" + \
        "meh := sse_wrap(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3, l0ratet4], n_states=1, n_epochs=4, seed_age=3.0, epoch_age_ends=[2.2, 1.2, 0.7])"
    script_str24 = script_str23 + "\ntr ~ discrete_sse(n=2, meh=meh, start_state=[0], stop=\"age\", stop_value=3.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str25 = "a <- [1,2,3]\nb <- a"
    script_str26 = "a <- [1,2,3]\nb <- [a, 4, 5, 6]"
    script_str27 = "rv9 ~ exponential(n=5, rate=[0.8, 0.9, 1.0, 1.1, 1.2], clamp=\"true\")"

    script_str28 = \
        "birth_rate := sse_rate(name=\"lambda\", value=1.0, states=[0,0,0], event=\"w_speciation\")\n" + \
        "death_rate := sse_rate(name=\"mu\", value=0.5, states=[0], event=\"extinction\")\n" + \
        "fossil_rate := sse_rate(name=\"psi\", value=0.8, states=[0], event=\"anc_sampling\")\n" + \
        "meh := sse_wrap(flat_rate_mat=[birth_rate,death_rate,fossil_rate], n_states=1, n_epochs=1)\n" + \
        "trs ~ discrete_sse(n=5, meh=meh, start_state=[0], stop=\"age\", stop_value=3.0, origin=\"true\", cond_spn=\"true\", cond_surv=\"true\")"

    script_str29 = \
        "w_birth_rate0 <- 1.0\n" + \
        "death_rate0 <- 0.2\n" + \
        "trans_rate_01 <- 0.5\n" + \
        "w_birth_rate1 <- 0.6\n" + \
        "death_rate1 <- 0.4\n" + \
        "trans_rate_10 <- 0.1\n" + \
        "b_birth_rate_201 <- 2.0\n" + \
        "b_birth_rate_202 <- 1.0\n" + \
        "b_birth_rate_212 <- 1.0\n" + \
        "trans_rate_20 <- 0.6\n" + \
        "trans_rate_21 <- 0.6\n" + \
        "det_w_birth_rate0 := sse_rate(name=\"birth_rate0\", value=w_birth_rate0, states=[0,0,0], event=\"w_speciation\")\n" + \
        "det_death_rate0 := sse_rate(name=\"death_rate0\", value=death_rate0, states=[0], event=\"extinction\")\n" + \
        "det_trans_rate_01 := sse_rate(name=\"q01\", value=trans_rate_01, states=[0,1], event=\"transition\")\n" + \
        "det_w_birth_rate1 := sse_rate(name=\"birth_rate0\", value=w_birth_rate0, states=[1,1,1], event=\"w_speciation\")\n" + \
        "det_death_rate1 := sse_rate(name=\"death_rate0\", value=death_rate0, states=[1], event=\"extinction\")\n" + \
        "det_trans_rate_10 := sse_rate(name=\"q10\", value=trans_rate_01, states=[1,0], event=\"transition\")\n" + \
        "det_b_birth_rate201 := sse_rate(name=\"birth_rate_20\", value=b_birth_rate_201, states=[2,0,1], event=\"bw_speciation\")\n" + \
        "det_b_birth_rate202 := sse_rate(name=\"birth_rate_201\", value=b_birth_rate_202, states=[2,0,2], event=\"bw_speciation\")\n" + \
        "det_b_birth_rate212 := sse_rate(name=\"birth_rate_201\", value=b_birth_rate_212, states=[2,1,2], event=\"bw_speciation\")\n" + \
        "det_trans_rate_20 := sse_rate(name=\"q20\", value=trans_rate_20, states=[2,0], event=\"transition\")\n" + \
        "det_trans_rate_21 := sse_rate(name=\"q21\", value=trans_rate_21, states=[2,1], event=\"transition\")\n" + \
        "meh := sse_wrap(flat_rate_mat=[det_w_birth_rate0, det_death_rate0, det_trans_rate_01, det_w_birth_rate1, det_death_rate1, det_trans_rate_10, det_b_birth_rate201, det_b_birth_rate202, det_b_birth_rate212, det_trans_rate_20, det_trans_rate_21], n_states=3, n_epochs=1)"

    # for copying and pasting in GUI:
    #
    # vectorized atomic rate parameter
    # l0rate := sse_rate(name="lambda", value=[1.0, 1.1], event="w_speciation")
    # meh := sse_wrap(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)

    # yule
    # l0 <- [0.9, 1.0, 1.1]
    # l0 ~ unif(n=3, min=0.9, max=1.1)
    # l0rate := sse_rate(name="lambda", value=l0, event="w_speciation")
    # meh := sse_wrap(flat_rate_mat=[l0rate], n_states=1, n_epochs=1)
    # trs ~ discrete_sse(n=3, meh=meh, start_state=[0], stop="age", stop_value=2.0, origin="true", cond_spn="false", cond_surv="true")
    # trs ~ discrete_sse(n=3, nr=2, meh=meh, start_state=[0], stop="age", stop_value=2.0, origin="true", cond_spn="false", cond_surv="true")

    # bisse
    # l0rate := sse_rate(name="lambda0", value=1.0, states=[0,0,0], event="w_speciation")
    # l1rate := sse_rate(name="lambda1", value=1.0, states=[1,1,1], event="w_speciation")
    # m0rate := sse_rate(name="mu0", value=0.8, states=[0,0,0], event="extinction")
    # m1rate := sse_rate(name="mu1", value=0.8, states=[1,1,1], event="extinction")
    # q01rate := sse_rate(name="q01", value=0.6, states=[0,1], event="transition")
    # q10rate := sse_rate(name="q10", value=0.6, states=[1,0], event="transition")
    # meh := sse_wrap(flat_rate_mat=[l0rate, l1rate, m0rate, m1rate, q01rate, q10rate], n_states=2, n_epochs=1)
    # trs ~ discrete_sse(n=5, meh=meh, start_state=[0], stop="age", stop_value=2.5, origin="true", cond_spn="false", cond_surv="true")

    # time-het yule
    # l0ratet1 := sse_rate(name="lambdat1", value=1.0, event="w_speciation")
    # l0ratet2 := sse_rate(name="lambdat2", value=0.25, event="w_speciation")
    # l0ratet3 := sse_rate(name="lambdat3", value=3.0, event="w_speciation")
    # l0ratet4 := sse_rate(name="lambdat4", value=0.4, event="w_speciation")
    # meh := sse_wrap(flat_rate_mat=[l0ratet1, l0ratet2, l0ratet3, l0ratet4], n_states=1, n_epochs=4, seed_age=3.0, epoch_age_ends=[2.2, 1.2, 0.7])
    # tr ~ discrete_sse(n=10, meh=meh, start_state=[0], stop="age", stop_value=3.0, origin="true", cond_spn="false", cond_surv="true")

    # simple hierarchical
    # rv1 ~ lognormal(n=2, mean=0.0, sd=0.25)
    # rv2 ~ lognormal(n=2, mean=-1.5, sd=0.25)
    # rv3 ~ normal(n=2, nr=10, mean=rv1, sd=rv2, clamp="true")

    # should throw exceptions
    # script_str20 = "1.0 + 1.0"
    # script_str21 = "~ lognormal(n=10, mean=0.0, sd=1.0)"
    # script_str22 = "a <-- 1.0"
    # script_str23 = "a <- 1.0 <- b"
    # script_str24 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ~"
    # script_str25 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) :="
    # script_str26 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ("
    # script_str27 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) ["
    # script_str28 = "a ~ lognormal(n=10, mean=0.0, sd=1.0) {"
    # script_str29 = "a <- 1.0.0"
    # script_str30 = "a <- 2.0"
    script_str31 = "a ~ lognormal(n=10, mean=[0.0, 0.01], sd=[1.0, 0.9])"
    script_str32 = "a ~ lognormal(n=2, mean=[0.0, 0.01], sd=[1.0, 0.9, 0.8])"
    script_str33 = "l0 ~ unif(n=3, nr=2, min=0.9, max=1.1)\nl0rate := sse_rate(name=\"lambda\", value=l0, event=\"w_speciation\")"
    script_str34 = script_str23 + "\ntr ~ discrete_sse(n=1, meh=meh, start_state=[0], stop=\"size\", stop_value=3.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"
    script_str35 = script_str23 + "\ntr ~ discrete_sse(n=1, meh=meh, start_state=[0], stop=\"age\", stop_value=2.0, origin=\"true\", cond_spn=\"false\", cond_surv=\"true\")"

    file_handle_exception = io.StringIO(script_str35)

    # script2pgm(file_handle_exception)

    # script2pgm(file_handle1)
    # script2pgm(file_handle2)
    # script2pgm(file_handle3)
    # script2pgm(file_handle4)
    # script2pgm(file_handle5)
    # script2pgm(file_handle6)
    # script2pgm(file_handle6_2)
    # script2pgm(file_handle7)
    # script2pgm(file_handle8)
    # script2pgm(file_handle9)
    # script2pgm(file_handle10)
    # script2pgm(file_handle11)
    # script2pgm(file_handle12)
    # script2pgm(file_handle13)
    # script2pgm(file_handle14)
    # script2pgm(file_handle15)
    # script2pgm(file_handle16)
    # script2pgm(file_handle17)
    # script2pgm(file_handle18)
    # script2pgm(file_handle19)
    # script2pgm(file_handle20)
    # script2pgm(file_handle21)
    # script2pgm(file_handle22)
    # script2pgm(file_handle23)
    # script2pgm(file_handle24)
    # script2pgm(file_handle25)
    # script2pgm(file_handle26) # will cause error if you generate rev script
    # script2pgm(file_handle27)
    script2pgm(script_str29, in_pj_file=False)