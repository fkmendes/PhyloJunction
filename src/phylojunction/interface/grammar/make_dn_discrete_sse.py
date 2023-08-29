from os import abort
import typing as ty

# pj imports
import phylojunction.distribution.dn_discrete_sse as dnsse
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def make_discrete_SSE_dn(
    dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
        -> pgm.DistributionPGM:

    ############################
    # Validating dn_param_dict #
    ############################
    if not dn_param_dict:
        raise ec.MissingSpecificationError(message="Cannot initialize discrete SSE distribution without specifications. Exiting...")

    #############################
    # IMPORTANT: Default values #
    #############################
    # default values for args that are not mandatory
    # all remaining args must be specified
    n_samples: int = 1
    n_repl: int = 1
    stash: sseobj.SSEStash = sseobj.SSEStash(
        sseobj.MacroevolEventHandler(
            sseobj.DiscreteStateDependentParameterManager([[]], 1)
        )
    )
    # event_handler: \
    #     sseobj.MacroevolEventHandler \
    #         = sseobj.MacroevolEventHandler(
    #             sseobj.DiscreteStateDependentParameterManager([[]], 1))
    # state_dep_prob_handler: \
    #     sseobj.DiscreteStateDependentProbabilityHandler \
    #         = sseobj.DiscreteStateDependentProbabilityHandler(
    #             sseobj.DiscreteStateDependentParameterManager([[]], 1))
    stop_str: str
    stop_values_list: ty.List[float] = []
    origin: bool = True
    cond_spn: bool = False
    cond_surv: bool = True

    # rejection sampling based on number
    # of observed taxa in reconstructed tree
    min_rec_taxa = 0
    max_rec_taxa = int(1e12)
    abort_at_obs = int(1e12)  # for trees out of control

    # if set to True by user, rejection sampling happens;
    # if False, reconstructed trees are printed whatever
    # they may be, upon writing (i.e., no rejection sampling)
    cond_obs_both_sides: bool = False

    start_states_list: ty.List[int] = []
    eps: float = 1e-12
    runtime_limit: int = 5  # 5 minutes

    # input validation (only things that are not already done by DnSSE)
    for arg, val in dn_param_dict.items():

        # val is a list!
        if val:

            # if element in val is string: remains unchanged
            #
            # if StochasticNodePGM: we get its string-fied value
            extracted_val = pgm.extract_value_from_nodepgm(val)

            ############################
            # Non-vectorized arguments #
            ############################

            # ... thus using only the first value!
            first_val: ty.Union[str, pgm.NodePGM]

            if len(extracted_val) >= 1:
                first_val = extracted_val[0]

            # if DeterministicNodePGM is in val
            # e.g., val = [pgm.DeterministicNodePGM]
            else:
                first_val = val[0]

            if arg == "stash" and isinstance(first_val, pgm.NodePGM):
                nodepgm_val = first_val.value

                if isinstance(nodepgm_val, sseobj.SSEStash):
                    stash = nodepgm_val
                # there should be only one event handler always, but it will be in list
                    # event_handler = nodepgm_val.get_meh()

                    # SSEStash will return None if prob_handler
                    # wasn't created by user through script
                    # prob_handler = nodepgm_val.get_prob_handler()

            elif arg in ("n", "nr", "runtime_limit", "min_rec_taxa",
                         "max_rec_taxa", "abort_at_obs"):
                try:
                    if isinstance(first_val, str):
                        # no vectorization allowed here
                        int_val = int(first_val)

                except:
                    raise ec.FunctionArgError(
                        arg, ("Was expecting an integer. Distribution "
                              "discrete_sse() could not be initialized."))

                # if user specified n or runtime_limit, we use it, otherwise defaults are used
                if arg == "n": n_samples = int_val
                if arg == "runtime_limit": runtime_limit = int_val
                if arg == "nr": n_repl = int_val
                if arg == "min_rec_taxa": min_rec_taxa = int_val
                if arg == "max_rec_taxa": max_rec_taxa = int_val
                if arg == "abort_at_obs": abort_at_obs = int_val

            elif arg == "stop":
                if isinstance(first_val, str):
                    stop_str = first_val.replace("\"", "")

            # defaults: True, False, True
            elif arg in ("origin", "cond_spn", "cond_surv", "cond_obs_both_sides"):
                if first_val != "\"true\"" and first_val != "\"false\"":
                    raise ec.FunctionArgError(arg, "Was expecting either \"true\" or \"false\". Distribution discrete_sse() could not be initialized.")

                else:
                    if isinstance(first_val, str):
                        parsed_val = first_val.replace("\"", "")

                    if arg == "origin" and parsed_val == "true": origin = True
                    elif arg == "origin" and parsed_val == "false": origin = False
                    elif arg == "cond_spn" and parsed_val == "true": cond_spn = True
                    elif arg == "cond_spn" and parsed_val == "false": cond_spn = False
                    elif arg == "cond_surv" and parsed_val == "true": cond_surv = True
                    elif arg == "cond_surv" and parsed_val == "false": cond_surv = False
                    elif arg == "cond_obs_both_sides" and parsed_val == "true": cond_obs_both_sides = True
                    elif arg == "cond_obs_both_sides" and parsed_val == "false": cond_obs_both_sides = False

            # if user specified epsilon, we use it, otherwise default is used
            elif arg == "eps":
                try:
                    if isinstance(first_val, str):
                        # no vectorization allowed here
                        float_val = float(first_val)

                except:
                    raise ec.FunctionArgError(arg, "Was expecting a double. Distribution discrete_sse() could not be initialized.")

                eps = float_val

            #####################
            # Vectorized values #
            #####################
            elif arg == "stop_value":
                _is_negative = False

                try:
                    stop_values_list = [float(v) for v in extracted_val if isinstance(v, str)]

                    for sv in stop_values_list:
                        if sv < 0.0:
                            _is_negative = True

                except:
                    if _is_negative:
                        raise ec.FunctionArgError(arg, "Stop value cannot be negative. Distribution discrete_sse() could not be initialized.")

                    raise ec.FunctionArgError(arg, "Was expecting number for stopping value(s). Distribution discrete_sse() could not be initialized.")

            elif arg == "start_state":
                try:
                    start_states_list = [int(v) for v in extracted_val if isinstance(v, str)]

                except:
                    raise ec.FunctionArgError(
                        arg,
                        ("Was expecting integers for starting states. "
                         "Distribution discrete_sse() could not be "
                         "initialized."))

    # making sure essential parameters of distribution have been specified
    if not any(stash.get_meh().state_dep_rate_manager.matrix_state_dep_params):
        raise ec.MissingParameterError("\"discrete_sse\"", "Parameter \"stash\" is not storing a valid macroevolutionary event handler.")

    if not stop_values_list:
        raise ec.MissingParameterError("\"discrete_sse\"", "Parameter \"stop_value\" is missing.")

    # TODO: have DnSSE take a prob_handler, and populate it if it's None upon initialization
    # all unit test should still work if this initialization is done correctly
    return dnsse.DnSSE(
        stash,
        stop_values_list,
        n=n_samples,
        n_replicates=n_repl,
        stop=stop_str,
        origin=origin,
        start_states_list=start_states_list,
        condition_on_speciation=cond_spn,
        condition_on_survival=cond_surv,
        condition_on_obs_both_sides_root=cond_obs_both_sides,
        min_rec_taxa=min_rec_taxa,
        max_rec_taxa=max_rec_taxa,
        abort_at_obs=abort_at_obs,
        epsilon=eps,
        runtime_limit=runtime_limit
    )


if __name__ == "__main__":
    # can be called from interface/grammar/
    # $ python3 make_dn_discrete_sse.py
    #
    # can also be called from phylojunction/
    # $ python3 interface/grammar/make_dn_discrete_sse.py
    # or
    # $ python3 -m interface.grammar.make_dn_discrete_sse
    #
    # can also be called from VS Code, if open folder is phylojuction/

    total_n_states = 2

    rates_t0_s0 = [
        sseobj.DiscreteStateDependentRate(
            name="lambda0",
            val=0.5,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[0, 0, 0]),
        sseobj.DiscreteStateDependentRate(
            name="mu0",
            val=0.25,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[0]),
        sseobj.DiscreteStateDependentRate(
            name="q01",
            val=0.75,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[0, 1])]

    rates_t0_s1 = [
        sseobj.DiscreteStateDependentRate(
            name="lambda1",
            val=1.5,
            event=sseobj.MacroevolEvent.W_SPECIATION,
            states=[1, 1, 1]),
        sseobj.DiscreteStateDependentRate(
            name="mu1",
            val=0.25,
            event=sseobj.MacroevolEvent.EXTINCTION,
            states=[1]),
        sseobj.DiscreteStateDependentRate(
            name="q10",
            val=0.75,
            event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION,
            states=[1, 0])]

    rates_t0 = rates_t0_s0 + rates_t0_s1

    # 1D: time slices (i) , 2D: all rates from all states in i-th time slice
    matrix_atomic_rate_params = [rates_t0]

    state_dep_rates_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_atomic_rate_params, total_n_states)

    event_handler = sseobj.MacroevolEventHandler(state_dep_rates_manager)

    sse_stash = sseobj.SSEStash(event_handler)

    # det_nd_pgm = pgm.DeterministicNodePGM(
        # "events", value=event_handler, parent_nodes=None)
    det_nd_pgm = pgm.DeterministicNodePGM(
        "sse_stash",
        value=sse_stash,
        parent_nodes=None)

    dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]] = dict()
    dn_param_dict["n"] = ["1"]
    dn_param_dict["nr"] = ["1"]
    dn_param_dict["stash"] = [det_nd_pgm]
    dn_param_dict["start_state"] = ["0"]
    dn_param_dict["stop"] = ["\"size\""]
    dn_param_dict["stop_value"] = ["10"]
    dn_param_dict["origin"] = ["\"true\""]
    dn_param_dict["cond_spn"] = ["\"false\""]
    dn_param_dict["cond_surv"] = ["\"true\""]

    discrete_sse_dn = make_discrete_SSE_dn(dn_param_dict)
    print(discrete_sse_dn.DN_NAME)
