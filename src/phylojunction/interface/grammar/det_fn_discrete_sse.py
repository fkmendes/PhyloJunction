import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

def extract_value_from_pgmnodes(pgm_node_list: ty.List[pgm.NodePGM]) \
        -> ty.List[float]:
        """_summary_

        Args:
            pgm_node_list (NodePGM): List of NodePGM objects (typing includes str because of type-safety)

        Raises:
            ec.NoPlatingAllowedError: _description_
            ec.StateDependentParameterMisspec: _description_

        Returns:
            ty.List[ty.List[float]]: _description_
        """
        many_nodes_pgm = len(pgm_node_list) > 1
        v_list: ty.List[float] = []

        for node_pgm in pgm_node_list:
            
            # so mypy won't complain
            if isinstance(node_pgm, pgm.NodePGM):
                
                # no plating supported
                if node_pgm.repl_size > 1:
                    raise ec.NoPlatingAllowedError("sse_rate", node_pgm.node_name)
                
                v = node_pgm.value # list (I think before I also allowed numpy.ndarray, but not anymore)
                
                # so mypy won't complain
                if isinstance(v, list):
                    if len(v) > 1 and many_nodes_pgm:
                        raise ec.StateDependentParameterMisspec(message="If many variables are passed as arguments to initialize another variable, each of these variables can contain only a single value. Exiting...")
                    elif len(v) == 1 and many_nodes_pgm:
                        v_list += v # making list longer (v should be a list, which is why I don't use append)
                    elif len(v) >= 1 and not many_nodes_pgm:
                        return v
        
        return v_list


def make_DiscreteStateDependentRate(det_fn_param_dict: \
    ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
    -> sseobj.DiscreteStateDependentRate:
    """Create SSEAtomicRate as prompted by deterministic function call

    Args:
        det_fn_param_dict (dict): dictionary containing parameter
            strings as keys, and lists of either strings or NodePGMs
            as value(s)
    """

    if not det_fn_param_dict:
        raise ec.NoSpecificationError(
            message="Cannot initialize state-dependent rate without " \
                + "specifications.")

    event_type = None
    value: ty.List[float] = []
    the_states: ty.List[int] = []

    ######################################
    # Checking for input correctness,    #
    # and preparing for instantiation of #
    # parameter objects                  #
    ######################################

    # val is a list
    for arg, val in det_fn_param_dict.items():
        if arg == "value":
            if not val:
                raise ec.StateDependentParameterMisspec(
                    message="Cannot initialize SSE rate without value. " \
                        "Exiting...")
            
            # val is a list of random variable objects
            # if type(val[0]) != str:
            if isinstance(val[0], pgm.NodePGM):
                # need to declare cast_val separately so mypy won't complain
                cast_val1: ty.List[pgm.NodePGM] = \
                    ty.cast(ty.List[pgm.NodePGM], val)
                value = extract_value_from_pgmnodes(cast_val1)
            
            # val is a list of strings
            elif isinstance(val[0], str):
                cast_val2: ty.List[str] = ty.cast(ty.List[str], val)
                value = [float(v) for v in cast_val2]

        elif arg == "name" and not val:
            raise ec.StateDependentParameterMisspec(
                message="Cannot initialize state-dependent rate without name.")
        
        elif arg == "event":
            if not val:
                raise ec.StateDependentParameterMisspec(
                    message="Cannot initialize state-dependent rate without " \
                        + "event type (e.g., \"w_speciation\"). Exiting...")
            
            # TODO: raise error if arg is not a string with quotes before and after

            # TODO: deal with vectorization later
            if val[0] == "\"w_speciation\"": event_type = sseobj.MacroevolEvent.W_SPECIATION
            elif val[0] == "\"bw_speciation\"": event_type = sseobj.MacroevolEvent.BW_SPECIATION
            elif val[0] == "\"asym_speciation\"": event_type = sseobj.MacroevolEvent.ASYM_SPECIATION
            elif val[0] == "\"extinction\"": event_type = sseobj.MacroevolEvent.EXTINCTION
            elif val[0] == "\"transition\"": event_type = sseobj.MacroevolEvent.ANAGENETIC_TRANSITION
            elif val[0] == "\"anc_sampling\"": event_type = sseobj.MacroevolEvent.ANCESTOR_SAMPLING

        elif arg == "states":
            # need to declare cast_val separately
            # so mypy won't complain
            the_states = ty.cast(ty.List[int], val)

    ##################################
    # Istantiating parameter objects #
    ##################################

    sse_rate_name = det_fn_param_dict["name"][0]
    
    if the_states:
        # TODO: deal with vectorization later
        if isinstance(sse_rate_name, str):
            # see 'invariance vs covariance' in mypy's doc for why list() is called here
            return sseobj.DiscreteStateDependentRate(name=sse_rate_name, val=list(value), event=event_type, states=the_states)
    
    else:
        # TODO: deal with vectorization later
        if isinstance(sse_rate_name, str):
            # see 'invariance vs covariance' in mypy's doc for why list() is called here
            return sseobj.DiscreteStateDependentRate(name=sse_rate_name, val=list(value), event=event_type)


def make_DiscreteStateDependentProbability(det_fn_param_dict: \
    ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
    -> sseobj.DiscreteStateDependentProbability:
    """
    Create DiscreteStateDependentProbability as prompted by
        deterministic function call

    Args:
        det_fn_param_dict (dict): dictionary containing parameter
            strings as keys, and lists of either strings or NodePGMs
            as value(s)
    """

    if not det_fn_param_dict:
        raise ec.NoSpecificationError(
            message="Cannot initialize state-dependent probability " \
                "without specifications.")

    value: ty.List[float] = []
    state: ty.Optional[int] = None

    ######################################
    # Checking for input correctness,    #
    # and preparing for instantiation of #
    # parameter objects                  #
    ######################################

     # val is a list
    for arg, val in det_fn_param_dict.items():
        if arg == "value":
            if not val:
                raise ec.StateDependentParameterMisspec(
                    message="Cannot initialize state-dependent probability" \
                        + " without value.")
            
            # val is a list of random variable objects
            # if type(val[0]) != str:
            if isinstance(val[0], pgm.NodePGM):
                # need to declare cast_val separately so mypy won't complain
                cast_val1: ty.List[pgm.NodePGM] = \
                    ty.cast(ty.List[pgm.NodePGM], val)
                value = extract_value_from_pgmnodes(cast_val1)
            
            # val is a list of strings
            elif isinstance(val[0], str):
                cast_val2: ty.List[str] = ty.cast(ty.List[str], val)
                value = [float(v) for v in cast_val2]

        elif arg == "name" and not val:
            raise ec.StateDependentParameterMisspec(
                message="Cannot initialize state-dependent probability" \
                    + "without name.")

        elif arg == "state":
            # need to declare cast_val separately
            # so mypy won't complain
            if len(val) > 1:
                raise ec.RequireScalarError(
                    "StateDependentProbability",
                    "state"
                )

            state = int(val[0])

    ##################################
    # Istantiating parameter objects #
    ##################################
    
    stade_dep_prob_name = det_fn_param_dict["name"][0]
    
    if isinstance(stade_dep_prob_name, str):
        if state == None:
            return sseobj.DiscreteStateDependentProbability(
                    name=stade_dep_prob_name,
                    val=list(value)
                    )

        else:
            return sseobj.DiscreteStateDependentProbability(
                name=stade_dep_prob_name,
                val=list(value),
                state=state
                )


def make_MacroevolEventHandler(
    det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
        -> sseobj.MacroevolEventHandler:

    if not det_fn_param_dict:
        raise ec.NoSpecificationError(message="Cannot initialize SSE stash without specifications. Exiting...")

    n_states: int = 1
    n_time_slices: int = 1
    time_slice_age_ends: ty.List[float] = []
    seed_age_for_time_slicing: ty.Optional[float] = None
    flat_state_dep_rate_mat: ty.List[pgm.DeterministicNodePGM]
    flat_state_dep_prob_mat: ty.List[pgm.DeterministicNodePGM] = []
    
    # val is a list
    for arg, val in det_fn_param_dict.items():
        # so mypy won't complain
        if isinstance(val[0], str):
            if arg == "n_states":
                try:
                    n_states = int(val[0]) # TODO: Forbid different number of states here, so if len(val[0]) > 1, throw exception)
                
                except:
                    raise ec.FunctionArgError(arg, "Was expecting an integer for \'n_states\'. Function was sse_stash().")

            if arg == "n_epochs":
                try:
                    n_time_slices = int(val[0]) # TODO: Forbid different number of epochs here, so if len(val[0]) > 1, throw exception)
                
                except:
                    raise ec.FunctionArgError(arg, "Was expecting an integer for \'n_epochs\'. Function was sse_wrap().")

            if arg == "seed_age":
                try:
                    # TODO: Forbid different seed ages, because I will forbid different time slice age ends
                    # (if you allow different seed ages, those might be younger than time slice ends
                    seed_age_for_time_slicing = float(val[0])
                
                except:
                    raise ec.FunctionArgError(arg,
                    "Was expecting a float. Function was sse_wrap().")

            if arg == "epoch_age_ends":
                try:
                    time_slice_age_ends = \
                        [float(v) for v in val if isinstance(v, str)]

                except: 
                    raise ec.FunctionArgError(arg,
                        "Could not convert epoch bound ages to floats. " \
                        + "Function was sse_wrap().")

                if len(val) != (n_time_slices - 1):
                    raise ec.FunctionArgsMismatchError(
                        "sse_wrap",
                        "\"sse_wrap\" expects that the number of epoch ends " \
                        + "is equal to the number of epochs minus 1")

            if arg == "flat_prob_mat":
                if det_fn_param_dict["flat_prob_mat"]:
                    flat_state_dep_prob_mat = \
                        [v for v in det_fn_param_dict["flat_prob_mat"] \
                            if isinstance(v, pgm.DeterministicNodePGM)]

                # total number of rates has to be divisible by number of slices
                if len(flat_state_dep_prob_mat) % n_time_slices != 0:
                    raise ec.WrongDimensionError(arg, len(flat_state_dep_prob_mat))

    try:
        flat_state_dep_rate_mat = \
            [v for v in det_fn_param_dict["flat_rate_mat"] \
                if isinstance(v, pgm.DeterministicNodePGM)] # list of NodePGM's

        # total number of rates has to be divisible by number of slices
        if len(flat_state_dep_rate_mat) % n_time_slices != 0:
            raise ec.WrongDimensionError(arg, len(flat_state_dep_rate_mat))
    
    except:
        raise ec.SSEStashMisspec(
            message="Cannot initialize SSE stash without a vector of " \
                + "state-dependent rates.")
        
    matrix_state_dep_rates: \
        ty.List[ty.List[sseobj.DiscreteStateDependentRate]] = []
    total_n_rate = len(flat_state_dep_rate_mat)
    n_rates_per_slice = int(total_n_rate / n_time_slices)

    matrix_state_dep_probs: \
        ty.List[ty.List[sseobj.DiscreteStateDependentProbability]] = []
    expected_n_prob_per_slice = n_states * n_time_slices
    total_n_prob = n_states * n_time_slices # initialize to expected
    n_prob_per_slice = int(total_n_prob / n_time_slices)

    # must provide one probability parameter per state
    # per time slice!
    probs_were_provided = False
    if len(flat_state_dep_prob_mat) > 0 and \
        (n_prob_per_slice != expected_n_prob_per_slice):
        probs_were_provided = True
        # TODO: write exception
        pass
    
    # populating state-dependent rate matrix #
    for i in range(0, total_n_rate, n_rates_per_slice):

        state_dep_rates: ty.List[sseobj.DiscreteStateDependentRate] = []
        for state_dep_rate_det_nd in flat_state_dep_rate_mat[i:(i+n_rates_per_slice)]:

            # so mypy won't complain
            if isinstance(state_dep_rate_det_nd.value, sseobj.DiscreteStateDependentRate):
                state_dep_rates.append(state_dep_rate_det_nd.value)
        
        # appending a list of DiscreteStateDependentRate
        # 1D: time slices, 2D: state-dep rate list
        matrix_state_dep_rates.append(state_dep_rates) 

    state_dep_rates_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_state_dep_rates,
            n_states,
            seed_age_for_time_slicing=seed_age_for_time_slicing,
            list_time_slice_age_ends=time_slice_age_ends
        )

    # populating state-dependent prob matrix #
    #
    # this iteration here is effectively over
    # time-slices if dimensions are indeed all
    # correct
    for i in range(0, total_n_prob, n_prob_per_slice):

        state_dep_probs: ty.List[sseobj.DiscreteStateDependentProbability] = []
        if probs_were_provided:
            for state_dep_prob_det_nd in flat_state_dep_prob_mat[i:(i+n_prob_per_slice)]:

                # so mypy won't complain
                if isinstance(state_dep_prob_det_nd.value,
                    sseobj.DiscreteStateDependentProbability):
                    state_dep_probs.append(state_dep_prob_det_nd.value)

        # else here means we will have i being
        # one integer per probability (per state
        # per slice)
        #
        # we will make all probabilities 1.0 by
        # default
        else:
            st = i % n_prob_per_slice
            state_dep_prob = sseobj.DiscreteStateDependentProbability(
                name="rho" + str(st) + "_t" + str(i // n_prob_per_slice),
                val=1.0, state=st)
                
            # so mypy won't complain
            if isinstance(state_dep_prob,
                sseobj.DiscreteStateDependentProbability):
                state_dep_probs.append(state_dep_prob)

            # appending a list of DiscreteStateDependentProbability
            # 1D: time slices, 2D: state-dep prob list
            matrix_state_dep_probs.append(state_dep_probs) 

    state_dep_probs_manager = \
        sseobj.DiscreteStateDependentParameterManager(
            matrix_state_dep_probs,
            n_states,
            seed_age_for_time_slicing=seed_age_for_time_slicing,
            list_time_slice_age_ends=time_slice_age_ends
        )

    # state_dep_probs_manager = \
    #     sseobj.DiscreteStateDependentParameterManager(
    #         # TODO
    #         n_states,
    #         seed_age_for_time_slicing=seed_age_for_time_slicing,
    #         list_time_slice_age_ends=time_slice_age_ends
    #     )
        
    # return sseobj.MacroevolEventHandler(
    #     state_dep_rates_manager,
    #     state_dep_probs_manager,

    # )

    return sseobj.MacroevolEventHandler(state_dep_rates_manager)