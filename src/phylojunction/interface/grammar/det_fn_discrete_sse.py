import sys
sys.path.extend(["../../", "../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/interface/grammar/)
import typing as ty

# pj imports
import utility.exception_classes as ec
import pgm.pgm as pgm
import calculation.discrete_sse as sseobj # type: ignore

def make_SSEAtomicRate(det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]):
    """Create SSEAtomicRate as prompted by deterministic function call. Still no support for vectorization.

    Args:
        det_fn_param_dict (dict): dictionary containing parameter strings as keys and argument strings/NodePGM as value(s)
    """

    def extract_value_from_pgmnodes(pgm_node_list: ty.List[ty.Union[str, pgm.NodePGM]]) -> ty.List[ty.List[float]]:
        """_summary_

        Args:
            pgm_node_list (NodePGM): List of NodePGM objects (typing includes str because of type-safety)

        Raises:
            ec.NoPlatingAllowedError: _description_
            ec.SSEAtomicRateMisspec: _description_

        Returns:
            ty.List[ty.List[float]]: _description_
        """
        many_nodes_pgm = len(pgm_node_list) > 1
        v_list: ty.List[ty.List[float]] = []

        for node_pgm in pgm_node_list:
            
            # so mypy won't complain
            if isinstance(node_pgm, pgm.NodePGM):
                
                # no plating supported
                if node_pgm.repl_size > 1:
                    raise ec.NoPlatingAllowedError("sse_rate", node_pgm.node_pgm_name)
                
                v = node_pgm.value # list (I think before I also allowed numpy.ndarray, but not anymore)
                
                # so mypy won't complain
                if isinstance(v, list):
                    if len(v) > 1 and many_nodes_pgm:
                        raise ec.SSEAtomicRateMisspec(message="If many variables are passed as arguments to initialize another variable, each of these variables can contain only a single value. Exiting...")
                    elif len(v) == 1 and many_nodes_pgm:
                        v_list += v # making list longer (v should be a list, which is why I don't use append)
                    elif len(v) >= 1 and not many_nodes_pgm:
                        return v
        
        return v_list


    if not det_fn_param_dict:
        raise ec.NoSpecificationError(message="Cannot initialize SSE rate without specifications. Exiting...")

    event_type = None
    the_states = None
    # val is a list
    for arg, val in det_fn_param_dict.items():
        if arg == "value":
            if not val: raise ec.SSEAtomicRateMisspec(message="Cannot initialize SSE rate without value. Exiting...")
            
            value: ty.List[ty.List[float]] = []
            # val is a list of random variable objects
            # if type(val[0]) != str:
            if isinstance(val[0], pgm.NodePGM):
                value = extract_value_from_pgmnodes(val)
            
            # val is a list of strings
            else: value = val

        elif arg == "name" and not val: raise ec.SSEAtomicRateMisspec(message="Cannot initialize SSE rate without name. Exiting...")
        
        elif arg == "event":
            if not val: raise ec.SSEAtomicRateMisspec(message="Cannot initialize SSE rate without event type (e.g., \"w_speciation\"). Exiting...")
            
            # TODO: raise error if arg is not a string with quotes before and after

            # TODO: deal with vectorization later
            if val[0] == "\"w_speciation\"": event_type = sseobj.MacroevolEvent.W_SPECIATION
            elif val[0] == "\"bw_speciation\"": event_type = sseobj.MacroevolEvent.BW_SPECIATION
            elif val[0] == "\"asym_speciation\"": event_type = sseobj.MacroevolEvent.ASYM_SPECIATION
            elif val[0] == "\"extinction\"": event_type = sseobj.MacroevolEvent.EXTINCTION
            elif val[0] == "\"transition\"": event_type = sseobj.MacroevolEvent.ANAGENETIC_TRANSITION

        elif arg == "states":
            # TODO: raise error if arg is not a string with quotes before and after
            # if type(val) != str: raise SSEAtomicRateMisspec(message="Cannot initialize SSEAtomicRate because argument for \'states\'" + \
            #     " must be a string bounded by quotes without, with states separated by \'-\' (e.g., \"2-0-1\"). Exiting...")

            # make it a list of integers
            if not event_type == sseobj.MacroevolEvent.EXTINCTION:
                the_states = val

    if the_states:
        # TODO: deal with vectorization later
        return sseobj.AtomicSSERateParameter(name=det_fn_param_dict["name"][0], val=value, event=event_type, states=the_states)
    
    else:
        # TODO: deal with vectorization later
        return sseobj.AtomicSSERateParameter(name=det_fn_param_dict["name"][0], val=value, event=event_type)


def make_MacroEvolEventHandler(det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) -> sseobj.MacroEvolEventHandler:

    if not det_fn_param_dict:
        raise ec.NoSpecificationError(message="Cannot initialize SSE stash without specifications. Exiting...")

    _n_states: int = 1
    _n_time_slices: int = 1
    _time_slice_age_ends: ty.List[str] = None
    _seed_age_for_time_slicing: float = None
    _flat_rate_mat = None
    # val is a list
    for arg, val in det_fn_param_dict.items():
        if arg == "n_states":
            try:
                _n_states = int(val[0]) # TODO: Forbid different number of states here, so if len(val[0]) > 1, throw exception)
            
            except:
                raise ec.FunctionArgError(arg, "Was expecting an integer for \'n_states\'. Function was sse_stash().")

        if arg == "n_epochs":
            try:
                _n_time_slices = int(val[0]) # TODO: Forbid different number of epochs here, so if len(val[0]) > 1, throw exception)
            except:
                raise ec.FunctionArgError(arg, "Was expecting an integer for \'n_epochs\'. Function was sse_stash().")

        if arg == "seed_age":
            try:
                # TODO: Forbid different seed ages, because I will forbid different time slice age ends
                # (if you allow different seed ages, those might be younger than time slice ends
                _seed_age_for_time_slicing = float(val[0])
            
            except:
                raise ec.FunctionArgError(arg, "Was expecting a float. Function was sse_stash().")

        if arg == "epoch_age_ends":
            try: 
                _time_slice_age_ends = [float(v) for v in val]

            except: 
                raise ec.FunctionArgError(arg, "Could not convert epoch bound ages to floats. Function was sse_stash().")

            if len(val) != (_n_time_slices - 1):
                print(len(val))
                print(_n_time_slices)
                raise ec.FunctionArgsMismatchError("\"sse_wrap\" expects that the number of epoch ends is equal to the number of epochs minus 1")
    
    try:
        _flat_rate_mat = det_fn_param_dict["flat_rate_mat"] # list of NodePGM's

        # number of rates has to be divisible by number of slices
        if len(_flat_rate_mat) % _n_time_slices != 0:
            raise ec.WrongDimensionError(arg, len(_flat_rate_mat))
    
    except:
        raise ec.SSEStashMisspec(message="Cannot initialize SSE stash without a vector of SSE rates. Exiting...")
        
    _matrix_det_node_atomic_rate_params: ty.List[pgm.DeterministicNodePGM] = []
    _matrix_atomic_rate_params: ty.List[ty.List[sseobj.AtomicSSERateParameter]] = []
    _total_n_rate_params = len(_flat_rate_mat)
    _n_rate_params_per_slice = int(_total_n_rate_params / _n_time_slices)
    
    # iterating over time slices
    for i in range(0, _total_n_rate_params, _n_rate_params_per_slice):

        atomic_rate_params = list()
        for atomic_rate_param_det_nd in _flat_rate_mat[i:(i+_n_rate_params_per_slice)]:
            atomic_rate_params.append(atomic_rate_param_det_nd.value)
        
        _matrix_atomic_rate_params.append(atomic_rate_params) # appending a list
                                                             # 1D: time slices, 2D: atomic SSE rate list
        
        _matrix_det_node_atomic_rate_params.append(_flat_rate_mat[i:(i+_n_rate_params_per_slice)]) # list of lists of DeterministicNodePGM's
    
    # matrix_atomic_rate_params = [det_node_atomic_rate.value for det_node_atomic_rate in _matrix_det_node_atomic_rate_params]
        
    # fig_rates_manager = FIGRatesManager(_matrix_det_node_atomic_rate_params, _n_states)
    fig_rates_manager = sseobj.FIGRatesManager(_matrix_atomic_rate_params, _n_states, seed_age_for_time_slicing=_seed_age_for_time_slicing, list_time_slice_age_ends=_time_slice_age_ends)
        
    return sseobj.MacroEvolEventHandler(fig_rates_manager)