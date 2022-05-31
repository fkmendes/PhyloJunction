import typing as ty

# pj imports
import phylojunction.distribution.dn_discrete_sse as dnsse
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

def make_discrete_SSE_dn(dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) -> pgm.DistributionPGM:

    ############################
    # Validating dn_param_dict #
    ############################
    if not dn_param_dict:
        raise ec.NoSpecificationError(message="Cannot initialize discrete SSE distribution without specifications. Exiting...")

    #############################
    # IMPORTANT: Default values #
    #############################
    # default values for args that are not mandatory
    # all remaining args must be specified
    _n: int = 1
    _n_repl: int = 1
    _event_handler: sseobj.MacroEvolEventHandler = sseobj.MacroEvolEventHandler(sseobj.FIGRatesManager([[]], 1))
    _stop: str
    _stop_values_list: ty.List[float] = []
    origin: bool = True
    cond_spn: bool = False
    cond_surv: bool = True
    _start_states_list: ty.List[int] = []
    eps: float = 1e-12
    runtime_limit: int = 5 # 5 minutes

    # input validation (only things that are not already done by DnSSE)
    # NOTE: val is a list!
    for arg, val in dn_param_dict.items():
        if val:
            _extracted_val = pgm.extract_value_from_nodepgm(val) # if element in val is string, it remains unchanged, if NodePGM, we get its string-fied value

            ############################
            # Non-vectorized arguments #
            ############################
            # ... thus using only the first value!
            _first_val: ty.Union[str, pgm.NodePGM]
            if len(_extracted_val) >= 1:
                _first_val = _extracted_val[0] 
            # DeterministicNodePGM is in val, e.g., val = [pgm.DeterministicNodePGM]
            else:
                _first_val = val[0]
            
            if arg == "meh" and isinstance(_first_val, pgm.NodePGM):
                nodepgm_val = _first_val.value
                
                if isinstance(nodepgm_val, sseobj.MacroEvolEventHandler):
                # there should be only one event handler always, but it will be in list
                    _event_handler = nodepgm_val

            elif arg in ("n", "nr", "runtime_limit"):
                try:
                    if isinstance(_first_val, str):
                        int_val = int(_first_val) # no vectorization allowed here
                except:
                    raise ec.FunctionArgError(arg, "Was expecting an integer. Distribution discrete_sse() could not be initialized.")

                # if user specified n or runtime_limit, we use it, otherwise defaults are used
                if arg == "n": _n = int_val
                if arg == "runtime_limit": runtime_limit = int_val
                if arg == "nr": _n_repl = int_val

            elif arg == "stop":
                if isinstance(_first_val, str):
                    _stop = _first_val.replace("\"", "")

            elif arg in ("origin", "cond_spn", "cond_surv"):
                if _first_val != "\"true\"" and _first_val != "\"false\"":
                    raise ec.FunctionArgError(arg, "Was expecting either \"true\" or \"false\". Distribution discrete_sse() could not be initialized.")

                else:
                    if isinstance(_first_val, str):
                        parsed_val = _first_val.replace("\"", "")

                    if arg == "origin" and parsed_val == "true": origin = True
                    elif arg == "origin" and parsed_val == "false": origin = False
                    elif arg == "cond_spn" and parsed_val == "true": cond_spn = True
                    elif arg == "cond_spn" and parsed_val == "false": cond_spn = False
                    elif arg == "cond_surv" and parsed_val == "true": cond_surv = True
                    elif arg == "cond_surv" and parsed_val == "false": cond_surv = False

            # if user specified epsilon, we use it, otherwise default is used
            elif arg == "eps":
                try:
                    if isinstance(_first_val, str):
                        float_val = float(_first_val) # no vectorization allowed here
                except:
                    raise ec.FunctionArgError(arg, "Was expecting a double. Distribution discrete_sse() could not be initialized.")

                eps = float_val

            #####################
            # Vectorized values #
            #####################
            elif arg == "stop_value":
                _is_negative = False
                try:
                    _stop_values_list = [float(v) for v in _extracted_val if isinstance(v, str)]
                    
                    for sv in _stop_values_list:
                        if sv < 0.0:
                            _is_negative = True

                except:
                    if _is_negative:
                        raise ec.FunctionArgError(arg, "Stop value cannot be negative. Distribution discrete_sse() could not be initialized.")

                    raise ec.FunctionArgError(arg, "Was expecting number for stopping value(s). Distribution discrete_sse() could not be initialized.")


            elif arg == "start_state":
                try:
                    _start_states_list = [int(v) for v in _extracted_val if isinstance(v, str)]
                except:
                    raise ec.FunctionArgError(arg, "Was expecting integers for starting states. Distribution discrete_sse() could not be initialized.")

    
    # making sure essential parameters of distribution have been specified
    if not any(_event_handler.fig_rates_manager.atomic_rate_params_matrix):
        raise ec.DnInitMisspec("\"discrete_sse\"", "Parameter \"meh\" is missing.")
    if not _stop_values_list:
        raise ec.DnInitMisspec("\"discrete_sse\"", "Parameter \"stop_value\" is missing.")
        
    return dnsse.DnSSE(
        _event_handler,
        _stop_values_list,
        n=_n,
        n_replicates=_n_repl,
        stop=_stop,
        origin=origin,
        start_states_list=_start_states_list,
        condition_on_speciation=cond_spn,
        condition_on_survival=cond_surv,
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

    rates_t0_s0 = [ sseobj.AtomicSSERateParameter(name="lambda0", val=0.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[0,0,0]),
                        sseobj.AtomicSSERateParameter(name="mu0", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[0]),
                        sseobj.AtomicSSERateParameter(name="q01", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[0,1]) ]

    rates_t0_s1 = [ sseobj.AtomicSSERateParameter(name="lambda1", val=1.5, event=sseobj.MacroevolEvent.W_SPECIATION, states=[1,1,1]),
                        sseobj.AtomicSSERateParameter(name="mu1", val=0.25, event=sseobj.MacroevolEvent.EXTINCTION, states=[1]),
                        sseobj.AtomicSSERateParameter(name="q10", val=0.75, event=sseobj.MacroevolEvent.ANAGENETIC_TRANSITION, states=[1,0]) ]

    rates_t0 = rates_t0_s0 + rates_t0_s1

    matrix_atomic_rate_params = [ rates_t0 ] # 1D: time slices (i) , 2D: all rates from all states in i-th time slice

    fig_rates_manager = sseobj.FIGRatesManager(matrix_atomic_rate_params, total_n_states)

    event_handler = sseobj.MacroEvolEventHandler(fig_rates_manager)

    det_nd_pgm = pgm.DeterministicNodePGM("events", value=event_handler, parent_nodes=None)

    dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]] = dict()
    dn_param_dict["n"] = ["1"]
    dn_param_dict["nr"] = ["1"]
    dn_param_dict["meh"] = [det_nd_pgm]
    dn_param_dict["start_state"] = ["0"]
    dn_param_dict["stop"] = ["\"size\""]
    dn_param_dict["stop_value"] = ["10"]
    dn_param_dict["origin"] = ["\"true\""]
    dn_param_dict["cond_spn"] = ["\"false\""]
    dn_param_dict["cond_surv"] = ["\"true\""]

    discrete_sse_dn = make_discrete_SSE_dn(dn_param_dict)
    print(discrete_sse_dn.DN_NAME)