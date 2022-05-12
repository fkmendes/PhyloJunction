import sys
sys.path.extend(["../../", "../", "../phylojunction"])
import typing as ty

# pj imports
import distribution.dn_discrete_sse as dnsse
import calculation.discrete_sse as sseobj
import pgm.pgm as pgm
import utility.exception_classes as ec


def make_discrete_SSE_dn(dn_param_dict):

    ############################
    # Validating dn_param_dict #
    ############################
    if not dn_param_dict:
        raise ec.NoSpecificationError(message="Cannot initialize discrete SSE distribution without specifications. Exiting...")

    # default values for args that are not mandatory
    # all remaining args must be specified
    n_repl = 1
    event_handler = None
    stop, stop_value = str(), str()
    origin, cond_spn, cond_surv = True, False, True
    start_states_list = list()
    eps = 1e-12
    runtime_limit = 5 # 5 minutes

    # input validation (only things that are not already done by DnSSE)
    # NOTE: val is a list!
    for arg, val in dn_param_dict.items():
        if val:
            if arg == "meh":
                meh_det_output_pgm = dn_param_dict["meh"][0] # there should be only one event handler always, but it will be in list
                event_handler = meh_det_output_pgm.value

            if arg in ("n", "nr", "runtime_limit"):
                try: int_val = int(val[0]) # no vectorization allowed here
                except: raise ec.FunctionArgError(arg, "Was expecting an integer. Distribution discrete_sse() could not be initialized.")

                # if user specified n or runtime_limit, we use it, otherwise defaults are used
                if arg == "n": n = int_val
                if arg == "runtime_limit": runtime_limit = int_val
                if arg == "nr": n_repl = int_val

            if arg == "stop":
                stop = val[0].replace("\"", "")

            if arg == "stop_value":
                try: stop_value = float(val[0])
                except: raise ec.FunctionArgError(arg, "Was expecting a float. Distribution discrete_sse() could not be initialized.")

                stop_value = float(val[0])

            # TODO deal with vectorization later
            if arg in ("origin", "cond_spn", "cond_surv"):
                if val[0] != "\"true\"" and val[0] != "\"false\"":
                    raise ec.FunctionArgError(arg, "Was expecting either \"true\" or \"false\". Distribution discrete_sse() could not be initialized.")

                else:
                    parsed_val = val[0].replace("\"", "")

                    if arg == "origin" and parsed_val == "true": origin = True
                    elif arg == "origin" and parsed_val == "false": origin = False
                    elif arg == "cond_spn" and parsed_val == "true": cond_spn = True
                    elif arg == "cond_spn" and parsed_val == "false": cond_spn = False
                    elif arg == "cond_surv" and parsed_val == "true": cond_surv = True
                    elif arg == "cond_surv" and parsed_val == "false": cond_surv = False

            if arg == "start_state":
                start_states_list = [int(v) for v in val]

            # if user specified epsilon, we use it, otherwise default is used
            if arg == "eps":
                try: float_val = float(val[0]) # no vectorization allowed here
                except:
                    raise ec.FunctionArgError(arg, "Was expecting a double. Distribution discrete_sse() could not be initialized.")

                eps = float_val

    return dnsse.DnSSE(
        n=n,
        n_replicates=n_repl,
        stop=stop,
        stop_value=stop_value,
        origin=origin,
        event_handler=event_handler,
        start_states_list=start_states_list,
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

    dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.DeterministicNodePGM]]]
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