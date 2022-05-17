import sys
sys.path.extend(["../../", "../", "../phylojunction"])
import typing as ty

# pj imports
import distribution.dn_parametric as dnpar
import pgm.pgm as pgm
import make_dn_discrete_sse as make_dnsse
import utility.exception_classes as ec
# from user_interface.dn_discrete_sse import make_discrete_SSE_dn # https://stackoverflow.com/questions/16981921/relative-imports-in-python-3

class DnGrammar():

    dn_grammar_dict: ty.Dict[str, ty.Tuple[str, ...]]

    ################################
    #  All available distributions #
    ################################
    dn_grammar_dict = {
        "lognormal":
            tuple(["n", "nr", "mean", "sd", "log_space"]),
        "normal":
            tuple(["n", "nr", "mean", "sd"]),
        "exponential":
            tuple(["n", "nr", "rate", "rate_parameterization"]),
        "gamma":
            tuple(["n", "nr", "shape", "scale"]),
        "unif":
            tuple(["n", "nr", "min", "max"]),
        "discrete_sse": tuple(["n", "nr", "stop", "stop_value", "origin", "meh",
                             "start_state", "eps", "runtime_limit",
                             "cond_spn", "cond_surv"]) # meh = MacroEvol Handler
    }

    def __init__(self) -> None:
        pass

    @classmethod
    def grammar_check(cls, dn_id: str, dn_param: str) -> bool:
        # we ignore the clamp parameter-arg pair; they are handled by cmd_parse()
        if dn_param in cls.dn_grammar_dict[dn_id] or dn_param == "clamp":
            return True
        return False

    @classmethod
    def extract_value_from_rvpgm(cls, val_list: ty.List[ty.Union[str, pgm.StochasticNodePGM]]) -> ty.List[str]:
        """
        Return copy of val_list if all elements are strings representing values.
        When elements are StochasticNodePGMs, replaces those objects by their values cast to string.
        If StochasticNodePGMs objects do not have .value field or if they cannot be string-fied, 
        raise exception.
        """
        extracted_val_list: ty.List[str] = []
        for v in val_list:
            if isinstance(v, str):
                extracted_val_list.append(v)
            
            elif isinstance(v, pgm.StochasticNodePGM):
                try:
                    extracted_val_list.append(str(v.value))
                except:
                    raise ec.VariableMisspec(str(v))

        return extracted_val_list

    @classmethod
    def init_return_discrete_SSE_dn(cls, det_fn_param_dict) -> pgm.DistributionPGM:
        # dn_param_dict is validated inside
        return make_dnsse.make_discrete_SSE_dn(det_fn_param_dict) # returns dnssewrap

    @classmethod
    def create_dn_obj(cls, dn_id: str, dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.StochasticNodePGM]]]) -> pgm.DistributionPGM:
        """_summary_

        Args:
            dn_id (_type_): _description_
            dn_param_dict (_type_): _description_

        Returns:
            _type_: _description_
        """

        _extracted_val: ty.List[str]

        #################################
        #  Parametric distributions #
        #################################
        if dn_id == "lognormal":
            # pars = [1, 1, [0.0], [1.0], True] # log-space True by default
            _ln_n_draws: int
            _ln_n_repl: int
            _ln_mean: ty.List[float]
            _ln_sd: ty.List[float]
            _ln_log_space: bool
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    # val = val[0] # TODO: deal with vectorization later
                    
                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name 
                    
                    _extracted_val = cls.extract_value_from_rvpgm(val) # if element in val is string, it remains unchanged, if NodePGM, we get its string-fied value

                    if not cls.grammar_check("lognormal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnLogNormal.DN_NAME, arg)
                        # pars[0] = int(val[0])
                        _ln_n_draws = int(_extracted_val[0]) # only one element always
                    elif arg == "nr":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnLogNormal.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        _ln_n_repl = int(_extracted_val[0]) # only one element always
                    elif arg == "mean":
                        # pars[2] = [float(v) for v in _extracted_val]
                        _ln_mean = [float(v) for v in _extracted_val]
                    elif arg == "sd":
                        # pars[3] = [float(v) for v in _extracted_val]
                        _ln_sd = [float(v) for v in _extracted_val]
                    elif arg == "log_space":
                        if _extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[4] = True
                            _ln_log_space = True
                        elif _extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[4] = False
                            _ln_log_space = False

            # return dnpar.DnLogNormal(pars, parent_node_tracker)
            return dnpar.DnLogNormal(_ln_n_draws, _ln_n_repl, _ln_mean, _ln_sd, _ln_log_space, parent_node_tracker)

        elif dn_id == "normal":
            # pars = [1, 1, [0.0], [1.0]]
            _norm_n_draws: int
            _norm_n_repl: int
            _norm_mean: ty.List[float]
            _norm_sd: ty.List[float]
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    
                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name
                    
                    _extracted_val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("normal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnNormal.DN_NAME, arg)
                        # pars[0] = int(extracted_val[0])
                        _norm_n_draws = int(_extracted_val[0])
                    elif arg == "nr":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnNormal.DN_NAME, arg)
                        # pars[1] = int(extracted_val[0])
                        _norm_n_repl = int(_extracted_val[0])
                    elif arg == "mean":
                        # pars[2] = [float(v) for v in extracted_val]
                        _norm_mean = [float(v) for v in _extracted_val]
                    elif arg == "sd":
                        # pars[3] = [float(v) for v in extracted_val]
                        _norm_sd = [float(v) for v in _extracted_val]

            # return dnpar.DnNormal(pars, parent_node_tracker)
            return dnpar.DnNormal(_norm_n_draws, _norm_n_repl, _norm_mean, _norm_sd, parent_node_tracker)

        elif dn_id == "exponential":
            # pars = [1, 1, [1.0], True] # rate_parameterization True by default
            _exp_n_draws: int
            _exp_n_repl: int
            _exp_scale_or_rate: ty.List[float]
            _exp_rate_parameterization: bool
            parent_node_tracker = dict() # { lambda: node_pgm1_name }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                    
                    _extracted_val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("exponential", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnExponential.DN_NAME, arg)
                        # pars[0] = int(_extracted_val[0])
                        _exp_n_draws = int(_extracted_val[0])
                    elif arg == "nr":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnExponential.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        _exp_n_repl = int(_extracted_val[0])
                    elif arg == "rate":
                        # pars[2] = [float(v) for v in _extracted_val]
                        _exp_scale_or_rate = [float(v) for v in _extracted_val]
                    elif arg == "rate_parameterization":
                        if _extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[3] = True
                            _exp_rate_parameterization = True
                        elif _extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[3] = False
                            _exp_rate_parameterization = False

            # return dnpar.DnExponential(pars, parent_node_tracker)
            return dnpar.DnExponential(_exp_n_draws, _exp_n_repl, _exp_scale_or_rate, _exp_rate_parameterization, parent_node_tracker)

        elif dn_id == "gamma":
            # pars = [1, 1, [1.0], [1.0], False] # rate_parameterization False by default
            _gamma_n_draws: int
            _gamma_n_repl: int
            _gamma_shape: ty.List[float]
            _gamma_scale_or_rate: ty.List[float]
            _gamma_rate_parameterization: bool
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name
                    
                    _extracted_val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("gamma", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnGamma.DN_NAME, arg)
                        # pars[0] = int(_extracted_val[0])
                        _gamma_n_draws = int(_extracted_val[0])
                    elif arg == "nr":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnGamma.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        _gamma_n_repl = int(_extracted_val[0])
                    elif arg == "shape":
                        # pars[2] = [float(v) for v in _extracted_val]
                        _gamma_shape = [float(v) for v in _extracted_val]
                    elif arg == "scale":
                        # pars[3] = [float(v) for v in _extracted_val]
                        _gamma_scale_or_rate = [float(v) for v in _extracted_val]
                    elif arg == "rate_parameterization":
                        if _extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[3] = True
                            _gamma_rate_parameterization = True
                        elif _extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[3] = False
                            _gamma_rate_parameterization = False

            # return dnpar.DnGamma(pars, parent_node_tracker)
            return dnpar.DnGamma(_gamma_n_draws, _gamma_n_repl, _gamma_shape, _gamma_scale_or_rate, _gamma_rate_parameterization, parent_node_tracker)

        elif dn_id == "unif":
            # pars: ty.List[ty.Union[int, ty.List[float]]] = [1, 1, [0.0], [1.0]]
            _unif_n_draws: int
            _unif_n_repl: int
            _unif_min: ty.List[float]
            _unif_max: ty.List[float]
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    
                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                    
                    _extracted_val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("unif", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnUnif.DN_NAME, arg)
                        # pars[0] = int(_extracted_val[0])
                        _unif_n_draws = int(_extracted_val[0])
                    elif arg == "nr":
                        if len(_extracted_val) > 1:
                            raise ec.RequireScalarError(dnpar.DnUnif.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        _unif_n_repl = int(_extracted_val[0])
                    elif arg == "min":
                        # pars[2] = [float(v) for v in _extracted_val]
                        _unif_min = [float(v) for v in _extracted_val]
                    elif arg == "max":
                        # pars[3] = [float(v) for v in _extracted_val]
                        _unif_max = [float(v) for v in _extracted_val]

            # return dnpar.DnUnif(pars, parent_node_tracker)
            return dnpar.DnUnif(_unif_n_draws, _unif_n_repl, _unif_min, _unif_max, parent_node_tracker)


        #################################
        #  Non-parametric distributions #
        #################################
        
        # discrete_sse
        # elif dn_id == "discrete_sse": (switch to this one when there's a new one)
        else:
            # dn_param_dict is validated inside
            return cls.init_return_discrete_SSE_dn(dn_param_dict)