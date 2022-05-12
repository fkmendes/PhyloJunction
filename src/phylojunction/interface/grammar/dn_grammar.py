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

    dn_grammar_dict: ty.Dict[str, ty.Tuple[str]]

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
    def extract_value_from_rvpgm(cls, rvpgm_list: ty.List[pgm.StochasticNodePGM]) -> ty.List[str]:
        try:
            return [str(i) for rvpgm in rvpgm_list for i in rvpgm.value]
        except:
            raise ec.VariableMisspec(str(rvpgm_list[0]))

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

        #################################
        #  Parametric distributions #
        #################################
        if dn_id == "lognormal":
            pars = [1, 1, [0.0], [1.0], True] # log-space True by default
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    # val = val[0] # TODO: deal with vectorization later
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                        val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("lognormal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[0] = int(val[0])
                    elif arg == "nr":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[1] = int(val[0])
                    elif arg == "mean":
                        pars[2] = [float(v) for v in val]
                        # pars[1] = float(val) # no vectorization
                    elif arg == "sd":
                        pars[3] = [float(v) for v in val]
                        # pars[2] = float(val) # no vectorization
                    elif arg == "log_space":
                        if val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            pars[4] = True
                        if val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            pars[4] = False

            return dnpar.DnLogNormal(pars, parent_node_tracker)

        if dn_id == "normal":
            pars = [1, 1, [0.0], [1.0]]
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                        val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("normal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[0] = int(val[0])
                    elif arg == "nr":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[1] = int(val[0])
                    elif arg == "mean":
                        pars[2] = [float(v) for v in val]
                    elif arg == "sd":
                        pars[3] = [float(v) for v in val]

            return dnpar.DnNormal(pars, parent_node_tracker)

        if dn_id == "exponential":
            pars = [1, 1, [1.0], True] # rate_parameterization True by default
            parent_node_tracker = dict() # { lambda: node_pgm1_name }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                        val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("exponential", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[0] = int(val[0])
                    elif arg == "nr":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[1] = int(val[0])
                    elif arg == "rate":
                        pars[2] = [float(v) for v in val]
                    elif arg == "rate_parameterization":
                        if val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            pars[3] = True
                        if val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            pars[3] = False

            return dnpar.DnExponential(pars, parent_node_tracker)

        if dn_id == "gamma":
            pars = [1, 1, [1.0], [1.0], False] # rate_parameterization False by default
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                        val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("gamma", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[0] = int(val[0])
                    elif arg == "nr":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[1] = int(val[0])
                    elif arg == "shape":
                        pars[2] = [float(v) for v in val]
                    elif arg == "scale":
                        pars[3] = [float(v) for v in val]

            return dnpar.DnGamma(pars, parent_node_tracker)

        if dn_id == "unif":
            pars = [1, 1, [0.0], [1.0]]
            parent_node_tracker = dict() # { mean: node_pgm1_name, sd: node_pgm2_name, ... }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_pgm_name # needed for building inference specifications
                        val = cls.extract_value_from_rvpgm(val)

                    if not cls.grammar_check("unif", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[0] = int(val[0])
                    elif arg == "nr":
                        if len(val) > 1:
                            raise ec.RequireScalarError(arg)
                        pars[1] = int(val[0])
                    elif arg == "min":
                        pars[2] = [float(v) for v in val]
                    elif arg == "max":
                        pars[3] = [float(v) for v in val]

            return dnpar.DnUnif(pars, parent_node_tracker)


        #################################
        #  Non-parametric distributions #
        #################################
        if dn_id == "discrete_sse":

            # dn_param_dict is validated inside
            return cls.init_return_discrete_SSE_dn(dn_param_dict)