import typing as ty

# pj imports
import phylojunction.distribution.dn_parametric as dnpar
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.grammar.make_dn_discrete_sse as make_dnsse
import phylojunction.utility.exception_classes as ec
# from user_interface.dn_discrete_sse import make_discrete_SSE_dn # https://stackoverflow.com/questions/16981921/relative-imports-in-python-3

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class PJDnGrammar():

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
            tuple(["n", "nr", "shape", "scale", "rate_parameterization"]),
        "unif":
            tuple(["n", "nr", "min", "max"]),
        "discrete_sse": tuple(["n", "nr", "stop", "stop_value", "origin",
                               "stash", "start_state", "eps", "runtime_limit",
                               "cond_spn", "cond_surv"])
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
    def init_return_discrete_SSE_dn(
        cls,
        dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> pgm.DistributionPGM:
        # dn_param_dict is validated inside

        return make_dnsse.make_discrete_SSE_dn(dn_param_dict)

    @classmethod
    def create_dn_obj(
        cls,
        dn_id: str,
        dn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> pgm.DistributionPGM:
        """Build and return sampling distribution object

        Args:
            dn_id (str): Name of the sampling distribution to create
            dn_param_dict (dict): Dictionary containing distribution
                parameter names (str) as keys and lists (of either
                strings or StochasticNodePGMs) as values

        Returns:
            DistributionPGM: a DistributionPGM instance
        """

        extracted_val: ty.List[str]

        #############################
        #  Parametric distributions #
        #############################

        if dn_id == "lognormal":
            #############################
            # IMPORTANT: Default values #
            #############################
            ln_n_draws: int = 1
            ln_n_repl: int = 1
            ln_mean: ty.List[float] = []
            ln_sd: ty.List[float] = []
            ln_log_space: bool = True
            parent_node_tracker = dict()

            # { mean: node_pgm1_name, sd: node_pgm2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    # val = val[0] # TODO: deal with vectorization later

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_name

                    # if element in val is string, it remains unchanged,
                    # if NodePGM, we get its string-fied value
                    extracted_val = pgm.extract_value_from_nodepgm(val)

                    if not cls.grammar_check("lognormal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(dnpar.DnLogNormal.DN_NAME, arg)
                        # only one element always
                        ln_n_draws = int(extracted_val[0])

                    elif arg == "nr":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(dnpar.DnLogNormal.DN_NAME, arg)
                        # only one element always
                        ln_n_repl = int(extracted_val[0])

                    elif arg == "mean":
                        # pars[2] = [float(v) for v in _extracted_val]
                        ln_mean = [float(v) for v in extracted_val]

                    elif arg == "sd":
                        # pars[3] = [float(v) for v in _extracted_val]
                        ln_sd = [float(v) for v in extracted_val]

                    elif arg == "log_space":
                        if extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[4] = True
                            ln_log_space = True

                        elif extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[4] = False
                            ln_log_space = False

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in ((ln_mean, "mean"), (ln_sd, "sd")):
                if not par_obj:
                    raise ec.MissingParameterError(par_name)

            # return dnpar.DnLogNormal(pars, parent_node_tracker)
            return dnpar.DnLogNormal(ln_n_draws, ln_n_repl, ln_mean, ln_sd, ln_log_space, parent_node_tracker)

        elif dn_id == "normal":
            #############################
            # IMPORTANT: Default values #
            #############################
            norm_n_draws: int = 1
            norm_n_repl: int = 1
            norm_mean: ty.List[float] = []
            norm_sd: ty.List[float] = []
            parent_node_tracker = dict()

            # { mean: node_pgm1_name, sd: node_pgm2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val = pgm.extract_value_from_nodepgm(val)

                    if not cls.grammar_check("normal", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnNormal.DN_NAME, arg)

                        norm_n_draws = int(extracted_val[0])

                    elif arg == "nr":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnNormal.DN_NAME, arg)

                        try:
                            norm_n_repl = int(extracted_val[0])

                        except ValueError:
                            raise ec.RequireIntegerError(
                                dnpar.DnNormal.DN_NAME, arg)

                    elif arg == "mean":
                        # pars[2] = [float(v) for v in extracted_val]
                        norm_mean = [float(v) for v in extracted_val]

                    elif arg == "sd":
                        # pars[3] = [float(v) for v in extracted_val]
                        norm_sd = [float(v) for v in extracted_val]

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in ((norm_mean, "mean"), (norm_sd, "sd")):
                if not par_obj:
                    raise ec.MissingParameterError(par_name)

            # return dnpar.DnNormal(pars, parent_node_tracker)
            return dnpar.DnNormal(norm_n_draws, norm_n_repl, norm_mean, norm_sd, parent_node_tracker)

        elif dn_id == "exponential":
            #############################
            # IMPORTANT: Default values #
            #############################
            exp_n_draws: int = 1
            exp_n_repl: int = 1
            exp_scale_or_rate: ty.List[float] = []
            exp_rate_parameterization: bool = True
            parent_node_tracker = dict()  # { lambda: node_pgm1_name }

            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        # needed for building inference specifications
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val = pgm.extract_value_from_nodepgm(val)

                    if not cls.grammar_check("exponential", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnExponential.DN_NAME, arg)
                        # pars[0] = int(_extracted_val[0])
                        exp_n_draws = int(extracted_val[0])

                    elif arg == "nr":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnExponential.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        exp_n_repl = int(extracted_val[0])

                    elif arg == "rate":
                        # pars[2] = [float(v) for v in _extracted_val]
                        exp_scale_or_rate = [float(v) for v in extracted_val]

                    elif arg == "rate_parameterization":
                        if extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[3] = True
                            exp_rate_parameterization = True

                        elif extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[3] = False
                            exp_rate_parameterization = False

            # making sure essential parameters of distribution have been specified
            if not exp_scale_or_rate:
                raise ec.MissingParameterError("rate")

            # return dnpar.DnExponential(pars, parent_node_tracker)
            return dnpar.DnExponential(
                exp_n_draws,
                exp_n_repl,
                exp_scale_or_rate,
                exp_rate_parameterization,
                parent_node_tracker)

        elif dn_id == "gamma":
            #############################
            # IMPORTANT: Default values #
            #############################
            gamma_n_draws: int = 1
            gamma_n_repl: int = 1
            gamma_shape: ty.List[float] = []
            gamma_scale_or_rate: ty.List[float] = []
            gamma_rate_parameterization: bool = False
            parent_node_tracker = dict()

            # { mean: node_pgm1_name, sd: node_pgm2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():

                    # needed for building inference specifications
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val = pgm.extract_value_from_nodepgm(val)

                    if not cls.grammar_check("gamma", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnGamma.DN_NAME, arg)
                        # pars[0] = int(_extracted_val[0])
                        gamma_n_draws = int(extracted_val[0])

                    elif arg == "nr":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnGamma.DN_NAME, arg)
                        # pars[1] = int(_extracted_val[0])
                        gamma_n_repl = int(extracted_val[0])

                    elif arg == "shape":
                        # pars[2] = [float(v) for v in _extracted_val]
                        gamma_shape = [float(v) for v in extracted_val]

                    elif arg == "scale":
                        # pars[3] = [float(v) for v in _extracted_val]
                        gamma_scale_or_rate = [float(v) for v in extracted_val]

                    elif arg == "rate_parameterization":
                        if extracted_val[0] in ("\"true\"", "\"T\"", "\"True\""):
                            # pars[3] = True
                            gamma_rate_parameterization = True

                        elif extracted_val[0] in ("\"false\"", "\"F\"", "\"False\""):
                            # pars[3] = False
                            gamma_rate_parameterization = False

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in \
                    ((gamma_shape, "shape"), (gamma_scale_or_rate, "scale")):
                if not par_obj:
                    raise ec.MissingParameterError(par_name)

            # return dnpar.DnGamma(pars, parent_node_tracker)
            return dnpar.DnGamma(
                gamma_n_draws,
                gamma_n_repl,
                gamma_shape,
                gamma_scale_or_rate,
                gamma_rate_parameterization,
                parent_node_tracker)

        elif dn_id == "unif":
            #############################
            # IMPORTANT: Default values #
            #############################
            unif_n_draws: int = 1
            unif_n_repl: int = 1
            unif_min: ty.List[float] = []
            unif_max: ty.List[float] = []
            parent_node_tracker = dict()

            # { mean: node_pgm1_name, sd: node_pgm2_name, ... }
            if dn_param_dict:
                # val is list
                for arg, val in dn_param_dict.items():
                    if isinstance(val[0], pgm.StochasticNodePGM):
                        # needed for building inference specifications
                        parent_node_tracker[arg] = val[0].node_name

                    extracted_val = pgm.extract_value_from_nodepgm(val)

                    if not cls.grammar_check("unif", arg):
                        raise ec.NotAParameterError(arg)

                    elif arg == "n":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnUnif.DN_NAME, arg)
                        unif_n_draws = int(extracted_val[0])

                    elif arg == "nr":
                        if len(extracted_val) > 1:
                            raise ec.RequireSingleValueError(
                                dnpar.DnUnif.DN_NAME, arg)
                        unif_n_repl = int(extracted_val[0])

                    elif arg == "min":
                        unif_min = [float(v) for v in extracted_val]

                    elif arg == "max":
                        unif_max = [float(v) for v in extracted_val]

            # making sure essential parameters of distribution have been specified
            for par_obj, par_name in ((unif_min, "min"), (unif_max, "max")):
                if not par_obj:
                    raise ec.MissingParameterError(par_name)

            # return dnpar.DnUnif(pars, parent_node_tracker)
            return dnpar.DnUnif(unif_n_draws, unif_n_repl, unif_min, unif_max, parent_node_tracker)

        #################################
        #  Non-parametric distributions #
        #################################

        # discrete_sse
        # elif dn_id == "discrete_sse": (switch to this one when there's a new one)
        else:
            # dn_param_dict is validated inside
            return cls.init_return_discrete_SSE_dn(dn_param_dict)
