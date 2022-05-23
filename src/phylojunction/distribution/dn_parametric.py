import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/distribution/)
import typing as ty
import math
import numpy as np
from scipy.stats import expon, lognorm, norm, gamma, uniform # type: ignore

# pj imports
import pgm.pgm as pgm
import calculation.math_utils as mu
import utility.exception_classes as ec
import utility.helper_functions as pjh
import inference.revbayes.rb_dn_parametric as rbpar

class DnLogNormal(pgm.DistributionPGM):
    
    DN_NAME = "Log-normal"

    n_draws: int
    n_repl: int
    ln_mean_list: ty.List[float]
    ln_sd_list: ty.List[float]
    ln_log_space: bool
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]] = dict()

    @staticmethod
    def draw_ln(n_draws: int, mean_param: float, sd_param: float, scale: float=1.0, log_space: bool=True) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from log-normal distribution.

        Args:
            n_draws (int): Number of draws (sample size).
            mean_param (float): Mean (location) of log-normal in log-space.
            sd_param (float): Std. deviation (shape) of log-normal distribution in log-space.
            scale (float): Median (scale) of log-normal distribution in log-space. Defaults to 1.0.
            log_space (bool, optional): Mean and std. deviation are provided in natural space. Defaults to True.

        Returns:
            list of float(s): Sample (list) from log-normal distribution.
        """

        if not log_space:
            sd_param = math.log(sd_param)
            mean_param = math.log(mean_param)

        # debugging
        # print("log-normal in log-space: mean = " + str(mean_param) + " sd = " + str(sd_param))

        return lognorm.rvs(s=math.exp(sd_param), loc=math.exp(mean_param), scale=1.0, size=n_draws)

    # validation of pars happens in phylojunction_grammar
    def __init__(self, n_draws: int, n_repl: int, ln_mean: ty.List[float], ln_sd: ty.List[float], ln_log_space: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:
        self.param_dict = dict()
        self.n_draws = n_draws
        self.n_repl = n_repl
        self.ln_mean_arg = ln_mean
        self.ln_sd_arg = ln_sd
        self.ln_log_space = ln_log_space

        check_sample_size_return = self.check_sample_size([self.ln_mean_arg, self.ln_sd_arg]) # one element per parameter
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        self.ln_mean_list = ty.cast(ty.List[float], self.vectorized_params[0])
        self.ln_sd_list = ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        if isinstance(parent_node_tracker, pgm.DeterministicNodePGM):
            raise ec.InvalidFunctionArgError(self.DN_NAME, "One of the arguments is a deterministic node. This is not allowed. Ignoring command...")

        self.parent_node_tracker = parent_node_tracker


    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()

        try:
            if len(self.ln_mean_list) == 1 and len(self.ln_sd_list) == 1:
                return ty.cast(ty.List[float], DnLogNormal.draw_ln(self.n_draws * self.n_repl, self.ln_mean_list[0], self.ln_sd_list[0], log_space=self.ln_log_space).tolist())

            # (DnLogNormal sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(self.ln_mean_list)):
                    # so mypy won't complain             
                    repl = ty.cast(ty.List, DnLogNormal.draw_ln(self.n_repl, self.ln_mean_list[i], self.ln_sd_list[i], log_space=self.ln_log_space).tolist())
                    sampled_values += repl

                return sampled_values

        except:
            raise ec.DnInitMisspec(self.DN_NAME, "A \'mean_param\', \'sd_param\', and \'log_space\' must be specified.")


    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)


    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return rbpar.get_ln_rev_inference_spec_info(self.n_draws, self.ln_mean_list, self.ln_sd_list, self.parent_node_tracker)

##############################################################################

class DnNormal(pgm.DistributionPGM):

    DN_NAME = "Normal"

    n_draws: int
    n_repl: int
    norm_mean_param_list: ty.List[float]
    norm_sd_param_list: ty.List[float]
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_normal(n_draws: int, mean_param: float, sd_param: float) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from normal distribution

        Args:
            n_draws (int): Number of draws (sample size)
            mean_param (float): Mean (location) of normal distribution
            sd_param (float): Std. deviation (scale) of normal distribution

        Returns:
            list of floats(s): Sample (list) from normal distribution
        """

        # debugging
        # print("normal: mean = " + str(mean_param) + " sd = " + str(sd_param))

        return norm.rvs(scale=sd_param, loc=mean_param, size=n_draws)

    # validation of pars happens in phylojunction_grammar
    # def __init__(self, pars: ty.List[ty.Union[int, ty.List[float]]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
    def __init__(self, n_draws: int, n_repl: int, norm_mean_param: ty.List[float], norm_sd_param: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        # if isinstance(pars[0], int) and isinstance(pars[1], int):
        #     self.n_draws = pars[0]
        #     self.n_repl = pars[1]
        self.n_draws = n_draws
        self.n_repl = n_repl
        self.norm_mean_param_arg = norm_mean_param
        self.norm_sd_param_arg = norm_sd_param

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        # check_sample_size_return = self.check_sample_size(pars[2:4]) # one element per parameter
        check_sample_size_return = self.check_sample_size([self.norm_mean_param_arg, self.norm_sd_param_arg]) # one element per parameter
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["mean_param"] = self.vectorized_params[0]
        # self.param_dict["sd_param"] = self.vectorized_params[1]
        self.norm_mean_param_list = ty.cast(ty.List[float], self.vectorized_params[0])
        self.norm_sd_param_list = ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker


    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            # mean_params = self.param_dict["mean_param"]
            # sd_params = self.param_dict["sd_param"]

            # so mypy won't complain
            # if isinstance(mean_params, list) and isinstance(sd_params, list):

            # use single provided mean and sd for all (n_draws * n_repl) draws
            # if len(mean_params) == 1 and len(sd_params) == 1:
            if len(self.norm_mean_param_list) == 1 and len(self.norm_sd_param_list) == 1:
                # so mypy won't complain
                # if isinstance(mean_params, float) and isinstance(sd_params, float):
                # return ty.cast(ty.List, DnNormal.draw_normal(self.n_draws * self.n_repl, mean_params[0], sd_params[0]).tolist())
                return ty.cast(ty.List, DnNormal.draw_normal(self.n_draws * self.n_repl, self.norm_mean_param_list[0], self.norm_sd_param_list[0]).tolist())


            # (DnNormal sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                # for i in range(len(mean_params)):
                for i in range(len(self.norm_mean_param_list)):
                    # so mypy won't complain             
                    # for some reason, just checking type like below won't work here
                    # ... no clue why (I'm still calling isinstance() just to be safe)
                    # mypy is only satisfied if I forcefully cast into float
                    # if isinstance(mean_params[i], (int, float)) and isinstance(sd_params[i], (int, float)):
                    #     mean_param = ty.cast(float, mean_params[i])
                    #     sd_param = ty.cast(float, sd_params[i])
                        # repl = ty.cast(ty.List, DnNormal.draw_normal(self.n_repl, mean_param, sd_param).tolist())
                    repl = ty.cast(ty.List, DnNormal.draw_normal(self.n_repl, self.norm_mean_param_list[i], self.norm_sd_param_list[i]).tolist())
                    sampled_values += repl

                return sampled_values

        except:
            raise ec.DnInitMisspec(self.DN_NAME, "A \'mean_param\' and \'sd_param\' must be specified.")


    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)


    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return rbpar.get_normal_rev_inference_spec_info(self.n_draws, self.norm_mean_param_list, self.norm_sd_param_list, self.parent_node_tracker)

##############################################################################

class DnExponential(pgm.DistributionPGM):

    DN_NAME = "Exponential"

    n_draws: int
    n_repl: int
    exp_scale_or_rate_list: ty.List[float]
    exp_rate_parameterization: bool = True
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_exp(n_draws: int, scale_or_rate_param: float, rate_parameterization: bool=True) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from exponential distribution.

        Args:
            n_draws (int): Number of draws (sample size).
            scale_or_rate_param (float): Scale (default) or rate of exponential distribution.
            rate_parameterization (bool, optional): Argument of 'scale_or_rate_param' is rate instead of scale. Defaults to True.

        Returns:
            float: Sample (list) from exponential distribution.
        """

        if rate_parameterization:
            scale_or_rate_param = 1.0 / scale_or_rate_param

        return expon.rvs(loc=0.0, scale=scale_or_rate_param, size=n_draws)

    # validation of pars happens in phylojunction_grammar
    # def __init__(self, pars: ty.List[ty.Union[int, ty.List[float], bool]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
    def __init__(self, n_draws: int, n_repl: int, scale_or_rate_param: ty.List[float], rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        # if isinstance(pars[0], int) and isinstance(pars[1], int) and isinstance(pars[3], bool):
        #     self.n_draws = pars[0]
        #     self.n_repl = pars[1]
        #     self.param_dict["rate_parameterization"] = pars[3] # True is the default

        self.n_draws = n_draws
        self.n_repl = n_repl
        self.exp_scale_or_rate_param_arg = scale_or_rate_param
        self.exp_rate_parameterization = rate_parameterization

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        # check_sample_size_return = self.check_sample_size([pars[2]]) # must convert to list because of method signature
        check_sample_size_return = self.check_sample_size([self.exp_scale_or_rate_param_arg]) # must convert to list because of method signature
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["scale_or_rate_param"] = self.vectorized_params[0]
        self.exp_scale_or_rate_param_list = ty.cast(ty.List[float], self.vectorized_params[0])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker


    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            # scale_or_rate_params = self.param_dict["scale_or_rate_param"]
            # rate_parameterization = self.param_dict["rate_parameterization"]

            # so mypy won't complain
            # if isinstance(scale_or_rate_params, list):

            # use single provided scale_or_rate_params for all (n_draws * n_repl) draws
            # if len(scale_or_rate_params) == 1:
            if len(self.exp_scale_or_rate_param_list) == 1:
                
                # so mypy won't complain
                # if isinstance(scale_or_rate_params[0], float) and isinstance(rate_parameterization, bool):
                    # return ty.cast(ty.List[float], DnExponential.draw_exp(self.n_draws * self.n_repl, scale_or_rate_params[0], rate_parameterization=rate_parameterization).tolist())
                    return ty.cast(ty.List[float], DnExponential.draw_exp(self.n_draws * self.n_repl, self.exp_scale_or_rate_param_list[0], rate_parameterization=self.exp_rate_parameterization).tolist())

            # (DnExponential sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                # for i in range(len(scale_or_rate_params)):
                for i in range(len(self.exp_scale_or_rate_param_list)):
                    # so mypy won't complain
                    # for some reason, just checking type like below won't work here
                    # ... no clue why (I'm still calling isinstance() just to be safe)
                    # mypy is only satisfied if I forcefully cast into float
                    # if isinstance(scale_or_rate_params[i], float) and isinstance(rate_parameterization, bool):
                    #     scale_or_rate_param = ty.cast(float, scale_or_rate_params[i])
                    #     rate_parameterization = ty.cast(bool, rate_parameterization)
                        # repl = ty.cast(ty.List[float], DnExponential.draw_exp(self.n_repl, scale_or_rate_param, rate_parameterization=rate_parameterization).tolist())
                        repl = ty.cast(ty.List[float], DnExponential.draw_exp(self.n_repl, self.exp_scale_or_rate_param_list[i], rate_parameterization=self.exp_rate_parameterization).tolist())
                        sampled_values += repl

                return sampled_values

        except:
            raise ec.DnInitMisspec(self.DN_NAME, "A \'scale_or_rate_param\' and \'rate_parameterization\' must be specified.")


    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)


    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return rbpar.get_exponential_rev_inference_spec_info(self.n_draws, self.exp_scale_or_rate_param_list, self.exp_rate_parameterization, self.parent_node_tracker)

##############################################################################

class DnGamma(pgm.DistributionPGM):

    DN_NAME = "Gamma"

    n_draws: int
    n_repl: int
    gamma_shape_param_list: ty.List[float]
    gamma_scale_or_rate_param_list: ty.List[float]
    gamma_rate_parameterization: bool
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]] = dict()
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_gamma(n_draws: int, shape_param: float, scale_or_rate_param: float, rate_parameterization: bool=False) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from gamma distribution

        Args:
            n_draws (int): Number of draws (sample size)
            shape_param (float): Gamma distribution shape parameter (represented by alpha or kappa sometimes)
            scale_or_rate_param (float): Gamma distribution scale or rate parameter.
            rate_parameterization (bool, optional): Argument of 'scale_or_rate_param' is rate instead of scale. Defaults to False.

        Returns:
            list of floats(s): Sample (list) from gamma distribution
        """

        if rate_parameterization:
            scale_or_rate_param = 1.0 / scale_or_rate_param

        return gamma.rvs(a=shape_param, scale=scale_or_rate_param, size=n_draws)

    # validation of pars happens in phylojunction_grammar
    # def __init__(self, pars: ty.List[ty.Union[int, ty.List[float], bool]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
    def __init__(self, n_draws: int, n_repl: int, shape_param: ty.List[float], scale_or_rate_param: ty.List[float], rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        # if isinstance(pars[0], int) and isinstance(pars[1], int):
        #     self.n_draws = pars[0]
        #     self.n_repl = pars[1]

        self.n_draws = n_draws
        self.n_repl = n_repl
        self.gamma_rate_parameterization = rate_parameterization
        self.gamma_shape_param_arg = shape_param
        self.gamma_scale_or_rate_param_arg = scale_or_rate_param

        # makes sure these parameters are lists, or converts them into lists
        # also multiplying element (when only one is provided) if necessary
        # check_sample_size_return = self.check_sample_size(pars[2:4]) # one element per parameter
        check_sample_size_return = self.check_sample_size([self.gamma_shape_param_arg, self.gamma_scale_or_rate_param_arg]) # one element per parameter
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["shape_param"] = self.vectorized_params[0]
        # self.param_dict["scale_or_rate_param"] = self.vectorized_params[1]
        # self.param_dict["rate_parameterization"] = bool(pars[4]) # False is the default
        self.gamma_shape_param_list = ty.cast(ty.List[float], self.vectorized_params[0])
        self.gamma_scale_or_rate_param_list = ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker


    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()
        
        try:
            # shape_params = self.param_dict["shape_param"]
            # scale_or_rate_params = self.param_dict["scale_or_rate_param"]
            # rate_parameterization = self.param_dict["rate_parameterization"]

            # so mypy won't complain
            # if isinstance(shape_params, list) and isinstance(scale_or_rate_params, list):
                # use single provided mean and sd for all (n_draws * n_repl) draws                
                # if len(shape_params) == 1 and len(scale_or_rate_params) == 1:
            if len(self.gamma_shape_param_list) == 1 and len(self.gamma_scale_or_rate_param_list) == 1:
                
                # so mypy won't complain
                # if isinstance(shape_params[0], float) and isinstance(scale_or_rate_params[0], float) and isinstance(rate_parameterization, bool):
                #    return ty.cast(ty.List[float], DnGamma.draw_gamma(self.n_draws * self.n_repl, shape_params[0], scale_or_rate_params[0], rate_parameterization=rate_parameterization).tolist())
                return ty.cast(ty.List[float], DnGamma.draw_gamma(self.n_draws * self.n_repl, self.gamma_shape_param_list[0], self.gamma_scale_or_rate_param_list[0], rate_parameterization=self.gamma_rate_parameterization).tolist())

            # (DnGamma should have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                # for i in range(len(shape_params)):
                for i in range(len(self.gamma_shape_param_list)):
                    # so mypy won't complain
                    # for some reason, just checking type like below won't work here
                    # ... no clue why (I'm still calling isinstance() just to be safe)
                    # mypy is only satisfied if I forcefully cast into float
                    # if isinstance(shape_params[i], float) and isinstance(scale_or_rate_params[i], float) and isinstance(rate_parameterization, bool):
                    #     shape_param = ty.cast(float, shape_params[i])
                    #     scale_or_rate_param = ty.cast(float, scale_or_rate_params[i])
                    #     rate_parameterization = ty.cast(bool, rate_parameterization)
                    #     repl = ty.cast(ty.List[float], DnGamma.draw_gamma(self.n_repl, shape_param, scale_or_rate_param, rate_parameterization=rate_parameterization).tolist())
                    repl = ty.cast(ty.List[float], DnGamma.draw_gamma(self.n_repl, self.gamma_shape_param_list[i], self.gamma_scale_or_rate_param_list[i], rate_parameterization=self.gamma_rate_parameterization).tolist())
                    sampled_values += repl

                return sampled_values

        except:
            raise ec.DnInitMisspec(self.DN_NAME, "A \'shape_param\' and \'scale_or_rate_param\' must be specified.")


    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)


    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return rbpar.get_gamma_rev_inference_spec_info(self.n_draws, self.gamma_shape_param_list, self.gamma_scale_or_rate_param_list, self.gamma_rate_parameterization, self.parent_node_tracker)

##############################################################################

class DnUnif(pgm.DistributionPGM):

    DN_NAME = "Uniform"

    n_draws: int
    n_repl: int
    min_param_list: ty.List[float]
    max_param_list: ty.List[float]
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_unif(n_draws: int, min_param: float, max_param: float) -> ty.Union[np.float64, np.ndarray]:
        """Return sample from log-normal distribution

        Args:
            n_draws (int): Number of draws (sample size)
            min_param (float): Minimum value
            max_param (float): Maximum value

        Returns:
            list of float(s): Sample (list) from uniform distribution
        """

        # see definition of uniform from scipy (need to subtract min_param)
        # if we want the "usual" behavior of a uniform distribution, where
        # scale behaves as a maximum
        return uniform.rvs(loc=min_param, scale=max_param - min_param, size=n_draws)

    # validation of pars happens in phylojunction_grammar
    # def __init__(self, pars: ty.List[ty.Union[int, ty.List[float]]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
    def __init__(self, n_draws: int, n_repl: int, min_param: ty.List[float], max_param: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]]) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        # if isinstance(pars[0], int) and isinstance(pars[1], int):
        #     self.n_draws = pars[0]
        #     self.n_repl = pars[1]
        self.n_draws = n_draws
        self.n_repl = n_repl
        self.min_param_arg = min_param
        self.max_param_arg = max_param

        # makes sure these parameters are lists, or converts them into lists
        # also multiplying element (when only one is provided) if necessary
        # check_sample_size_return = self.check_sample_size(pars[2:4]) # one element per parameter
        check_sample_size_return = self.check_sample_size([self.min_param_arg, self.max_param_arg]) # one element per parameter
        if isinstance(check_sample_size_return, list):
            self.vectorized_params = check_sample_size_return

        # params
        # self.param_dict["min_param"] = self.vectorized_params[0]
        # self.param_dict["max_param"] = self.vectorized_params[1]
        self.min_param_list = ty.cast(ty.List[float], self.vectorized_params[0])
        self.max_param_list = ty.cast(ty.List[float], self.vectorized_params[1])

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[ty.Union[float, ty.Any]] = []

        try:
            # min_params = self.param_dict["min_param"]
            # max_params = self.param_dict["max_param"]

            # so mypy won't complain
            # if isinstance(min_params, list) and isinstance(max_params, list):
                # use single provided mean and sd for all (n_draws * n_repl) draws
                # if len(min_params) == 1 and len(max_params) == 1:
            if len(self.min_param_list) == 1 and len(self.max_param_list) == 1:
                # so mypy won't complain
                # if isinstance(min_params, float) and isinstance(max_params, float):
                    # return ty.cast(ty.List[float], DnUnif.draw_unif(self.n_draws * self.n_repl, float(min_params[0]), float(max_params[0])).tolist())
                return ty.cast(ty.List[float], DnUnif.draw_unif(self.n_draws * self.n_repl, self.min_param_list[0], self.max_param_list[0]).tolist())

            # (DnUniform sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                # for i in range(len(min_params)):
                for i in range(len(self.min_param_list)):
                    # so mypy won't complain                    
                    # for some reason, just checking type like below won't work here
                    # ... no clue why (I'm still calling isinstance() just to be safe)
                    # mypy is only satisfied if I forcefully cast into float
                    # if isinstance(min_params[i], (int, float)) and isinstance(max_params[i], (int, float)):
                    #     min_param = ty.cast(float, min_params[i])
                    #     max_param = ty.cast(float, max_params[i])
                    #     repl = ty.cast(ty.List[float], DnUnif.draw_unif(self.n_repl, min_param, max_param).tolist())
                        repl = ty.cast(ty.List[float], DnUnif.draw_unif(self.n_repl, self.min_param_list[i], self.max_param_list[i]).tolist())
                        sampled_values += repl

                return sampled_values

        # mypy doesn't understand this except is mandatory and would necessarily
        # quit the program
        except:
            raise ec.DnInitMisspec(self.DN_NAME, "A \'min_param\' and \'max_param\' must be specified. Exiting...")


    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)


    def get_rev_inference_spec_info(self) -> ty.List[str]:
        return rbpar.get_unif_rev_inference_spec_info(self.n_draws, self.min_param_list, self.max_param_list, self.parent_node_tracker)

##############################################################################

if __name__ == "__main__":
    # can be called from distribution/
    # $ python3 dn_parametric.py
    # 
    # can also be called from phylojunction/
    # $ python3 distribution/dn_parametric.py
    # or
    # $ python3 -m distribution.dn_parametric
    #
    # can also be called from VS Code, if open folder is phylojuction/

    # ln1 = DnLogNormal([1, [10,9,8], [7,6,5], True]) # should throw exception
    n_draws = 3
    n_repl = 1

    dummy_parent_node_tracker: ty.Dict[str, str] = dict()

    ln2 = DnLogNormal(n_draws, n_repl, [-2, -2.1, -2.2], [0.1, 0.2, 0.15], True, dummy_parent_node_tracker)
    ln3 = DnLogNormal(n_draws, n_repl, [-2], [0.1], True, dummy_parent_node_tracker)
    ln4 = DnLogNormal(n_draws, n_repl, [-2], [0.1], True, dummy_parent_node_tracker)

    print(ln2.generate())
    print(ln3.generate())
    print(ln4.generate())