import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/distribution/)
import typing as ty
import math
import numpy as np
from scipy.stats import expon, lognorm, norm, gamma, uniform # type: ignore

# from tpsimulator_classes import *
# from tpsimulator_utils import *

# pj imports
import pgm.pgm as pgm
import calculation.math_utils as mu
import utility.exception_classes as ec
import utility.helper_functions as pjh

class DnLogNormal(pgm.DistributionPGM):
    
    DN_NAME = "Log-normal"

    n_draws: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.List[ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

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
    def __init__(self, pars: ty.List[ty.Union[int, ty.List[float], bool]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        if isinstance(pars[0], int) and isinstance(pars[1], int):
            self.n_draws = pars[0]
            self.n_repl = pars[1]

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        self.vectorized_params = self.check_sample_size(pars[2:4])

        # params
        self.param_dict["mean_param"] = self.vectorized_params[0]
        self.param_dict["sd_param"] = self.vectorized_params[1]
        self.param_dict["log_space"] = pars[4]

        # for inference, we need to keep track of parent node names
        if isinstance(parent_node_tracker, pgm.DeterministicNodePGM):
            raise ec.InvalidFunctionArgError(self.DN_NAME, "One of the arguments is a deterministic node. This is not allowed. Ignoring command...")

        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()

        try:
            mean_params = self.param_dict["mean_param"]
            sd_params = self.param_dict["sd_param"]
            log_space = self.param_dict["log_space"]

            # use single provided mean and sd for all (n_draws * n_repl) draws
            if len(mean_params) == 1 and len(sd_params) == 1:
                return DnLogNormal.draw_ln(self.n_draws * self.n_repl, float(mean_params[0]), float(sd_params[0]), log_space=log_space).tolist()

            # (DnLogNormal sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(mean_params)):
                    repl = DnLogNormal.draw_ln(self.n_repl, float(mean_params[i]), float(sd_params[i]), log_space=log_space).tolist()
                    sampled_values += repl

                return sampled_values

        except:
            exit("Drawing from log-normal distribution failed." + \
                " A \'mean_param\', \'sd_param\', and \'log_space\' must be specified. Exiting...")

        # return mu.draw_dist(self.n_draws, ParametricDistribution.LOGNORMAL, self.param_dict, n_repl=self.n_repl)

    def check_sample_size(self, param_list) -> ty.List[ty.List[ty.Union[int, float, str]]]:
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        log_space = self.param_dict["log_space"]
        real_mean_list = self.param_dict["mean_param"] # one per simulation
        real_sd_list = self.param_dict["sd_param"] # one per simulation

        if log_space:
            real_mean_list = [math.exp(i) for i in real_mean_list]
            real_sd_list = [math.exp(i) for i in real_sd_list]

        # real_mean_list and real_sd_list will have n_sim values inside, even if they are all the same
        for ith_sim in range(self.n_draws):
            ith_sim_str = "dnLognormal(mean="

            # if we can find a node that holds the value of the mean, we use it
            try:
                ith_sim_str += self.parent_node_tracker["mean"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_mean_list[ith_sim])

            ith_sim_str += ", sd="

            # if we can find a node that holds the value of the sd, we use it
            try:
                ith_sim_str += self.parent_node_tracker["sd"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_sd_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

        return rev_str_list

##############################################################################

class DnNormal(pgm.DistributionPGM):

    DN_NAME = "Normal"

    n_draws: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.List[ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]

    @staticmethod
    def draw_normal(n_draws, mean_param: float, sd_param: float) -> ty.Union[np.float64, np.ndarray]:
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
    def __init__(self, pars: ty.List[ty.Union[int, ty.List[float]]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        if isinstance(pars[0], int) and isinstance(pars[1], int):
            self.n_draws = pars[0]
            self.n_repl = pars[1]

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        self.vectorized_params = self.check_sample_size(pars[2:4]) # one element per parameter

        # params
        self.param_dict["mean_param"] = self.vectorized_params[0]
        self.param_dict["sd_param"] = self.vectorized_params[1]

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            mean_params = self.param_dict["mean_param"]
            sd_params = self.param_dict["sd_param"]

            # use single provided mean and sd for all (n_draws * n_repl) draws
            if len(mean_params) == 1 and len(sd_params) == 1:
                return DnNormal.draw_normal(self.n_draws * self.n_repl, float(mean_params[0]), sd_params[0]).tolist()

            # (DnNormal sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(mean_params)):
                    repl = DnNormal.draw_normal(self.n_repl, mean_params[i], sd_params[i]).tolist()
                    sampled_values += repl

                return sampled_values

        except:
            exit("Drawing from normal distribution failed." + \
                " A \'mean_param\' and \'sd_param\' must be specified. Exiting...")

        # return mu.draw_dist(self.n_draws, ParametricDistribution.NORMAL, self.param_dict, n_repl=self.n_repl)

    def check_sample_size(self, param_list):
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        real_mean_list = self.param_dict["mean_param"] # one per simulation
        real_sd_list = self.param_dict["sd_param"] # one per simulation

        for ith_sim in range(self.n_draws):
            ith_sim_str = "dnNormal(mean="

            # if we can find a node that holds the value of the mean, we use it
            try:
                ith_sim_str += self.parent_node_tracker["mean"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_mean_list[ith_sim])

            ith_sim_str += ", sd="

            # if we can find a node that holds the value of the sd, we use it
            try:
                ith_sim_str += self.parent_node_tracker["sd"] # returns NodePGM, and we grab its name
            except:
                ith_sim_str += str(real_sd_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

        return rev_str_list

##############################################################################

class DnExponential(pgm.DistributionPGM):

    DN_NAME = "Exponential"

    n_draws: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.List[ty.List[ty.Union[int, float, str]]]]
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
    def __init__(self, pars: ty.List[ty.Union[int, ty.List[float], bool]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        if isinstance(pars[0], int) and isinstance(pars[1], int):
            self.n_draws = pars[0]
            self.n_repl = pars[1]

        # makes sure these parameters are lists, or converts them into list
        # also multiplying element (when only one is provided) if necessary
        self.vectorized_params = self.check_sample_size(pars[2])

        # params
        self.param_dict["scale_or_rate_param"] = self.vectorized_params[0]
        self.param_dict["rate_parameterization"] = pars[3] # True is the default

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            scale_or_rate_params = self.param_dict["scale_or_rate_param"]
            rate_parameterization = self.param_dict["rate_parameterization"]

            # use single provided scale_or_rate_params for all (n_draws * n_repl) draws
            if len(scale_or_rate_params) == 1:
                # so mypy won't complain
                if isinstance(scale_or_rate_params[0], float) and isinstance(rate_parameterization, bool):
                    return DnExponential.draw_exp(self.n_draws * self.n_repl, scale_or_rate_params[0], rate_parameterization=rate_parameterization).tolist()

            # (DnExponential sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(scale_or_rate_params)):
                    # so mypy won't complain
                    if isinstance(scale_or_rate_params[i], float) and isinstance(rate_parameterization, bool):
                        repl = DnExponential.draw_exp(self.n_repl, scale_or_rate_params[i], rate_parameterization=rate_parameterization).tolist()
                        sampled_values += repl

                return sampled_values

        except:
            exit("Drawing from exponential distribution failed." + \
                " A \'scale_or_rate_param\' and \'rate_parameterization\' must be specified. Exiting...")

        # return mu.draw_dist(self.n_draws, ParametricDistribution.EXPONENTIAL, self.param_dict, n_repl=self.n_repl)

    def check_sample_size(self, param_list):
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        scale_or_rate_list = self.param_dict["scale_or_rate_param"]
        rate_parameterization = self.param_dict["rate_parameterization"] # True is the default in PJ

        # if user wants scale parameterization in PJ, we need to invert it
        # for use with rev
        if not rate_parameterization:
            scale_or_rate_list = [1.0/v for v in scale_or_rate_list]

        for ith_sim in range(self.n_draws):
            ith_sim_str = "dnExponential(lambda="

            # if we can find a node that holds the value of the rate, we use it
            try:
                ith_sim_str += self.parent_node_tracker["rate"] # key: arg in PJ syntax, value: NodePGM name passed as arg
            except:
                ith_sim_str += str(scale_or_rate_list[ith_sim])
            ith_sim_str += ")"

            rev_str_list.append(ith_sim_str)

        return rev_str_list

##############################################################################

class DnGamma(pgm.DistributionPGM):

    DN_NAME = "Gamma"

    n_draws: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.List[ty.List[ty.Union[int, float, str]]]]
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
    def __init__(self, pars: ty.List[ty.Union[int, ty.List[float], bool]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        if isinstance(pars[0], int) and isinstance(pars[1], int):
            self.n_draws = pars[0]
            self.n_repl = pars[1]

        # makes sure these parameters are lists, or converts them into lists
        # also multiplying element (when only one is provided) if necessary
        self.vectorized_params = self.check_sample_size(pars[2:4])

        # params
        self.param_dict["shape_param"] = self.vectorized_params[0]
        self.param_dict["scale_or_rate_param"] = self.vectorized_params[1]
        self.param_dict["rate_parameterization"] = pars[4] # False is the default

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = list()
        
        try:
            shape_params = self.param_dict["shape_param"]
            scale_or_rate_params = self.param_dict["scale_or_rate_param"]
            rate_parameterization = self.param_dict["rate_parameterization"]

            # use single provided mean and sd for all (n_draws * n_repl) draws
            if len(shape_params) == 1 and len(scale_or_rate_params) == 1:
                # so mypy won't complain
                if isinstance(shape_params[0], float) and isinstance(scale_or_rate_params[0], float) and isinstance(rate_parameterization, bool):
                    return DnGamma.draw_gamma(self.n_draws * self.n_repl, shape_params[0], scale_or_rate_params[0], rate_parameterization=rate_parameterization).tolist()

            # (DnGamma should have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(shape_params)):
                    # so mypy won't complain
                    if isinstance(shape_params[i], float) and isinstance(scale_or_rate_params[i], float) and isinstance(rate_parameterization, bool):
                        repl = DnGamma.draw_gamma(self.n_repl, shape_params[i], scale_or_rate_params[i], rate_parameterization=rate_parameterization).tolist()
                        sampled_values += repl

                return sampled_values

        except:
            exit("Drawing from gamma distribution failed." + \
                " A \'shape_param\' and \'scale_or_rate_param\' must be specified. Exiting...")
        
        # return mu.draw_dist(self.n_draws, ParametricDistribution.GAMMA, self.param_dict, n_repl=self.n_repl)

    def check_sample_size(self, param_list):
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)

##############################################################################

class DnUnif(pgm.DistributionPGM):

    DN_NAME = "Uniform"

    n_draws: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.List[ty.List[ty.Union[int, float, str]]]]
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
    def __init__(self, pars: ty.List[ty.Union[int, ty.List[float]]], parent_node_tracker: ty.Optional[ty.Dict[str, str]]=None) -> None:
        self.param_dict = dict()

        # so mypy won't complain...
        if isinstance(pars[0], int) and isinstance(pars[1], int):
            self.n_draws = pars[0]
            self.n_repl = pars[1]

        # makes sure these parameters are lists, or converts them into lists
        # also multiplying element (when only one is provided) if necessary
        self.vectorized_params = self.check_sample_size(pars[2:4])

        # params
        self.param_dict["min_param"] = self.vectorized_params[0]
        self.param_dict["max_param"] = self.vectorized_params[1]

        # for inference, we need to keep track of parent node names
        self.parent_node_tracker = parent_node_tracker

    def generate(self) -> ty.List[float]:
        sampled_values: ty.List[float] = []

        try:
            min_params = self.param_dict["min_param"]
            max_params = self.param_dict["max_param"]

            # use single provided mean and sd for all (n_draws * n_repl) draws
            if len(min_params) == 1 and len(max_params) == 1:
                return DnUnif.draw_unif(self.n_draws * self.n_repl, float(min_params[0]), float(max_params[0])).tolist()

            # (PJUniform sould have taken care of checking dimensions)
            # otherwise, use one mean and sd for each n_repl draws for each simulation
            else:
                for i in range(len(min_params)):
                    repl = DnUnif.draw_unif(self.n_repl, float(min_params[i]), float(max_params[i])).tolist()
                    sampled_values += repl # [0] because the output comes in a numpy nd.array

                return sampled_values

        except:
            exit("Drawing from uniform distribution failed." + \
                " A \'min_param\' and \'max_param\' must be specified. Exiting...")

        # print("inside PJUnif, generate")
        # print("n_draws = " + str(self.n_draws))
        # print("n_repl = " + str(self.n_repl))
        # print("param_dict = ")
        # print(self.param_dict)

        # return mu.draw_dist(self.n_draws, ParametricDistribution.UNIFORM, self.param_dict, n_repl=self.n_repl)

    def check_sample_size(self, param_list):
        return pjh.verify_or_convert2_vector(param_list, self.DN_NAME, size_to_grow=self.n_draws)

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        min_list = self.param_dict["min_param"]
        max_list = self.param_dict["max_param"]

        for ith_sim in range(self.n_draws):
            rev_str_list.append("dnUniform(lower=" + str(min_list[ith_sim]) + ", upper=" + str(max_list[ith_sim]) + ")")

        return rev_str_list

##############################################################################

if __name__ == "__main__":
    # can be called from distribution/
    # $ python3 dn_parametric.py
    # 
    # can also be called from phylojunction/
    # $ python3 distribution/dn_parametric.py
    # or
    # $ python3 -m distribution.dn_parametric.py
    #
    # can also be called from VS Code, if open folder is phylojuction/

    # ln1 = DnLogNormal([1, [10,9,8], [7,6,5], True]) # should throw exception
    n_draws = 3
    n_repl = 1
    ln2 = DnLogNormal([n_draws, n_repl, [-2, -2.1, -2.2], [0.1, 0.2, 0.15], True])
    ln3 = DnLogNormal([n_draws, n_repl, [-2], [0.1], True])
    ln4 = DnLogNormal([n_draws, n_repl, [-2], [0.1], True])

    print(ln2.generate())
    print(ln3.generate())
    print(ln4.generate())