import numpy as np
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm

class DnLogNormal(pgm.DistributionPGM):
    DN_NAME: str
    n_samples: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]
    @staticmethod
    def draw_ln(n_samples: int, mean_param: float, sd_param: float, scale: float = ..., log_space: bool = ...) -> ty.Union[np.float64, np.ndarray]: ...
    def __init__(self, n_samples: int, n_repl: int, ln_mean_param: ty.List[float], ln_sd_param: ty.List[float], ln_log_space: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) -> None: ...
    def generate(self) -> ty.List[float]: ...
    def init_check_vectorize_sample_size(self, param_list: ty.Optional[ty.List[ty.Any]]=None) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...

class DnNormal(pgm.DistributionPGM):
    DN_NAME: str
    n_samples: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]
    @staticmethod
    def draw_normal(n_samples: int, mean_param: float, sd_param: float) -> ty.Union[np.float64, np.ndarray]: ...
    def __init__(self, n_samples: int, n_repl: int, norm_mean_param: ty.List[float], sd_param: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) -> None: ...
    def generate(self) -> ty.List[float]: ...
    def init_check_vectorize_sample_size(self, param_list: ty.Optional[ty.List[ty.Any]]=None) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...

class DnExponential(pgm.DistributionPGM):
    DN_NAME: str
    n_samples: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]
    @staticmethod
    def draw_exp(n_samples: int, scale_or_rate_param: float, rate_parameterization: bool = ...) -> ty.Union[np.float64, np.ndarray]: ...
    def __init__(self, n_samples: int, n_repl: int, rate_param: ty.List[float], rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) -> None: ...
    def generate(self) -> ty.List[float]: ...
    def init_check_vectorize_sample_size(self, param_list: ty.Optional[ty.List[ty.Any]]=None) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...

class DnGamma(pgm.DistributionPGM):
    DN_NAME: str
    n_samples: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]
    @staticmethod
    def draw_gamma(n_samples: int, shape_param: float, scale_or_rate_param: float, rate_parameterization: bool = ...) -> ty.Union[np.float64, np.ndarray]: ...
    def __init__(self, n_samples: int, n_repl: int, shape_param: ty.List[float], scale_or_rate_param: ty.List[float], rate_parameterization: bool, parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) -> None: ...
    def generate(self) -> ty.List[float]: ...
    def init_check_vectorize_sample_size(self, param_list: ty.Optional[ty.List[ty.Any]]=None) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...

class DnUnif(pgm.DistributionPGM):
    DN_NAME: str
    n_samples: int
    n_repl: int
    vectorized_params: ty.List[ty.List[ty.Union[int, float, str]]]
    param_dict: ty.Dict[str, ty.Union[bool, ty.List[ty.Union[int, float, str]]]]
    parent_node_tracker: ty.Optional[ty.Dict[str, str]]
    @staticmethod
    def draw_unif(n_samples: int, min_param: float, max_param: float) -> ty.Union[np.float64, np.ndarray]: ...
    def __init__(self, n_samples: int, n_repl: int, min_param: ty.List[float], max_param: ty.List[float], parent_node_tracker: ty.Optional[ty.Dict[str, str]] = None) -> None: ...
    def generate(self) -> ty.List[float]: ...
    def init_check_vectorize_sample_size(self, param_list: ty.Optional[ty.List[ty.Any]]=None) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
