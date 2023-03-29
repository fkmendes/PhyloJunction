import typing as ty
from _typeshed import Incomplete

MacroevolEvent: Incomplete

class DiscreteStateDependentRate:
    value: ty.List[float]
    name: Incomplete
    state_tuple: Incomplete
    departing_state: Incomplete
    arriving_states: Incomplete
    event: Incomplete
    str_representation: str
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], event: MacroevolEvent, name: str = ..., states: ty.List[int] = ...) -> None: ...
    def sample(self) -> None: ...
    def get_length(self): ...
    def get_gcf(self) -> None: ...

class DiscreteStateDependentProbability:
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], name: str="", state: int=0) -> None: ...
    def _initialize_str_representation(self) -> None: ...

class DiscreteStateDependentParameter:
    val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]]
    name: str
    state: int
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], name: str="", state: int=0) -> None: ...
    def _initialize_str_representation(self) -> None: ...

class DiscreteStateDependentParameterManager:
    matrix_state_dep_params: ty.List[ty.List[DiscreteStateDependentParameter]]
    state_count: int
    seed_age_for_time_slicing: ty.Optional[float]
    n_time_slices: int
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    state_dep_params_dict: Incomplete
    epsilon: float
    def __init__(self, matrix_atomic_rate_params: ty.List[ty.List[DiscreteStateDependentParameter]], total_state_count: int, seed_age_for_time_slicing: ty.Optional[float] = ..., list_time_slice_age_ends: ty.Optional[ty.List[float]] = ..., epsilon: float = ...) -> None: ...
    def init_matrix_state_dep_params_dict(self) -> None: ...
    def state_dep_params_at_time(self, params_matrix, a_time: float): ...

class MacroevolEventHandler:
    state_dep_par_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: float
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    str_representation: str
    def __init__(self, a_fig_rates_manager: DiscreteStateDependentParameterManager) -> None: ...
    def total_rate(self, a_time: float, state_representation_dict: ty.Dict[int, ty.Set[str]], value_idx: int = ..., departing_state: ty.Optional[int] = ..., debug: bool = ...): ...
    def sample_event_atomic_parameter(self, denominator: float, a_time: float, state_indices: ty.List[int], value_idx: int = ..., a_seed: ty.Optional[float] = ..., debug: bool = ...): ...
    def get_gcf(self) -> None: ...
    def get_length(self) -> None: ...
