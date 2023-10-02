import typing as ty
from _typeshed import Incomplete

MacroevolEvent: Incomplete
DiscreteStateDependentParameterType: Incomplete

MacroevolEventValue = ty.Literal[0, 1, 2, 3, 4, 5]
DiscreteStateDependentParameterTypeValue: ty.Literal[0, 1]

class DiscreteStateDependentRate:
    value: ty.List[float]
    name: Incomplete
    states: ty.List[int]
    state_tuple: ty.Tuple[int]
    departing_state: int
    arriving_states: ty.Union[int, ty.Tuple[int]]
    event: MacroevolEvent
    str_representation: str
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], event: MacroevolEvent, name: str = ..., states: ty.List[int]=[]) -> None: ...
    def sample(self) -> None: ...
    def get_length(self): ...
    def get_gcf(self) -> None: ...

class DiscreteStateDependentProbability:
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], name: str="", state: int=0) -> None: ...
    def _initialize_str_representation(self) -> None: ...

class DiscreteStateDependentParameter:
    value: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]]
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
    param_type: DiscreteStateDependentParameterType
    epsilon: float
    def __init__(self, matrix_atomic_rate_params: ty.List[ty.List[DiscreteStateDependentParameter]], total_state_count: int, seed_age_for_time_slicing: ty.Optional[float] = ..., list_time_slice_age_ends: ty.Optional[ty.List[float]] = ..., epsilon: float = ...) -> None: ...
    def _check_single_and_init_param_type(self) -> None: ...
    def _check_all_states_in_all_time_slices(self) -> None: ...
    def init_matrix_state_dep_params_dict(self) -> None: ...
    def state_dep_params_at_time(self, a_time: float, params_matrix: ty.Optional[ty.List[ty.List[DiscreteStateDependentParameter]]]=None) -> ty.List[DiscreteStateDependentParameter]: ...

class MacroevolEventHandler:
    state_dep_rate_manager: DiscreteStateDependentParameterManager
    state_dep_prob_manager: ty.Optional[DiscreteStateDependentParameterManager]
    state_count: int
    n_time_slices: int
    seed_age: float
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    str_representation: str
    def __init__(self, state_dep_rate_manager: DiscreteStateDependentParameterManager, state_dep_prob_manager: ty.Optional[DiscreteStateDependentParameter]=None) -> None: ...
    def total_rate(self, a_time: float, state_representation_dict: ty.Dict[int, ty.Set[str]], value_idx: int = ..., departing_state: ty.Optional[int] = ..., debug: bool = ...): ...
    def sample_event_atomic_parameter(self, denominator: float, a_time: float, state_indices: ty.List[int], value_idx: int = ..., a_seed: ty.Optional[float] = ..., debug: bool = ...): ...
    def get_gcf(self) -> None: ...
    def get_length(self) -> None: ...

class DiscreteStateDependentProbabilityHandler:
    state_dep_prob_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.List[ty.Optional[float]]
    str_representation: str

    def __init__(self, state_dep_prob_manager: DiscreteStateDependentParameterManager) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def _state_dep_prob_at_time(self, a_time) -> ty.List[DiscreteStateDependentParameter]: ...
    def randomly_decide_taxon_sampling_at_time_at_state(self, a_time, state_idx) -> bool: ...

class SSEStash:
    meh: MacroevolEventHandler
    prob_handler: DiscreteStateDependentProbabilityHandler
    def __init__(self, macroevol_event_handler: MacroevolEventHandler, state_dep_prob_handler: ty.Optional[DiscreteStateDependentProbabilityHandler]=None) -> None: ...
    def get_meh() -> MacroevolEventHandler: ...
    def get_prob_handler() -> DiscreteStateDependentProbabilityHandler: ...