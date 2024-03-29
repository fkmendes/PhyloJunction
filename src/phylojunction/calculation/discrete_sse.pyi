import typing as ty
import enum
from _typeshed import Incomplete

class MacroevolEvent(enum.Enum): Incomplete
class DiscreteStateDependentParameterType(enum.Enum): Incomplete

MacroevolEventValue = ty.Literal[0, 1, 2, 3, 4, 5]
DiscreteStateDependentParameterTypeValue: ty.Literal[0, 1]

class DiscreteStateDependentParameter:
    value: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]]
    name: str
    state: int
    epoch_idx: int
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], name: str = "", state: int = 0) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def __len__(self) -> int: ...

class DiscreteStateDependentRate(DiscreteStateDependentParameter):
    departing_state: int
    arriving_state: int
    state_tuple: ty.Tuple[int]
    event: MacroevolEvent
    str_representation: str
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], event: MacroevolEvent, name: str = "", states: ty.List[int] = [], epoch_idx: int = 1) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...

class DiscreteStateDependentProbability(DiscreteStateDependentParameter):
    str_representation: str
    def __init__(self, val: ty.Union[int, float, str, ty.List[ty.Union[int, float, str]]], name: str = "", state: int = 0, epoch_idx: int = 1) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
    def __len__(self) -> int: ...

class DiscreteStateDependentParameterManager:
    matrix_state_dep_params: ty.List[ty.List[DiscreteStateDependentParameter]]
    state_dep_params_dict: ty.Dict[int, ty.List[ty.List[DiscreteStateDependentParameter]]]
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    param_type: DiscreteStateDependentParameterType
    epsilon: float
    def __init__(self, matrix_state_dep_params: ty.List[ty.List[DiscreteStateDependentParameter]], total_state_count: int, seed_age_for_time_slicing: ty.Optional[float] = None, list_time_slice_age_ends: ty.Optional[ty.List[float]] = None, epsilon: float = 1e-12) -> None: ...
    def _init_matrix_state_dep_params_dict(self) -> None: ...
    def _init_update_slice_age_and_t_ends(self, list_time_slice_age_ends) -> None: ...
    def _init_check_single_and_init_param_type(self) -> None: ...
    def _init_check_correct_number_params_per_time_slice(self) -> None: ...
    def _init_check_repeated_parameter_in_time_slice(self) -> None: ...
    def state_dep_params_at_time(self, a_time: float, params_matrix: ty.Optional[ty.List[ty.List[DiscreteStateDependentParameter]]] = None) -> ty.List[DiscreteStateDependentParameter]: ...
    def __len__(self) -> int: ...

class MacroevolEventHandler:
    sse_rate_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    seed_age: ty.Optional[float]
    str_representation: str
    def __init__(self, sse_rate_manager: DiscreteStateDependentParameterManager) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def total_rate(self, a_time: float, state_representation_dict: ty.Dict[int, ty.Set[str]], value_idx: int = 0, departing_state: ty.Optional[int] = None, debug: ty.Optional[bool] = False) -> ty.Union[float, ty.Tuple[float, ty.List[float]]]: ...
    def sample_event_sse_rate_param(self, denominator: float, a_time: float, state_indices: ty.List[int], value_idx: int = 0, a_seed: ty.Optional[float] = None, debug: ty.Optional[bool] = False) -> ty.List[DiscreteStateDependentRate]: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class DiscreteStateDependentProbabilityHandler:
    state_dep_prob_manager: DiscreteStateDependentParameterManager
    state_count: int
    n_time_slices: int
    seed_age: ty.Optional[float]
    slice_age_ends: ty.List[float]
    slice_t_ends: ty.Optional[ty.List[float]]
    str_representation: str
    def __init__(self, state_dep_prob_manager: DiscreteStateDependentParameterManager) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def _state_dep_prob_at_time(self, a_time) -> ty.List[DiscreteStateDependentParameter]: ...
    def randomly_decide_taxon_sampling_at_time_at_state(self, a_time, state_idx) -> bool: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class SSEStash:
    meh: MacroevolEventHandler
    prob_handler: DiscreteStateDependentProbabilityHandler
    def __init__(self, macroevol_event_handler: MacroevolEventHandler, state_dep_prob_handler: ty.Optional[DiscreteStateDependentProbabilityHandler]=None) -> None: ...
    def _initialize_str_representation(self) -> None: ...
    def _initialize_missing_prob_handler(self) -> None: ...
    def get_meh() -> MacroevolEventHandler: ...
    def get_prob_handler() -> DiscreteStateDependentProbabilityHandler: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...

class StateIntoPatternConverter:
    n_char: int
    n_states_per_char: int
    n_states: int
    int2set_dict: ty.Dict[str, str]
    set2int_dict: ty.Dict[str, str]
    def __init__(self, n_characters: int, n_states_per_char: int) -> None: ...
    def _initialize_dicts(self) -> None: ...
    
