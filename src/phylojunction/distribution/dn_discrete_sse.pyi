import dendropy as dp # type: ignore
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
from phylojunction.data.tree import AnnotatedTree

class DnSSE(pgm.DistributionPGM):
    DN_NAME: str
    n_sim: int
    n_repl: int
    with_origin: bool
    stop: str
    stop_val: float
    condition_on_speciation: bool
    condition_on_survival: bool
    events: sseobj.MacroevolEventHandler
    start_states: ty.List[int]
    state_count: int
    n_time_slices: int
    slice_t_ends: ty.List[float]
    seed_age: float
    seeds: ty.List[int]
    epsilon: float
    runtime_limit: int
    debug: bool
    def __init__(self, event_handler: sseobj.MacroevolEventHandler = ..., stop_value: ty.List[float] = ..., n: int = ..., n_replicates: int = ..., stop: ty.Optional[str] = ..., origin: bool = ..., start_states_list: ty.List[int] = ..., condition_on_speciation: bool = ..., condition_on_survival: bool = ..., seeds_list: ty.Optional[ty.List[int]] = ..., epsilon: float = ..., runtime_limit: int = ..., debug: bool = ...) -> None: ...
    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    def get_next_event_time(self, total_rate: float, a_seed: ty.Optional[int] = ...) -> float: ...
    def execute_birth(self, tr_namespace: dp.TaxonNamespace, chosen_node: dp.Node, state_representation_dict: ty.Dict[int, ty.Set[str]], untargetable_node_set, cumulative_node_count: int, last_node2speciate: dp.Node, macroevol_atomic_param: sseobj.AtomicSSERateParameter, debug: bool = ...) -> ty.Tuple[dp.Node, int]: ...
    def execute_death(self, chosen_node: dp.Node, state_representation_dict: ty.Dict[int, ty.Set[str]], untargetable_node_set, debug: bool = ...) -> dp.Node: ...
    def execute_anatrans(self, chosen_node: dp.Node, state_representation_dict: ty.Dict[int, ty.Set[str]], macroevol_rate_param: sseobj.AtomicSSERateParameter, debug: bool = ...) -> None: ...
    def execute_event(self, tr_namespace, macroevol_rate_param: sseobj.AtomicSSERateParameter, chosen_node: dp.Node, state_representation_dict: ty.Dict[int, ty.Set[str]], untargetable_node_set, cumulative_node_count, last_chosen_node, debug: bool = ...) -> ty.Tuple[dp.Node, int]: ...
    def simulate(self, a_start_state: int, value_idx: int = ..., a_seed: ty.Optional[int] = ...) -> AnnotatedTree: ...
    def generate(self) -> ty.List[AnnotatedTree]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
    def is_tr_ok(self, ann_tr: AnnotatedTree) -> bool: ...
