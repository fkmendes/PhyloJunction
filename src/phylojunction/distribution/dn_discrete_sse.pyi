import dendropy as dp # type: ignore
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat
from phylojunction.data.tree import AnnotatedTree

class DnSSE(pgm.DistributionPGM):
    DN_NAME: str
    n_sim: int
    n_repl: int
    with_origin: bool
    root_is_born: bool
    tree_died: bool
    stop: str
    stop_val: float
    condition_on_speciation: bool
    condition_on_survival: bool
    condition_on_obs_both_sides_root: bool
    min_rec_taxa: int
    max_rec_taxa: int
    abort_at_obs: int
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
    def __init__(self,
                 event_handler: sseobj.MacroevolEventHandler = ...,
                 stop_value: ty.List[float] = ...,
                 n: int = ...,
                 n_replicates: int = ...,
                 stop: ty.Optional[str] = ...,
                 origin: bool = ...,
                 start_states_list: ty.List[int] = ...,
                 condition_on_speciation: bool = ...,
                 condition_on_survival: bool = ...,
                 condition_on_obs_both_sides_root: bool = ...,
                 min_rec_taxa: int = ...,
                 max_rec_taxa: int = ...,
                 abort_at_obs: int = ...,
                 seeds_list: ty.Optional[ty.List[int]] = ...,
                 epsilon: float = ...,
                 runtime_limit: int = ...,
                 debug: bool = ...) -> None: ...
    def _initialize_missing_prob_handler(self) -> None: ...
    def _check_sample_size(self) -> None: ...
    def get_next_event_time(self, total_rate: float, a_seed: ty.Optional[int] = ...) -> float: ...
    def execute_birth(self,
                      tr_namespace: dp.TaxonNamespace,
                      chosen_node: dp.Node,
                      state_representation_dict: ty.Dict[int, ty.Set[str]],
                      sa_lineage_dict: ty.Dict[str, ty.List[pjsa.SampledAncestor]],
                      untargetable_node_set: ty.Set[str],
                      cumulative_node_count: int,
                      macroevol_atomic_param: sseobj.DiscreteStateDependentRate,
                      event_t: float,
                      debug: bool = ...) -> ty.Tuple[dp.Node, int]: ...
    def execute_death(self,
                      chosen_node: dp.Node,
                      state_representation_dict: ty.Dict[int, ty.Set[str]],
                      sa_lineage_dict: ty.Dict[str, ty.List[pjsa.SampledAncestor]],
                      untargetable_node_set: ty.Set[str],
                      event_t: float,
                      debug: bool = ...) -> dp.Node: ...
    def execute_anatrans(self,
                         tr_namespace: dp.TaxonNamespace,
                         chosen_node: dp.Node,
                         state_representation_dict: ty.Dict[int, ty.Set[str]],
                         state_transition_dict: ty.Dict[str, ty.List[pjat.AttributeTransition]],
                         sa_lineage_dict: ty.Dict[str, ty.List[pjsa.SampledAncestor]],
                         untargetable_node_set: ty.Set[str],
                         macroevol_rate_param: sseobj.DiscreteStateDependentRate,
                         event_t: float,
                         debug: bool = ...) -> None: ...
    def execute_sample_ancestor(self,
                                tr_namespace: dp.TaxonNamespace,
                                chosen_node: dp.Node,
                                state_representation_dict: ty.Dict[int, ty.Set[str]],
                                sa_lineage_dict: ty.Dict[str, ty.List[pjsa.SampledAncestor]],
                                untargetable_node_set: ty.Set[str],
                                cumulative_sa_count: int,
                                event_t: float,
                                debug: bool = ...) -> int: ...
    def update_sa_lineage_dict(self, a_time: float, sa_lineage_dict: ty.Dict[str, ty.List[pjsa.SampledAncestor]], sa_lineage_node_labels: ty.List[str], debug: bool=False) -> None: ...
    def execute_event(self,
                      tr_namespace,
                      macroevol_rate_param: sseobj.DiscreteStateDependentRate,
                      chosen_node: dp.Node,
                      state_representation_dict: ty.Dict[int, ty.Set[str]],
                      state_transition_dict: ty.Dict[str, ty.List[pjat.AttributeTransition]],
                      untargetable_node_set: ty.Set[str],
                      cumulative_node_count: int,
                      cumulative_sa_count: int,
                      last_chosen_node,
                      event_t: float,
                      debug: bool = ...) -> ty.Tuple[dp.Node, int, int]: ...
    def simulate(self, a_start_state: int, value_idx: int = ..., a_seed: ty.Optional[int] = ...) -> AnnotatedTree: ...
    def generate(self) -> ty.List[AnnotatedTree]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
    def is_tr_ok(self, ann_tr: AnnotatedTree) -> bool: ...
