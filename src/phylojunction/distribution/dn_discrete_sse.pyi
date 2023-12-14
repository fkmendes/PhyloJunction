import dendropy as dp # type: ignore
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.data.attribute_transition as pjat
from phylojunction.data.tree import AnnotatedTree
from phylojunction.data.sampled_ancestor \
    import SampledAncestor  # type: ignore
from phylojunction.data.attribute_transition \
    import AttributeTransition  # type: ignore

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
    sse_stash: sseobj.SSEStash
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
                 sse_stash: sseobj.SSEStash = ...,
                 n: int = ...,
                 n_replicates: int = ...,
                 origin: bool = ...,
                 start_states_list: ty.List[int] = ...,
                 stop: str = "",
                 stop_value: ty.List[float] = [],
                 condition_on_speciation: bool = ...,
                 condition_on_survival: bool = ...,
                 condition_on_obs_both_sides_root: bool = ...,
                 min_rec_taxa: int = ...,
                 max_rec_taxa: int = ...,
                 abort_at_obs: int = ...,
                 epsilon: float = 1e-12,
                 runtime_limit: int = 5,
                 rng_seed: ty.Optional[int] = None,
                 debug: ty.Optional[bool] = False) -> None: ...
    def _initialize_missing_prob_handler(self) -> None: ...
    def _check_sample_size(self) -> None: ...
    def _germinate_tree(self,
                        a_start_state: int,
                        with_origin: bool,
                        state_representation_dict: ty.Dict[int, ty.Set[str]],
                        untargetable_node_set: ty.Set[str]) -> ty.Tuple[dp.Tree, dp.Node, int, bool]: ...
    def _update_sa_lineage_dict(self,
                               a_time: float,
                               sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                               sa_lineage_node_labels: ty.List[str],
                               debug: bool=False) -> None: ...
    def _execute_birth(self,
                       tr_namespace: dp.TaxonNamespace,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                       untargetable_node_set: ty.Set[str],
                       cumulative_node_count: int,
                       sse_birth_rate_object: sseobj.DiscreteStateDependentRate,
                       event_t: float,
                       debug=False) -> ty.Tuple[dp.Node, int]: ...
    def _execute_death(self,
                       tr_namespace: dp.TaxonNamespace,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       untargetable_node_set: ty.Set[str],
                       event_t: float,
                       debug=False) -> dp.Node: ...
    def _execute_anatrans(self,
                          tr_namespace: dp.TaxonNamespace,
                          chosen_node: dp.Node,
                          state_representation_dict: ty.Dict[int, ty.Set[str]],
                          state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                          untargetable_node_set: ty.Set[str],
                          sse_anatrans_rate_object: sseobj.DiscreteStateDependentRate,
                          event_t: float,
                          debug: bool = False) -> None: ...
    def _execute_sample_ancestor(self,
                                 tr_namespace: dp.TaxonNamespace,
                                 chosen_node: dp.Node,
                                 state_representation_dict: ty.Dict[int, ty.Set[str]],
                                 sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                                 untargetable_node_set: ty.Set[str],
                                 cumulative_sa_count: int,
                                 event_t: float,
                                 debug: bool = False) -> int: ...
    def _execute_event(self,
                       tr_namespace,
                       sse_rate_object: sseobj.DiscreteStateDependentRate,
                       chosen_node: dp.Node,
                       state_representation_dict: ty.Dict[int, ty.Set[str]],
                       state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                       sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                       untargetable_node_set: ty.Set[str],
                       cumulative_node_count: int,
                       cumulative_sa_count: int,
                       event_t: float,
                       debug: bool = False) -> ty.Tuple[dp.Node, int, int]: ...
    def _annotate_sampled(self,
                          a_time: float,
                          living_nodes: ty.List[dp.Node],
                          sample_idx: int) -> None: ...
    def _get_next_event_time(total_rate: float) -> float: ...
    def simulate(self, a_start_state: int, value_idx: int = ..., a_seed: ty.Optional[int] = ...) -> AnnotatedTree: ...
    def generate(self) -> ty.List[AnnotatedTree]: ...
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...
    def is_tr_ok(self, ann_tr: AnnotatedTree) -> bool: ...
