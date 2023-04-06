import dendropy as dp # type: ignore
import matplotlib.pyplot as plt # type: ignore
import typing as ty
from _typeshed import Incomplete
from matplotlib.figure import Figure as Figure # type: ignore

import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat

class AnnotatedTree(dp.Tree):
    with_origin: bool
    origin_node: ty.Optional[dp.Node]
    origin_age: ty.Optional[float]
    root_node: ty.Optional[dp.Node]
    root_age: ty.Optional[float]
    tree_read_as_newick_by_dendropy: bool
    tree: dp.Tree
    state_count: int
    seed_age: float
    epsilon: float
    tree_died: ty.Optional[bool]
    tree_invalid: ty.Optional[bool]
    no_event: bool
    state_count_dict: ty.Dict[int, int]
    alive_state_count_dict: ty.Dict[int, int]
    alive_sampled_state_count_dict: ty.Dict[int, int]
    dead_state_count_dict: ty.Dict[int, int]
    node_heights_dict: ty.Dict[str, float]
    node_ages_dict: ty.Dict[str, float]
    node_attr_dict: ty.Dict[str, ty.Dict[str, ty.Any]]
    slice_t_ends: ty.List[ty.Optional[float]]
    slice_age_ends: ty.Optional[ty.List[float]]
    sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]]
    at_dict: ty.Optional[ty.Dict[str, ty.List[pjat.AttributeTransition]]]
    n_extant_terminal_nodes: int
    n_extinct_terminal_nodes: int
    n_extant_sampled_terminal_nodes: int
    n_extant_sampled_terminal_nodes: int
    n_sa_obs_nodes: int
    extant_terminal_nodes_labels: ty.Tuple[str, ...]
    extinct_terminal_nodes_labels: ty.Tuple[str, ...]
    extant_sampled_terminal_nodes_labels: ty.Tuple[str, ...]
    root_edge_length: Incomplete
    def __init__(self,
                 a_tree: dp.Tree,
                 total_state_count: int,
                 start_at_origin: bool=False,
                 max_age: ty.Optional[float]=None,
                 slice_t_ends: ty.List[ty.Optional[float]]=[],
                 slice_age_ends: ty.Optional[ty.List[float]]=None,
                 sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]]=None,
                 tree_died: ty.Optional[bool]=None,
                 tree_invalid: ty.Optional[bool]=None,
                 epsilon: float=1e-12) -> None: ...
    def count_sampled_ancestors(self) -> None: ...
    def count_observable_nodes(self) -> None: ...
    def count_observable_node_states(self) -> None: ...
    def name_internal_nodes(self) -> None: ...
    def populate_node_age_height_dicts(self, unit_branch_lengths: bool = ...) -> None: ...
    def prepare_taxon_namespace_for_nexus_printing(self) -> None: ...
    def find_if_extant_or_sa_on_both_sides(self, a_node: dp.Node) -> bool: ...
    def extract_reconstructed_tree(self) -> dp.Tree: ...
    def populate_nd_attr_dict(self, attrs_of_interest_list) -> None: ...
    def plot_node(self, axes: plt.Axes, node_attr: ty.Optional[str]=None, **kwargs) -> None: ...
    def get_taxon_states_str(self, nexus: bool=False) -> str: ...

color_map: Incomplete

def get_node_name(nd): ...
def get_x_coord_from_nd_heights(ann_tr, use_age: bool = ..., unit_branch_lengths: bool = ...): ...
def get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin: bool = ...): ...
def plot_ann_tree(ann_tr: AnnotatedTree,
                  axes: plt.Axes,
                  use_age: bool=False,
                  start_at_origin: bool=False,
                  attr_of_interest: str = ...,
                  sa_along_branches: bool=True) -> None: ...
