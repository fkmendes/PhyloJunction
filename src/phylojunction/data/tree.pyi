import dendropy as dp # type: ignore
import matplotlib.pyplot as plt # type: ignore
import typing as ty
from _typeshed import Incomplete
from matplotlib.figure import Figure as Figure # type: ignore

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
    tree_died: bool
    state_count_dict: ty.Dict[int, int]
    node_heights_dict: ty.Dict[str, float]
    node_ages_dict: ty.Dict[str, float]
    node_attr_dict: ty.Dict[str, ty.Dict[str, ty.Any]]
    slice_t_ends: ty.List[ty.Optional[float]]
    slice_age_ends: ty.Optional[ty.List[float]]
    n_extant_obs_nodes: int
    n_extinct_obs_nodes: int
    n_sa_obs_nodes: int
    extant_obs_nodes_labels: ty.Tuple[str, ...]
    extinct_obs_nodes_labels: ty.Tuple[str, ...]
    root_edge_length: Incomplete
    def __init__(self, a_tree: dp.Tree, total_state_count: int, start_at_origin: bool = ..., max_age: ty.Optional[float] = ..., slice_t_ends: ty.List[ty.Optional[float]] = ..., slice_age_ends: ty.Optional[ty.List[float]] = ..., epsilon: float = ...) -> None: ...
    def count_sampled_ancestors(self) -> None: ...
    def count_observable_nodes(self) -> None: ...
    def count_observable_node_states(self) -> None: ...
    def name_internal_nodes(self) -> None: ...
    def populate_node_age_height_dicts(self, unit_branch_lengths: bool = ...) -> None: ...
    def populate_nd_attr_dict(self, attrs_of_interest_list) -> None: ...
    def get_gcf(self, axes, node_attr: Incomplete | None = ..., **kwargs) -> None: ...

color_map: Incomplete

def get_node_name(nd): ...
def get_x_coord_from_nd_heights(ann_tr, use_age: bool = ..., unit_branch_lengths: bool = ...): ...
def get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin: bool = ...): ...
def get_gcf_ann_tree(ann_tr: AnnotatedTree, axes: plt.Axes, use_age: bool = ..., start_at_origin: bool = ..., attr_of_interest: str = ...) -> None: ...
