from __future__ import annotations
import typing as ty
import abc
import matplotlib.pyplot as plt # type: ignore
import numpy as np
from _typeshed import Incomplete
from abc import ABC, abstractmethod

# pj imports
# from data.tree import AnnotatedTree

R = ty.TypeVar('R')
class DummyAttribute:
    pass
def abstract_attribute(obj: ty.Callable[[ty.Any], R] = None) -> R:
    _obj = ty.cast(ty.Any, obj)
    if obj is None:
        _obj = DummyAttribute()
    _obj.__is_abstract_attribute__ = True
    return ty.cast(R, _obj)

class ProbabilisticGraphicalModel:
    node_dict: ty.Dict[NodePGM, ty.Any]
    node_name_val_dict: ty.Dict[str, NodePGM]
    n_nodes: int
    sample_size: int
    def __init__(self) -> None: ...
    def add_node(self, node_pgm: NodePGM) -> None: ...
    def get_node_pgm_by_name(self, node_name): ...
    def get_display_str_by_name(self, node_name, sample_idx: Incomplete | None = ..., repl_size: int = ...): ...
    def get_sorted_node_pgm_list(self) -> ty.List[NodePGM]: ...

class DistributionPGM(ABC, metaclass=abc.ABCMeta):

    @property
    @abstractmethod
    def DN_NAME(self):
        pass

    @abstract_attribute
    def n_draws(self):
        pass

    @abstract_attribute
    def n_repl(self):
        pass

    @abstractmethod
    def __init__(self): ...
    @abstractmethod
    def generate(self) -> ty.List[ty.Any]: ...
    @abstractmethod
    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: ...
    @abstractmethod
    def get_rev_inference_spec_info(self) -> ty.List[str]: ...

class NodePGM(ABC, metaclass=abc.ABCMeta):
    node_pgm_name: str
    value: ty.Optional[ty.List[ty.Any]]
    sample_size: int
    repl_size: int
    call_order_idx: int
    is_sampled: bool
    is_deterministic: bool
    is_clamped: bool
    parent_nd_list: Incomplete
    param_of: Incomplete
    
    @abstractmethod
    def __init__(self, node_name: str, sample_size: int, value: ty.Optional[ty.List[ty.Any]] = ..., replicate_size: int = ..., call_order_idx: ty.Optional[int] = ..., sampled: bool = ..., deterministic: bool = ..., clamped: bool = ..., parent_nodes: ty.Optional[ty.List[NodePGM]] = ...): ...
    def flatten_and_extract_values(self) -> None: ...
    def get_start2end_str(self, start: int, end: int) -> str: ...
    def __str__(self) -> str: ...
    def __hash__(self): ...
    def __eq__(self, other) -> bool: ...
    def __lt__(self, other) -> bool: ...
    def __len__(self) -> int: ...
    
    @abstractmethod
    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: int=0, repl_size: int=1, branch_attr: ty.Optional[str]="state") -> None: ...
    
    @abstractmethod
    def populate_operator_weight(self): ...

class StochasticNodePGM(NodePGM):
    is_sampled: bool
    sampling_dn: Incomplete
    operator_weight: float
    
    def __init__(self, node_pgm_name: str, sample_size: int, sampled_from: Incomplete | None = ..., value: ty.Optional[ty.List[ty.Any]] = ..., replicate_size: int = ..., call_order_idx: ty.Optional[int] = ..., deterministic: bool = ..., clamped: bool = ..., parent_nodes: ty.Optional[ty.List[ty.Any]] = ...) -> None: ...
    value: Incomplete
    def sample(self) -> None: ...
    def __lt__(self, other): ...
    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: int=0, repl_size: int=1, branch_attr: ty.Optional[str]="state") -> None: ...
    def populate_operator_weight(self) -> None: ...

class DeterministicNodePGM(NodePGM):
    is_sampled: bool
    
    def __init__(self, node_pgm_name, value: Incomplete | None = ..., call_order_idx: Incomplete | None = ..., deterministic: bool = ..., parent_nodes: Incomplete | None = ...) -> None: ...
    def __lt__(self, other) -> bool: ...
    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: int=0, repl_size: int=1, branch_attr: ty.Optional[str]="state") -> None: ...
    def populate_operator_weight(self) -> None: ...

def extract_value_from_nodepgm(val_list: ty.List[ty.Union[str, NodePGM]]) -> ty.List[str]: ...

def get_histogram_gcf(axes: plt.Axes, values_list: ty.List[float], sample_idx: ty.Optional[int] = ..., repl_size: int = ...) -> None: ...
