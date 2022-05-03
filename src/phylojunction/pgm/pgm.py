import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/pgm/)
import typing as ty
import numpy as np
from abc import ABC, abstractmethod

# pj imports
import utility.exception_classes as ec
from data.tree import AnnotatedTree

class DistributionPGM(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def generate(self) -> ty.Union[ty.Any, ty.List[ty.Any]]:
        pass

    @abstractmethod
    def check_sample_size(self) -> None: 
        """Check sample size against number of provided parameter values

        This is the function behind the vectorization functionality
        """
        pass

    @abstractmethod
    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass

class NodePGM(ABC):
    # for later, I think
    # value: ty.Union[float, ty.List[ty.Union[float,T]]] = None

    # don't want to allow NodePGM to be initialized
    @abstractmethod
    def __init__(self, node_name: str, sample_size: int, value = None, replicate_size: int=1, call_order_idx=None, sampled: bool=False, deterministic: bool=False, clamped: bool=False, parent_nodes=None):
        self.node_pgm_name = node_name
        self.value = value
        self.sample_size = sample_size
        self.repl_size = replicate_size
        self.call_order_idx = call_order_idx
        self.is_sampled = sampled
        self.is_deterministic = deterministic
        self.is_clamped = clamped
        self.parent_nd_list = parent_nodes
        # note that when the pgm_obj adds this to its list of nodes,
        # value will be None (value is populated when we call .sample()),
        # and self.length will be = 1; we nonetheless add this call here
        # for completion (useful in debugging and testing
        #
        # as we build the PGF through a script/gui, self.populate_length()
        # is in fact called from outside through method get_length()

        # self.full_length = self.length * self.repl_size # total number of values
        self.param_of = None

        if type(self.value) in (list, np.ndarray) and not isinstance(self.value[0], pjc.AnnotatedTree):
            self.flatten_and_extract_values()

    # side-effect updates self.value
    def flatten_and_extract_values(self) -> None:
        values: ty.List[float] = []
        for v in self.value:
            if not type(v) in (int, float, str, np.float64):
                values_inside_nodes = v.value

                for val in values_inside_nodes:
                    if not type(val) in (int, float, str, np.float64):
                        raise ec.VariableAssignmentError(self.node_pgm_name)

                    values.append(val)

            else: values.append(v)

        self.value = values

    def get_start2end_str(self, start:int, end:int) -> str:
        if type(self.value) == np.ndarray:
            self.value = ", ".join(str(v) for v in self.value.tolist()[start:end])

        if type(self.value) == list:
            if len(self.value) >= 2:
                if type(self.value[0]) in (int, float, str, np.float64):
                    return ", ".join(str(v) for v in self.value[start:end])
                else:
                    return "\n".join(str(v) for v in self.value[start:end])
            else:
                return str(self.value[0])

        return str(self.value)

    def __hash__(self):
        return hash(self.node_pgm_name)

    def __eq__(self, other) -> bool:
        return other.node_pgm_name == self.node_pgm_name

    # stringify _all_ values
    def __str__(self) -> str:
        try:
            return self.get_start2end_str(0, len(self.value))
        except:
            return str(self.value)

    def __lt__(self, other) -> bool:
        return self.call_order_idx < other.call_order_idx

    def __len__(self) -> int:
        try:
            n_values = len(self.value)
            n_repls = self.repl_size

            if n_values >= 1:
                if n_values % n_repls == 0:
                    return int(n_values / n_repls)
                else:
                    raise ec.ReplicateNumberError(node_pgm_name=self.node_pgm_name)
        except:
            return 1

    @abstractmethod
    def get_gcf(self):
        pass

    @abstractmethod
    def populate_operator_weight(self):
        pass

class DeterministicNodePGM(NodePGM):
    def __init__(self, node_pgm_name, value=None, call_order_idx=None, deterministic=True, parent_nodes=None):
        super().__init__(node_pgm_name, sample_size=None, value=value, call_order_idx=call_order_idx, deterministic=deterministic, parent_nodes=parent_nodes)
        self.is_sampled = False

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other) -> bool:
        return super().__lt__(other)

    def get_gcf(self, axes, idx_value=0, **kwargs):
        # return get_blank_gcf(axes)
        pass

    def populate_operator_weight(self):
        pass

##############################################################################

if __name__ == "__main__":
    # can be called from pgm/
    # $ python3 pgm.py
    # 
    # can also be called from phylojunction/
    # $ python3 pgm/pgm.py
    # or
    # $ python3 -m pgm.pgm
    #
    # can also be called from VS Code, if open folder is phylojuction/

    pass