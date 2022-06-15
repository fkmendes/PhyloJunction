from __future__ import annotations
import typing as ty
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import matplotlib.ticker as mticker # type: ignore
from abc import ABC, abstractmethod

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.data.tree as pjtr

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

# code for @abstract attribute
R = ty.TypeVar('R')
def abstract_attribute(obj: ty.Callable[[ty.Any], R] = None) -> R:
    class DummyAttribute:
        pass
    _obj = ty.cast(ty.Any, obj)
    if obj is None:
        _obj = DummyAttribute()
    _obj.__is_abstract_attribute__ = True
    return ty.cast(R, _obj)

class ProbabilisticGraphicalModel():
    node_dict: ty.Dict[NodePGM, ty.Any]
    node_name_val_dict: ty.Dict[str, NodePGM]
    n_nodes: int
    sample_size: int

    def __init__(self):
        self.node_dict = dict() # keys are proper PGM nodes, values are their values
        self.node_name_val_dict = dict() # keys are StochasticNodePGM names, vals are StochasticNodePGM objects
        self.n_nodes = 0
        self.sample_size = 0 # how many simulations will be run

    def add_node(self, node_pgm: NodePGM) -> None:
        # check that nodes carry the right number of values (the number of simulations)
        if isinstance(node_pgm, StochasticNodePGM):
            
            #############
            # Important #
            #############

            # note how only sampled nodes have any business in setting
            # the sample size of a PGM object; this means we let the users
            # fool around with nodes with assigned (fixed, clamped) values
            # through '<-'
            if node_pgm.is_sampled:
            # if the pgm's sample size is still 0,
            # or if we started off with a scalar node but then added sampled node,
            # we update the pgm's sample size
                if not self.sample_size or (self.sample_size == 1 and node_pgm.sample_size > 1):
                    self.sample_size = node_pgm.sample_size
            
            # if the number of values in a node is 1, it gets vectorized, so this is allowed;
            # but if the node is sampled and the number of values is > 1 and < than that of other nodes,
            # we have a problem;
                elif self.sample_size != node_pgm.sample_size and node_pgm.sample_size > 1:
                    print("self.sample_size = " + str(self.sample_size))
                    print("node_pgm.sample_size = " + str(node_pgm.sample_size))
                    raise RuntimeError("Number of simulations did not match. Exiting...")

        replacing_node = False
        if node_pgm in self.node_dict:
            # removes existing node with that node_pgm_name
            call_order_idx = node_pgm.call_order_idx
            del self.node_dict[node_pgm]
            replacing_node = True

        # maintaining call order
        if replacing_node:
            node_pgm.call_order_idx = call_order_idx
        
        # new rv, add to call order
        else:
            self.n_nodes += 1
            node_pgm.call_order_idx = self.n_nodes
        
        # if node is not deterministic and has a value and node_pgm_name
        try:
            self.node_dict[node_pgm] = node_pgm.value
        except:
            self.node_dict[node_pgm] = None
            
        self.node_name_val_dict[node_pgm.node_pgm_name] = node_pgm
        

    def get_node_pgm_by_name(self, node_name):
        if node_name in self.node_name_val_dict:
            return self.node_name_val_dict[node_name]


    def get_display_str_by_name(self, node_name, sample_idx=None, repl_size=1):
        if node_name in self.node_name_val_dict:
            return str(self.node_name_val_dict[node_name]) # calls __str__() of NodePGM


    def get_sorted_node_pgm_list(self) -> ty.List[NodePGM]:
        node_pgm_list: ty.List[NodePGM] = [node_pgm for node_pgm in self.node_dict]
        node_pgm_list.sort()

        return node_pgm_list

##############################################################################

class DistributionPGM(ABC):

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
    def __init__(self):
        pass

    @abstractmethod
    def generate(self) -> ty.List[ty.Any]:
        pass

    @abstractmethod
    def check_sample_size(self, param_list: ty.List[ty.Any]=[]) -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]: 
        """Check sample size against number of provided parameter values

        This is the function behind the vectorization functionality
        """
        pass

    @abstractmethod
    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass

##############################################################################

class NodePGM(ABC):
    # for later, I think
    # value: ty.Union[float, ty.List[ty.Union[float,T]]] = None

    # don't want to allow NodePGM to be initialized
    @abstractmethod
    def __init__(self, node_name: str, sample_size: int, value: ty.Optional[ty.List[ty.Any]]=None, replicate_size: int=1, call_order_idx: ty.Optional[int]=None, sampled: bool=False, deterministic: bool=False, clamped: bool=False, parent_nodes: ty.Optional[ty.List[NodePGM]]=None):
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

        if isinstance(self.value, (list, np.ndarray)) and not isinstance(self.value[0], pjtr.AnnotatedTree):
            self.flatten_and_extract_values()

    # side-effect updates self.value
    def flatten_and_extract_values(self) -> None:
        values: ty.List[ty.Any] = []
        
        # so mypy won't complain
        if isinstance(self.value, list):
            
            for v in self.value:
                if not isinstance(v, (int, float, str, np.float64)):
                    values_inside_nodes = v.value

                    # so mypy won't complain
                    if isinstance(values_inside_nodes, list):

                        for val in values_inside_nodes:
                            if not isinstance(val, (int, float, str, np.float64)):
                                raise ec.VariableAssignmentError(self.node_pgm_name)

                            values.append(val)

                else: values.append(v)

        self.value = values

    def get_start2end_str(self, start: int, end: int) -> str:
        if isinstance(self.value, np.ndarray):
            self.value = ", ".join(str(v) for v in self.value.tolist()[start:end])

        if isinstance(self.value, list):
            if len(self.value) >= 2:
                if isinstance(self.value[0], (int, float, str, np.float64)):
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
            if isinstance(self.value, list):
                return self.get_start2end_str(0, len(self.value))
        except:
            return str(self.value)

        # cosmetic return required by mypy
        return ""

    def __lt__(self, other) -> bool:
        return self.call_order_idx < other.call_order_idx

    def __len__(self) -> int:
        try:
            if isinstance(self.value, list):
                n_values = len(self.value)
            n_repls = self.repl_size

            if n_values >= 1:
                if n_values % n_repls == 0:
                    return int(n_values / n_repls)
                else:
                    raise ec.ReplicateNumberError(node_pgm_name=self.node_pgm_name)
        except:
            return 1

        # cosmetic return required by mypy
        return 1

    @abstractmethod
    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: int=0, repl_size: int=1, branch_attr: ty.Optional[str]="state") -> None:
        pass

    @abstractmethod
    def populate_operator_weight(self):
        pass

##############################################################################

class StochasticNodePGM(NodePGM):

    def __init__(self, node_pgm_name: str, sample_size: int, sampled_from: ty.Optional[DistributionPGM]=None, value: ty.Optional[ty.List[ty.Any]]=None, replicate_size: int=1, call_order_idx: ty.Optional[int]=None, deterministic: bool=False, clamped: bool=False, parent_nodes: ty.Optional[ty.List[ty.Any]]=None):
        super().__init__(node_pgm_name, sample_size=sample_size, value=value, replicate_size=replicate_size, call_order_idx=call_order_idx, deterministic=deterministic, clamped=clamped, parent_nodes=parent_nodes)
        
        self.is_sampled = False
        self.sampling_dn = sampled_from # dn object
        self.operator_weight = 0.0 # used for MCMC move/operator setup during inference
        if not value:
            self.sample()

        self.populate_operator_weight()
        
        # pgm specs
        if self.sampling_dn:
            self.is_sampled = True # r.v. value is sampled

    def sample(self):
        if self.sampling_dn:
            self.value = self.sampling_dn.generate()
        else:
            raise RuntimeError("exiting...")

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other):
        return super().__lt__(other)

    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: int=0, repl_size: int=1, branch_attr: ty.Optional[str]="state") -> None:
        """_summary_

        Args:
            axes (_type_): _description_
            sample_idx (int, optional): Which sample to plot. Defaults to 0.
            repl_idx (int, optional): Which tree replicate to plot (one tree is plotted at a time). Defaults to 0.
            repl_size (int, optional): How many scalar random variables to plot at a time. Defaults to 1.
            branch_attr (str, optional): Which discrete attribute associated to a branch length to color by. Defaults to "state".
        """
        # if list 
        if super().__len__() >= 1 and isinstance(self.value, list):
            # if tree, we only plot one
            if isinstance(self.value[0], pjtr.AnnotatedTree):
                
                # so mypy won't complain
                if isinstance(sample_idx, int) and isinstance(self.value[sample_idx*repl_size + repl_idx], pjtr.AnnotatedTree):
                    
                    if self.value[0].state_count > 1:
                        self.value[sample_idx*repl_size + repl_idx].get_gcf(axes, node_attr=branch_attr)
                    else:
                        self.value[sample_idx*repl_size + repl_idx].get_gcf(axes)
            
            # not a tree
            else:
                # so mypy won't complain
                hist_vals = [ty.cast(float, v) for v in self.value]
                
                # one sample
                if not sample_idx == None and self.sampling_dn:
                    get_histogram_gcf(axes, hist_vals, sample_idx=sample_idx, repl_size=repl_size)
                # all samples
                else:
                    get_histogram_gcf(axes, hist_vals, repl_size=repl_size)
        
        # if scalar
        elif isinstance(self.value, (int, float, str)):
            return get_histogram_gcf(axes, [float(self.value)])

    def populate_operator_weight(self):
        if isinstance(self.value, (list, np.ndarray)):
            if isinstance(self.value[0], (int, float, str, np.float64)):
                self.operator_weight = 1 # this is a scalar random variable
            # has objects like MacroevolStateDependentRateParameter inside list of values
            else:
                # TODO: later see how rev moves 2D-arrays and tree nodes
                # print("value has objects inside")
                # print(self.value)
                pass
        else:
            # TODO: later see how rev moves 2D-arrays and tree nodes
            raise RuntimeError("Could not determine dimension of StochasticNodePGM when figuring out operator weight. Exiting...")

##############################################################################

class DeterministicNodePGM(NodePGM):
    def __init__(self, node_pgm_name, value=None, call_order_idx=None, deterministic=True, parent_nodes=None):
        super().__init__(node_pgm_name, sample_size=None, value=value, call_order_idx=call_order_idx, deterministic=deterministic, parent_nodes=parent_nodes)
        self.is_sampled = False

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other) -> bool:
        return super().__lt__(other)

    def get_gcf(self, axes: plt.Axes, sample_idx: ty.Optional[int]=None, repl_idx: ty.Optional[int]=0, repl_size: ty.Optional[int]=1, branch_attr: ty.Optional[str]="state") -> None:
        # return get_blank_gcf(axes)
        pass

    def populate_operator_weight(self):
        pass


############
# Plotting #
############
def get_histogram_gcf(axes: plt.Axes, values_list: ty.List[float], sample_idx: ty.Optional[int]=None, repl_size: int=1) -> None:

    values_list_to_plot = list()
    # if not sample_idx == None:
    if isinstance(sample_idx, int):
        start = sample_idx * repl_size
        end = start + repl_size 
        values_list_to_plot = [float(v) for v in values_list[start:end]]
    else:
        values_list_to_plot = [float(v) for v in values_list]

    # figure canvas was created outside main loop in GUI
    axes.cla() 
    counts, bins, _ = axes.hist(values_list_to_plot, histtype='stepfilled', color="steelblue", edgecolor="none", alpha=0.3) 
    
    axes.spines['left'].set_visible(True)
    axes.spines['bottom'].set_visible(True)
    # axes.patch.set_alpha(0.0)
    
    # bins = axes.get_xticks().tolist()
    
    axes.xaxis.set_major_locator(mticker.FixedLocator(bins))
    
    # label_format = '{:,.0f}'
    # axes.set_xticklabels([label_format.format(x) for x in bins])
    
    plt.xticks(bins)
    plt.yticks(np.linspace(0, max(counts.tolist()), 5))
    
    axes.axes.xaxis.set_ticklabels([np.round(i,decimals=2) for i in bins])


####################
# Helper functions #
####################
def extract_value_from_nodepgm(val_list: ty.List[ty.Union[str, NodePGM]]) -> ty.List[str]:
    """
    Return copy of val_list if all elements are strings representing values.
    When elements are StochasticNodePGMs, replaces those objects by their values cast to string (their values must be within a list).
    If StochasticNodePGMs objects do not have .value field or if they cannot be string-fied, 
    raise exception.
    """
    extracted_val_list: ty.List[str] = []
    for v in val_list:
        if isinstance(v, str):
            extracted_val_list.append(v)
        
        elif isinstance(v, StochasticNodePGM) and v.value:
            try:
                extracted_val_list.extend([str(i) for i in v.value])
            except:
                raise ec.VariableMisspec(str(v))

    return extracted_val_list # will be empty if DeterministicNodePGM is in val_list


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