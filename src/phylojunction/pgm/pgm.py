from __future__ import annotations
import typing as ty
import numpy as np
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.ticker as mticker  # type: ignore
import statistics as stat  # type: ignore
import pandas as pd  # type: ignore
from tabulate import tabulate  # type: ignore
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
        # keys are proper PGM nodes, values are their values
        self.node_dict = dict()

        # keys are StochasticNodePGM names
        # vals are StochasticNodePGM objects
        self.node_name_val_dict = dict()

        self.n_nodes = 0
        self.sample_size = 0  # how many simulations will be run

    def add_node(self, node_pgm: NodePGM) -> None:
        # check that nodes carry the right number of values
        # (the number of simulations)
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
                # or if we started off with a scalar node but then
                # added sampled node, we update the pgm's sample size
                if not self.sample_size or \
                        (self.sample_size == 1 and node_pgm.sample_size > 1):
                    self.sample_size = node_pgm.sample_size

                # if the number of values in a node is 1, it gets vectorized,
                # so this is allowed; but if the node is sampled and the number
                # of values is > 1 and < than that of other nodes, we have a
                # problem
                elif self.sample_size != node_pgm.sample_size and \
                        node_pgm.sample_size > 1:

                    raise ec.DAGCannotAddNodeError(
                        node_pgm.node_name,
                        ("Specified number of simulations was both "
                         "> 1 and different from that in a previous "
                         "command."))

        replacing_node = False
        if node_pgm in self.node_dict:
            # removes existing node with that node_name
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

        # if node is not deterministic and has a value and node_name
        try:
            self.node_dict[node_pgm] = node_pgm.value

        # have to double-check this can ever be a problem
        except AttributeError:
            print(("\n\n\nSomehow we could not grab a node's value."
                  "Come back to this line of code!"))
            self.node_dict[node_pgm] = None

        self.node_name_val_dict[node_pgm.node_name] = node_pgm

    def get_node_pgm_by_name(self, node_name):
        if node_name in self.node_name_val_dict:
            return self.node_name_val_dict[node_name]

    def get_display_str_by_name(
            self,
            node_name,
            sample_idx=None,
            repl_size=1):

        if node_name in self.node_name_val_dict:
            # calls __str__() of NodePGM
            return str(self.node_name_val_dict[node_name])

    def get_sorted_node_pgm_list(self) -> ty.List[NodePGM]:
        node_pgm_list: ty.List[NodePGM] = [node_pgm for node_pgm in self.node_dict]
        node_pgm_list.sort()

        return node_pgm_list

##############################################################################


class ValueGenerator(ABC):

    @abstract_attribute
    def n_samples(self):
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
    def init_check_vectorize_sample_size(
        self,
        param_list: ty.List[ty.Any] = []) \
            -> ty.Optional[ty.List[ty.List[ty.Union[int, float, str]]]]:
        """Check sample size against number of provided parameter values

        This is the function behind the vectorization functionality
        """
        pass

    @abstractmethod
    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass


class DistributionPGM(ValueGenerator):
    @property
    @abstractmethod
    def DN_NAME(self):
        pass


class ConstantFn(ValueGenerator):
    @property
    @abstractmethod
    def CT_FN_NAME(self):
        pass

##############################################################################


class NodePGM(ABC):
    # for later, I think
    # value: ty.Union[float, ty.List[ty.Union[float,T]]] = None

    # don't want to allow NodePGM to be initialized
    @abstractmethod
    def __init__(self,
                 node_name: str,
                 sample_size: int,
                 value: ty.Optional[ty.List[ty.Any]] = None,
                 replicate_size: int = 1,
                 call_order_idx: ty.Optional[int] = None,
                 sampled: bool = False,
                 deterministic: bool = False,
                 clamped: bool = False,
                 parent_nodes: ty.Optional[ty.List[NodePGM]] = None):

        self.node_name = node_name
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

        # total number of values
        # self.full_length = self.length * self.repl_size
        self.param_of = None

        if isinstance(self.value, (list, np.ndarray)) and not \
           isinstance(self.value[0], pjtr.AnnotatedTree):
            self._flatten_and_extract_values()

    # side-effect updates self.value
    def _flatten_and_extract_values(self) -> None:
        values: ty.List[ty.Any] = []

        # so mypy won't complain
        if isinstance(self.value, list):
            for v in self.value:
                if not isinstance(v, (int, float, str, np.float64)):
                    values_inside_nodes = v.value

                    # so mypy won't complain
                    if isinstance(values_inside_nodes, list):
                        for val in values_inside_nodes:
                            if not isinstance(val,
                                              (int, float, str, np.float64)):
                                raise ec.VariableAssignmentError(self.node_name)

                            values.append(val)

                else:
                    values.append(v)

        self.value = values

    # called by GUI
    def get_start2end_str(self,
                          start: int,
                          end: int,
                          repl_idx: int = 0,
                          is_tree: bool = False) -> str:        
        if isinstance(self.value, np.ndarray):
            self.value = \
                ", ".join(str(v) for v in self.value.tolist()[start:end])

        if isinstance(self.value, list):
            if len(self.value) >= 2:
                if not is_tree:
                    if isinstance(self.value[0],
                                  (int, float, str, np.float64)):
                        return ", ".join(
                            str(v) for v in self.value[start:end])

                    else:
                        return "\n".join(
                            str(v) for v in self.value[start:end])

                # not a tree
                else:
                    return str(self.value[start + repl_idx])

            # single element in value
            else:
                return str(self.value[0])

        return str(self.value)

    def __hash__(self):
        return hash(self.node_name)

    def __eq__(self, other) -> bool:
        return other.node_name == self.node_name

    # stringify _all_ values
    def __str__(self) -> str:
        # try:
        if isinstance(self.value, list):
            return self.get_start2end_str(0, len(self.value))

        # if this blows up at some point, find out
        # what kind of error is possible and
        # reimplement the except making it non-bare
        # except:
        else:
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
                    raise ec.ReplicateNumberError(node_name=self.node_name)

        # just being ultra-safe...
        # self.value has no len()
        # n_repls is None
        except TypeError:
            return 1

        # cosmetic return required by mypy
        return 1

    @abstractmethod
    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: ty.Optional[int] = None,
                  repl_idx: int = 0,
                  repl_size: int = 1,
                  branch_attr: ty.Optional[str] = "state") -> None:
        pass

    @abstractmethod
    def populate_operator_weight(self):
        pass

    @abstractmethod
    def get_node_stats_str(self,
                           start: int,
                           end: int,
                           repl_idx: int) -> str:
        stats_str = ""

        if isinstance(self.value, list):
            # could be (i) >= 1 samples with 1 replicate
            #          (ii) >= 1 samples with > 1 replicates
            # and we do (end - start) > 1 to make sure (i) doesn't go through
            if len(self.value) >= 2 and (end - start) > 1:
                if isinstance(self.value[0], (str, int, float, np.float64)):
                    # do stuff with the below
                    repl_values = [float(v) for v in self.value[start:end]]
                    repl_value_avg = stat.mean(repl_values)
                    repl_value_sd = stat.stdev(repl_values)

                    # so that it is very explicit what we are averaging over
                    avg_str, std_str = str(), str()

                    # all samples, all replicates
                    if (end - start) == len(self.value):
                        avg_str = "Sample average   = "
                        std_str = "Sample std. dev. = "

                    # one sample (all replicates)
                    else:
                        avg_str = "Replicate average   = "
                        std_str = "Replicate std. dev. = "

                    stats_str = \
                        avg_str + str(repl_value_avg) + "\n" + \
                        std_str + str(repl_value_sd)

                # values are objects (e.g., AnnotatedTree), not scalars
                else:
                    # do stuff with the below
                    repl_objs = self.value[start:end]
                    repl_all_stats_dict: ty.Dict[str, ty.List[float]] = dict()

                    # focal replicate given by repl_idx above
                    focal_repl_stats_dict: \
                        ty.Dict[str, ty.Union[int, float]] = dict()
                    # all replicates
                    repl_value_avg_dict: ty.Dict[str, float] = dict()
                    # all replicates
                    repl_value_sd_dict: ty.Dict[str, float] = dict()

                    # reading over all obj replicates and collecting
                    # their stats
                    for idx, repl_obj in enumerate(repl_objs):
                        obj_stat_dict = repl_obj.get_stats_dict()

                        for st, stat_v in obj_stat_dict.items():
                            float_stat_v: float = 0.0

                            if stat_v != "None":
                                try:
                                    float_stat_v = float(stat_v)

                                except ValueError:
                                    raise ec.NodePGMNodeStatCantFloatError(
                                        self.node_name)

                            try:
                                repl_all_stats_dict[st].append(float_stat_v)

                            except KeyError:
                                repl_all_stats_dict[st] = [float_stat_v]

                            # getting focal replicate stats
                            if idx == repl_idx:
                                focal_repl_stats_dict = obj_stat_dict

                    # now getting means and sds
                    for st, stat_v_list in repl_all_stats_dict.items():
                        repl_value_avg_dict[st] = stat.mean(stat_v_list)
                        repl_value_sd_dict[st] = stat.stdev(stat_v_list)

                    # focal replicate
                    focal_obj_repl_value_df = \
                        pd.DataFrame(focal_repl_stats_dict.items())
                    # all replicates
                    obj_repl_value_avg_df = \
                        pd.DataFrame(repl_value_avg_dict.items())
                    # all replicates
                    obj_repl_value_sd_df = \
                        pd.DataFrame(repl_value_sd_dict.items())

                    # putting it all together
                    obj_repl_stats_df = \
                        pd.concat([focal_obj_repl_value_df,
                                   obj_repl_value_avg_df.iloc[:, 1],
                                   obj_repl_value_sd_df.iloc[:, 1]],
                                  ignore_index=True, axis=1)
                    repl_stats_str = tabulate(
                        obj_repl_stats_df,
                        headers=["Summary stat.",
                                 "This replicate",
                                 "Replicates avg.",
                                 "Replicates std. dev."],
                        tablefmt="plain",
                        showindex=False
                    ).lstrip()

                    stats_str = repl_stats_str

            # looking at single replicate (from any one sample)
            elif (end - start) == 1:
                if isinstance(self.value[0], (str, int, float, np.float64)):
                    return stats_str

                # if it's an object (e.g., AnnotatedTree, or SequenceAlignment)
                else:
                    obj_stats_dict = self.value[start:end][0].get_stats_dict()
                    obj_stats_df = pd.DataFrame(obj_stats_dict.items())
                    stats_str = tabulate(
                        obj_stats_df,
                        headers=["", ""],
                        tablefmt="plain",
                        showindex=False
                    ).lstrip()

        return stats_str

##############################################################################


class StochasticNodePGM(NodePGM):

    random_value: ty.List[ty.Any]

    def __init__(self,
                 node_name: str,
                 sample_size: int,
                 sampled_from: ty.Optional[DistributionPGM] = None,
                 returned_from: ty.Optional[ConstantFn] = None,
                 value: ty.Optional[ty.List[ty.Any]] = None,
                 replicate_size: int = 1,
                 call_order_idx: ty.Optional[int] = None,
                 deterministic: bool = False,
                 clamped: bool = False,
                 parent_nodes: ty.Optional[ty.List[ty.Any]] = None):

        self.is_sampled = False
        self.sampling_dn = sampled_from  # dn object
        self.constant_fn = returned_from  # constant fn object

        if not value:
            self.random_value = self.get_value()
        
        else:
            self.random_value = value

        super().__init__(node_name,
                         sample_size=sample_size,
                         value=self.random_value,
                         replicate_size=replicate_size,
                         call_order_idx=call_order_idx,
                         deterministic=deterministic,
                         clamped=clamped,
                         parent_nodes=parent_nodes)

        # used for MCMC #move/operator setup during inference
        self.operator_weight = 0.0            

        self.populate_operator_weight()

        # pgm specs
        if self.sampling_dn:
            self.is_sampled = True  # r.v. value is sampled

    def get_value(self) -> ty.List[ty.Any]:
        if self.sampling_dn:
            return self.sampling_dn.generate()
        
        elif self.constant_fn:
            return self.constant_fn.generate()

        else:
            raise RuntimeError("exiting...")

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other):
        return super().__lt__(other)

    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: ty.Optional[int] = None,
                  repl_idx: int = 0,
                  repl_size: int = 1,
                  branch_attr: str = "state") -> None:
        """_summary_

        Args:
            axes (matplotlib.pyplot.Axes): _description_
            sample_idx (int, optional): Which sample to plot.
                Defaults to 0.
            repl_idx (int, optional): Which tree replicate to plot
                (one tree is plotted at a time). Defaults to 0.
            repl_size (int, optional): How many scalar random variables
                to plot at a time. Defaults to 1.
            branch_attr (str, optional): Which discrete attribute
                associated to a branch length to color by. Defaults to "state".
        """

        # if list
        if super().__len__() >= 1 and isinstance(self.value, list):
            # if tree, we only plot one
            if isinstance(self.value[0], pjtr.AnnotatedTree):
                # so mypy won't complain
                if isinstance(sample_idx, int) and \
                    isinstance(self.value[sample_idx * repl_size + repl_idx],
                               pjtr.AnnotatedTree):
                    if self.value[0].state_count > 1:
                        self.value[sample_idx * repl_size + repl_idx] \
                            .plot_node(axes, node_attr=branch_attr)

                    else:
                        self.value[sample_idx * repl_size + repl_idx] \
                            .plot_node(axes)

            # not a tree
            else:
                # so mypy won't complain
                hist_vals = [ty.cast(float, v) for v in self.value]

                # one sample
                if sample_idx is not None and self.sampling_dn:
                    plot_node_histogram(axes,
                                        hist_vals,
                                        sample_idx=sample_idx,
                                        repl_size=repl_size)

                # all samples
                else:
                    plot_node_histogram(axes,
                                        hist_vals,
                                        repl_size=repl_size)

        # if scalar
        elif isinstance(self.value, (int, float, str)):
            # return
            plot_node_histogram(axes, [float(self.value)])

    def populate_operator_weight(self):
        if isinstance(self.value, (list, np.ndarray)):
            if isinstance(self.value[0], (int, float, str, np.float64)):
                self.operator_weight = 1  # this is a scalar random variable

            # has objects like DiscreteStateDependentRate inside list of values
            else:
                # TODO: later see how rev moves 2D-arrays and tree nodes
                # print("value has objects inside")
                # print(self.value)
                pass

        else:
            # TODO: later see how rev moves 2D-arrays and tree nodes
            raise RuntimeError(
                ("Could not determine dimension of StochasticNodePGM when"
                 " figuring out operator weight. Exiting..."))

    def get_node_stats_str(self, start: int, end: int, repl_idx: int) -> str:
        return super().get_node_stats_str(start, end, repl_idx)

##############################################################################


class DeterministicNodePGM(NodePGM):
    def __init__(self,
                 node_name,
                 value=None,
                 call_order_idx=None,
                 deterministic=True,
                 parent_nodes=None):

        super().__init__(node_name,
                         sample_size=None,
                         value=value,
                         call_order_idx=call_order_idx,
                         deterministic=deterministic,
                         parent_nodes=parent_nodes)

        self.is_sampled = False

    def __str__(self) -> str:
        return super().__str__()

    def __lt__(self, other) -> bool:
        return super().__lt__(other)

    # deterministic nodes have all sorts of members, len() here should have no meaning
    def __len__(self) -> int:
        return 0

    def plot_node(self,
                  axes: plt.Axes,
                  sample_idx: ty.Optional[int] = None,
                  repl_idx: ty.Optional[int] = 0,
                  repl_size: ty.Optional[int] = 1,
                  branch_attr: ty.Optional[str] = "state") -> None:

        plot_blank(axes)

    def populate_operator_weight(self):
        pass

    def get_node_stats_str(self, start: int, end: int, repl_idx: int) -> str:
        return ""


############
# Plotting #
############
def plot_node_histogram(axes: plt.Axes,
                        values_list: ty.List[float],
                        sample_idx: ty.Optional[int] = None,
                        repl_size: int = 1) -> None:

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
    counts, bins, _ = \
        axes.hist(values_list_to_plot,
                  histtype='stepfilled',
                  color="steelblue",
                  edgecolor="none",
                  alpha=0.3)

    axes.spines['left'].set_visible(True)
    axes.spines['bottom'].set_visible(True)
    # axes.patch.set_alpha(0.0)
    # bins = axes.get_xticks().tolist()
    axes.xaxis.set_major_locator(mticker.FixedLocator(bins))
    # label_format = '{:,.0f}'
    # axes.set_xticklabels([label_format.format(x) for x in bins])
    plt.xticks(bins)
    plt.yticks(np.linspace(0, max(counts.tolist()), 5))
    axes.axes.xaxis.set_ticklabels([np.round(i, decimals=2) for i in bins])


def plot_blank(axes: plt.Axes) -> None:
    axes.cla()
    axes.patch.set_alpha(0.0)
    axes.xaxis.set_ticks([])
    axes.yaxis.set_ticks([])
    axes.spines['left'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)


####################
# Helper functions #
####################
def extract_value_from_nodepgm(
        val_list: ty.List[ty.Union[str, NodePGM]]) -> ty.List[str]:
    """
    Return list of values

    If all elements are strings, returns copy of 'val_list'.
    When elements are StochasticNodePGMs, replaces those objects
    by their values after casting to string (their values must be
    within a list).

    If StochasticNodePGMs objects do not have .value field or if
    they cannot be string-fied, raise exception.
    """
    extracted_val_list: ty.List[str] = []
    for v in val_list:
        if isinstance(v, str):
            extracted_val_list.append(v)

        elif isinstance(v, StochasticNodePGM) and v.value:
            try:
                extracted_val_list.extend([str(i) for i in v.value])

            except (AttributeError, TypeError) as e:
                raise ec.VariableMisspec(str(v))

    # will be empty if DeterministicNodePGM is in val_list
    return extracted_val_list


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
