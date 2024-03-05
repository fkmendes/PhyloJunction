import typing as ty
import dendropy as dp  # type: ignore
import matplotlib  # type: ignore
import copy  # type: ignore
import numpy as np
import pandas as pd
import collections
# from dendropy import Node, Tree, Taxon

# plotting tree
import matplotlib.pyplot as plt  # type: ignore
import matplotlib.collections as mpcollections  # type: ignore
from matplotlib.figure import Figure  # type: ignore
from matplotlib import colors

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
import phylojunction.data.sampled_ancestor as pjsa
import phylojunction.data.attribute_transition as pjat

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class AnnotatedTree(dp.Tree):
    """Tree annotated with discrete states. 

    Parameters:
        tree (dendropy.Tree): Main class member, holding the full tree.
        tree_reconstructed (dendropy.Tree): Reconstructed tree
            produced by pruning the full tree from non-observed
            phylogenetic paths.
        origin_node (dendropy.Node): Origin node. It is 'None' if no
            origin node.
        root_node (dendropy.Node): Root node. It is 'None' if no root
            node.
        brosc_node (dendropy.Node): Node that occurs before the root
            node when tree starts from origin node. The brosc_node can
            be a terminal or internal node, and is replaced by the root
            node if a speciation event happens.
        alternative_root_label (str): An AnnotatedTree may be read in
            with a constant function rather than be simulated by
            discrete_sse(). In this case, the root node is named 'root'
            so that AnnotatedTree is initialized properly, but we want
            to keep track of an alternative label for the root node,
            which may appear in a stochastic map table, or in a table
            with tip states. Defaults to ''.
        with_origin (bool): Flag indicating that tree starts at origin
            node. Upon instantiation of class, defaults to 'False'. 
        tree_read_as_newick (bool): Flag indicating that tree was not
            simulated in PJ, but read as a Newick string instead.
            Upon instantiation of class, defaults to 'False'.
        condition_on_obs_both_sides_root (bool): Flag indicating if
            root node, when simulation starts from it (!), has at least
            one living on either side. If not specified by user upon
            instantiation of class, defaults to 'False'.
        tree_died (bool): Flag indicating if full tree went extinct
            before reaching the specified age stopping condition (i.e.,
            the reconstructed tree should be empty).
        tree_invalid (bool): Flag indicating (if tree was simulated)
            that rejection sampling deemed the tree invalid. If not
            specified upon class instantiation, tree is assumed to have
            been read as a string, and flag is assigned 'True'.
        seed_age (float): Age of seed node (either origin or root).
        max_age (float): Maximum age of the tree when it was simulated
            with an age stopping condition. If not provided by user
            upon class instantiation, it is 'None'.
        origin_age (float): Age of origin node, if there is one,
            otherwise 'None'.
        origin_edge_length (float): Length of the branch connecting the
            origin node to its single child (brosc or root) or to the
            present. If there is a root, this length does not care
            about intervening direct (sampled) ancestor between the
            origin and the root. Will be 0.0 if no origin node.
        root_age (float): Age of root node, if there is one, otherwise
            0.0.
        node_heights_dict (dict): Dictionary that holds the heights of
            all nodes in the tree. Keys are node labels, values are
            floats.
        node_ages_dict (dict): Dictionary that holds the ages of all
            nodes in the tree. Keys are node labels, values are
            floats.
        slice_t_ends (float): List of floats with the end times for
            specified time slices (epochs). If not provided by user
            Upon instantiation of class, will be 'None'.
        slice_age_ends (dict): List of floats with the end ages for
            specified time slices (epochs). If not provided by user
            Upon instantiation of class, will be 'None'.
        state_count (int): How many states there are.
        state_count_dict (dict): Dictionary tabulating how many
            terminal nodes in the full tree are in each state. Keys are
            integers representing states and values are their counts.
        extant_terminal_state_count_dict (dict): Dictionary tabulating
            how many living (terminal) nodes are in each state. Keys
            are integers representing states and values are their
            counts.
        extant_sampled_terminal_state_count_dict (dict): Dictionary
            tabulating how many living (terminal) and sampled nodes are
            in each state. Keys are integers representing states and
            values are their counts.
        extinct_sampled_state_count_dict (dict): Dictionary tabulating
            how many extinct (terminal) nodes are in each state. Keys
            are integers representing states and values are their
            counts.
        sa_node_state_count_dict (dict): Dictionary tabulating
            how many direct (sampled) ancestor nodes are in each state.
            Keys are integers representing states and values are their
            counts.
        node_attr_dict (dict): Nested dictionaries. Keys of the outer
            dictionary are node labels, values are inner dictionaries.
            Inner dictionary's key are attribute names (str) and
            values are the attribute's values (Any).
        n_extant_terminal_nodes (int): Count of living terminal nodes.
        n_extinct_terminal_nodes (int): Count of dead terminal nodes.
        n_extant_sampled_terminal_nodes (int): Count of living and
            sampled terminal nodes.
        n_sa_nodes (int): Count of direct (sampled) ancestors. These
            are nodes that are both sampled and observed by
            definition.
        epsilon (float): Float threshold to determine if a
            tiny decimal number is to be considered 0.0 or not. In
            other words, if the difference between a tiny value 'x'
            and 0.0 is smaller than epsilon, then 'x' is set to 0.0.
            If not provided by user upon initialization of class,
            defaults to 1e-12.
    """
    
    tree: dp.Tree
    tree_reconstructed: dp.Tree

    # tree members #
    origin_node: ty.Optional[dp.Node]
    root_node: ty.Optional[dp.Node]
    brosc_node: ty.Optional[dp.Node]
    alternative_root_label: str

    # flags
    with_origin: bool
    tree_read_as_newick: bool
    condition_on_obs_both_sides_root: bool  # for rec tree
    tree_died: ty.Optional[bool]
    tree_invalid: ty.Optional[bool]

    # age related #
    seed_age: float
    max_age: ty.Optional[float]
    origin_age: ty.Optional[float]
    origin_edge_length: float
    root_age: float  # can be None
    node_heights_dict: ty.Dict[str, float]
    node_ages_dict: ty.Dict[str, float]
    slice_t_ends: ty.Optional[ty.List[float]]
    slice_age_ends: ty.Optional[ty.List[float]]

    # state related #
    state_count: int
    state_count_dict: ty.Dict[int, int]
    extant_terminal_state_count_dict: ty.Dict[int, int]
    extant_terminal_sampled_state_count_dict: ty.Dict[int, int]
    extinct_terminal_state_count_dict: ty.Dict[int, int]
    # { node_name: { attr: value }}

    # tree size related #
    node_attr_dict: ty.Dict[str, ty.Dict[str, ty.Any]]
    n_extant_terminal_nodes: int
    n_extinct_terminal_nodes: int
    n_extant_sampled_terminal_nodes: int
    n_extinct_sampled_terminal_nodes: int
    n_sa_nodes: int

    # node labels #
    extant_terminal_nodes_labels: ty.Tuple[str, ...]
    extinct_terminal_nodes_labels: ty.Tuple[str, ...]
    extant_sampled_terminal_nodes_labels: ty.Tuple[str, ...]
    sa_obs_nodes_labels: ty.Tuple[str, ...]

    # for plotting #
    sa_lineage_dict: \
        ty.Optional[
            ty.Dict[str, ty.List[pjsa.SampledAncestor]]]  # can be None
    # attribute transitions (along branches, including the ones simultaneous
    # with speciation)
    at_dict: \
        ty.Optional[ty.Dict[
            str, ty.List[pjat.AttributeTransition]]]  # can be None
    clado_at_dict: \
        ty.Optional[ty.Dict[str, pjat.AttributeTransition]]  # can be None

    # to deal with effectively zero floats
    epsilon: float

    def __init__(
            self,
            a_tree: dp.Tree,
            total_state_count: int,
            start_at_origin: bool = False,
            alternative_root_label: str = "",
            condition_on_obs_both_sides_root: bool = False,
            max_age: ty.Optional[float] = None,
            slice_t_ends: ty.Optional[ty.List[float]] = [],
            slice_age_ends: ty.Optional[ty.List[float]] = None,
            sa_lineage_dict:
            ty.Optional[
                ty.Dict[str, ty.List[pjsa.SampledAncestor]]] = None,
            at_dict:
            ty.Optional[
                ty.Dict[str, ty.List[pjat.AttributeTransition]]] = None,
            clado_at_dict:
            ty.Optional[ty.Dict[str, pjat.AttributeTransition]] = None,
            tree_died: ty.Optional[bool] = None,
            tree_invalid: ty.Optional[bool] = None,
            read_as_newick_string: bool = False,
            epsilon: float = 1e-12):

        # using during initialization
        self.epsilon = epsilon

        # trees
        self.tree = a_tree
        self.tree_reconstructed = None

        # check if simulated, or specified by user
        self.tree_read_as_newick = read_as_newick_string
        # but if simulated, nodes must have .alive member, if not,
        # we update flag member
        if hasattr(self.tree.seed_node, "alive"):
            self.tree.seed_node.alive

        else:
            self.tree_read_as_newick = True
        
        # initializing tree members
        self.origin_node = None
        self.root_node = None
        self.brosc_node = None
        self.alternative_root_label = alternative_root_label

        # flags
        self.with_origin = start_at_origin
        self.tree_died = tree_died
        self.condition_on_obs_both_sides_root = \
            condition_on_obs_both_sides_root
        # rejection sampling when stop condition is age
        self.tree_invalid = tree_invalid
        if not isinstance(self.tree_invalid, bool):
            # if flag not passed, tree was not simulated, assume valid
            self.tree_invalid = False

        # checking input health before doing more
        self._check_input_health()

        # age and time related
        self.origin_edge_length = 0.0
        self.seed_age = self.tree.max_distance_from_root()
        self.max_age = max_age
        self.node_heights_dict = dict()
        self.node_ages_dict = dict()
        self.slice_t_ends = slice_t_ends
        self.slice_age_ends = slice_age_ends

        # state related
        self.state_count = total_state_count
        self.state_count_dict = \
            dict((int(s), 0) for s in range(0, self.state_count, 1))
        self.extant_terminal_state_count_dict = \
            dict((int(s), 0) for s in range(0, self.state_count, 1))
        self.extant_terminal_sampled_state_count_dict = \
            dict((int(s), 0) for s in range(0, self.state_count, 1))
        self.extinct_terminal_state_count_dict = \
            dict((int(s), 0) for s in range(0, self.state_count, 1))
        self.sa_state_count_dict = \
            dict((int(s), 0) for s in range(0, self.state_count, 1))
        # TODO: later deal with this

        # TODO: add argument for node_attr_dict (when attrs are passed in
        # by user when reading newick tree) -- in this case, no need to
        # initialize the dictionary here, and also no need to populate it
        # in other methods
        # 
        # all attributes
        # populated on demand by populate_nd_attr_dict,
        # cannot be pickled for some reason...
        # self.node_attr_dict = pjh.autovivify(2)
        self.node_attr_dict = \
            collections.defaultdict(
                pjh.create_str_defaultdict)  # this one can be pickled

        # initializes (side-effect):
        # self.origin_age
        # self.origin_edge_length
        # self.root_age
        # self.tree_died
        # self.brosc_node
        self._init_and_update_origin_root_members()

        # for plotting and counting nodes
        self.sa_lineage_dict = sa_lineage_dict
        self.at_dict = at_dict
        self.clado_at_dict = clado_at_dict

        # node counting
        self.n_extant_terminal_nodes = 0
        self.n_extinct_terminal_nodes = 0
        self.n_extant_sampled_terminal_nodes = 0
        self.n_sa_nodes = 0

        # initializes (side-effect):
        # (i)   self.n_extant_terminal_nodes
        # (ii)  self.n_extant_sampled_terminal_nodes
        # (iii) self.n_extinct_sampled_terminal_nodes
        # (iv)  self.extant_terminal_nodes_labels
        # (v)   self.extant_sampled_terminal_nodes_labels
        # (iv)  self.extinct_sampled_terminal_nodes_labels
        self._count_terminal_nodes()

        # initializes (side-effect):
        # (i)  self.sa_obs_nodes_labels
        # (ii) self.n_sa_nodes
        self._count_sampled_ancestors()

        # initializes (side-effect):
        # (i)   self.state_count_dict
        # (ii)  self.extant_terminal_state_count_dict
        # (iii) self.extant_terminal_sampled_state_count_dict
        # (iv)  self.sa_state_count_dict
        self._count_node_states()

        # initializes (side-effect):
        # (i)  self.node_heights_dict
        # (ii) self.node_ages_dict
        self._populate_node_age_height_dicts()

        # prepare for dendropy's Nexus printing
        self._prepare_taxon_namespace_for_nexus_printing()

    def _check_input_health(self) -> None:
        """Check the validity of some of initialization arguments.
        
        Raises:
            AnnotatedTreeMisspec: Is raised if tree is said to have an
                origin node, but instead starts from the root, or if
                it is said to start from the root, but its oldest node
                has a single child.
        """

        if self.with_origin:
            if len(self.tree.seed_node.child_nodes()) > 1:
                raise ec.AnnotatedTreeMisspecError(
                        ("Argument 'with_origin' evaluated to 'True', but "
                        "seed node had two children (it was a root node)."))

        elif len(self.tree.seed_node.child_nodes()) == 1:
                raise ec.AnnotatedTreeMisspecError(
                    ("Argument 'with_origin' evaluated to 'False', but "
                    "seed node had a single child (it was an origin node)."))
        
        # making sure we have all attributes necessary for
        # successfully initializing AnnotatedTree
        for nd in self.tree.preorder_node_iter():
            # first we check that we have all attributes in place
            if not hasattr(nd, "state"):
                raise ec.AnnotatedTreeNodeMissingAttrError(
                    nd.label,
                    "state",
                    "Issue happened when initializing AnnotatedTree"
                )
            
            if not hasattr(nd, "alive"):
                raise ec.AnnotatedTreeNodeMissingAttrError(
                    nd.label,
                    "alive",
                    "Issue happened when initializing AnnotatedTree"
                )
            
            if not hasattr(nd, "sampled"):
                raise ec.AnnotatedTreeNodeMissingAttrError(
                    nd.label,
                    "sampled",
                    "Issue happened when initializing AnnotatedTree"
                )
            
            if not hasattr(nd, "is_sa"):
                raise ec.AnnotatedTreeNodeMissingAttrError(
                    nd.label,
                    "is_sa",
                    "Issue happened when initializing AnnotatedTree"
                )

    def _has_tree_died(self, origin_or_root_node: dp.Node) -> bool:
        """Check if tree has died (all terminals died).

        This method will return 'False' if the tree was read as a
        Newick string. The assumption is that the user would not
        enter a tree if it knew all lineages died.
        
        If the tree was simulated by PJ, then every node must have
        attribute '.alive'. This method then returns 'True' if it
        verifies that for all terminals '.alive == False'.

        Args:
            origin_or_root_node (dendropy.Node): Either the origin or
                the root node

        Returns:
            (bool): Whether or not the tree died
        """

        # if tree simulated, we scan all tips
        if not self.tree_read_as_newick:
            for nd in origin_or_root_node.leaf_iter():
                if not nd.is_sa and nd.alive:
                    return False
        
        # if tree not simulated (read as newick string),
        # we assume it is alive
        else:
            return False

        # if tree was simulated (not read as newick string),
        # and none of the non-SA tips is alive, then it must
        # be dead
        return True

    def _recursively_find_node_time(
            self,
            a_node: dp.Node,
            reference_time: float = 0.0) -> float:
        """
        Find a node's time recursively and return it.
        
        This method travels a node's path all the way to the start of
        the process (root or origin node) and returns path length.

        Args:
            a_node (dendropy.Node): The node whose age we want to get.
            running_time_sum (float): Point of reference from which to
                compute time (0.0 if the start of the process). This
                can also be seen as a running time sum. Defaults to
                0.0.

        Returns:
            (float): The node's time relative to a reference time.
        """

        # stop recursion
        if a_node == self.tree.seed_node:
            return reference_time

        # recur (reference_time is a running time sum)
        return self._recursively_find_node_time(
            a_node.parent_node,
            reference_time = reference_time + a_node.edge_length)

    def _init_and_update_origin_root_members(self) -> None:
        """Initialize class members related to the root and origin.

        This method initializes class members related to the origin and
        root node. It also flags if the tree died. The initialization
        side-effects are on:
            (i)   origin_age
            (ii)  origin_edge_length (the distance between the origin
                  node and its single child, i.e., either the brosc
                  node or the root, or between the root and the present
                  moment). If there is a root, the length of this edge
                  does not stop at direct (sampled) ancestors.
            (iii) root_age
            (iv)  tree_died (all terminal nodes went extinct)
            (v)   brosc_node
        """

        if self.with_origin:
            self.origin_node = self.tree.seed_node
            self.root_age = 0.0
            origin_children: ty.List[dp.Node] = list()

            for nd in self.tree.preorder_node_iter():
                if nd.label != "origin":
                    origin_children.append(nd)
                
                    if nd.label == "brosc":
                        self.brosc_node = nd
            
            origin_children_labels: ty.Tuple[str] = \
                tuple([nd.label for nd in origin_children])
            
            # origin_children = \
            #     [nd for nd in self.tree.nodes() if nd.label != "origin"]
            # origin_children_labels = \
            #     tuple([nd.label for nd in origin_children])

            # # if there is a brosc node, we get it
            # try:
            #     self.brosc_node = \
            #         [nd for nd in origin_children if nd.label == "brosc"][0]

            # except Exception as e:
            #     # print("Exception 1 inside tree.py: ",
            #     #       type(e).__name__, " - ", e)
            #     pass  # self.brosc_node will remain None

            # Case (a): no events took place, tree may or not have gone extinct
            #
            # [origin] --- [brosc] ... [max age] (died)
            # [origin] ------ [brosc at max age] (survived)
            if self.origin_node.num_child_nodes() == 1 and \
                    origin_children_labels[0] == "brosc":
                # tree is just the origin edge, which is also the origin age
                self.origin_age = \
                    self.origin_edge_length = origin_children[0].edge_length

            # Case (b): at least one event took place (root may or not have
            # been born) and tree may or not have died
            elif len(origin_children) > 1:
                # Case (b.1) There is a root
                #
                # b.1.1 (died):
                # [origin] - [root + children] ............. [stop condition]
                # b.1.2 (survived):
                # [origin] -------------- [root + children at stop condition]
                if "root" in origin_children_labels or \
                        self.alternative_root_label in origin_children_labels:
                    self.origin_age = self.seed_age

                    # because the origin may have immediate children who are
                    # not the root, we need to specifically ask for the root
                    self.root_node = \
                        [nd for nd in origin_children \
                         if nd.label in ("root", self.alternative_root_label)][0]

                    # the origin edge is by definition the edge between the
                    # origin and the root; so by finding the time of the root
                    # we get the origin edge length
                    self.origin_edge_length = \
                        self._recursively_find_node_time(self.root_node)

                    # if there is a root but the tree dies, the root age will
                    # be 0.0 (the origin age and the origin edge length are
                    # the same)
                    self.root_age = \
                        self.origin_age - self.origin_edge_length

                # Case (b.2) There is no root, so event(s) must have been
                # ancestor sampling
                #
                # b.2.1 (died):
                # [origin] - [anc. sampling(s)] - [brosc] ... [max age]
                # b.2.2 (survived):
                # [origin] - [anc. sampling(s)] ---- [brosc at max age]
                else:
                    for nd in self.tree.seed_node.leaf_iter():
                        if not nd.is_sa and abs(nd.edge_length) > self.epsilon:
                            self.origin_age = \
                                self._recursively_find_node_time(nd)

                    self.origin_edge_length = \
                        self._recursively_find_node_time(self.brosc_node)

            # updating tree_died member
            # (it was not specified upon initialization)
            if not isinstance(self.tree_died, bool):
                self.tree_died = True

                # with max_age
                # (origin_age should have been set above but we double
                # check here)
                if isinstance(self.max_age, float) and \
                    isinstance(self.origin_age, float) and \
                        (self.max_age - self.origin_age) <= self.epsilon:
                    self.tree_died = False

                # no max_age or no origin_age
                else:
                    # but there is a brosc_node
                    # if self.brosc_node is not None:

                
                    # there is a brosc node
                    #
                    # who is alive
                    if isinstance(self.brosc_node, dp.Node):
                        if self.brosc_node.alive:
                            self.tree_died = False

                        # or dead
                        else:
                            self.tree_died = True

                    # no brosc_node
                    elif self.brosc_node is None:
                        # _has_tree_died will return 
                        self.tree_died = self._has_tree_died(self.origin_node)

                        # will get here if node doesn't have .alive
                        # except Exception as e:
                            # print("Exception 2 inside tree.py: ", type(e).__name__, " - ", e)
                            # pass

                    # no way to tell, we assume tree died
                    # else:
                    #     self.tree_died = True

            # a bit of cleaning for printing the tree
            if self.tree_read_as_newick and len(self.tree.leaf_nodes()) == 1 \
                    and isinstance(self.root_node, dp.Node):
                self.root_node.label = \
                    str(self.root_node.taxon).replace("\'", "")

        # starting at root
        else:
            # if tree starts at root node, and that node has an
            # edge_length, it must mean that there is an origin and
            # that the tree was probably read as a newick string, and
            # somehow 'start_at_origin' was not set to 'True' upon
            # initialization
            #
            # below we must then create an origin node and update
            # everything accordingly
            if len(self.tree.seed_node.child_nodes()) == 2 and \
                    self.tree.seed_node.edge_length:
                self.root_node = self.tree.seed_node
                self.origin_edge_length = self.root_node.edge_length
                self.origin_node = dp.Node(taxon=dp.Taxon(label="origin"),
                                           label="origin",
                                           edge_length=0.0)
                self.origin_node.alive = False
                self.origin_node.add_child(self.root_node)
                self.tree.seed_node = self.origin_node
                self.seed_age = self.tree.max_distance_from_root()
                self.origin_age = self.seed_age
                self.root_age = self.origin_age - self.origin_edge_length
                
                # tree is alive since it was read as a Newick string,
                # unless tree_died was passed upon initialization)
                if not isinstance(self.tree_died, bool):
                    self.tree_died = False

            # tree was built by simulator, no origin node
            else:
                self.origin_age = None
                self.origin_edge_length = 0.0
                self.root_node = self.tree.seed_node
                self.root_age = self.tree.max_distance_from_root()

                # updating tree_died member
                # (it was not specified upon initialization)
                if not isinstance(self.tree_died, bool):
                    self.tree_died = True

                    # with max_age
                    if isinstance(self.max_age, float) and \
                        isinstance(self.root_age, float) and self.max_age \
                        and (self.max_age - self.root_age) \
                            <= self.epsilon:
                        self.tree_died = False

                    # no max_age and no root age
                    else:
                        # try:
                        self.tree_died = self._has_tree_died(self.root_node)

                        # will get here if node doesn't have .alive
                        # except Exception as e:
                            # print("Exception 3 inside tree.py: ", type(e).__name__, " - ", e)
                            # pass

                # fixes root label if tree has only origin + root
                if self.tree_read_as_newick and \
                        len(self.tree.leaf_nodes()) == 1:
                    self.root_node.label = self.root_node.taxon

    def _count_terminal_nodes(self) -> None:
        """Count extant and extinct terminal nodes.

        This method counts extant and extinct terminal nodes and
        updates the respective class members. By definition, direct
        (sampled) ancestors are not terminal nodes, so they are
        ignored.

        This method depends on tree_died having been updated
        correctly.

        The initialization side-effects are on:
            (i)   self.n_extant_terminal_nodes
            (ii)  self.n_extant_sampled_terminal_nodes
            (iii) self.n_extinct_terminal_nodes
            (iv)  self.extant_terminal_nodes_labels
            (v)   self.extant_sampled_terminal_nodes_labels
            (vi)  self.extinct_terminal_nodes_labels
        """

        # nd.distance_from_root() gives distance to seed!
        extant_terminal_nd_labels_list: ty.List[str] = []
        extinct_terminal_nd_labels_list: ty.List[str] = []
        extant_sampled_terminal_nd_labels_list: ty.List[str] = []

        # tree died!
        #
        # no extant taxa, so counting extinct terminal nodes only
        if self.tree_died:
            # no root, single lineage from origin died
            if self.brosc_node:
                self.n_extinct_terminal_nodes = 1  # brosc node
                extinct_terminal_nd_labels_list \
                    .append(self.brosc_node.label)

            # with root, every lineage died
            else:
                for nd in self.tree.leaf_node_iter():
                    # making sure it is not a direct (sampled) ancestor
                    if not nd.is_sa and not nd.alive:
                        self.n_extinct_terminal_nodes += 1
                        extinct_terminal_nd_labels_list.append(nd.label)

        # tree did not die!
        #
        # there can be both extant and extinct terminal nodes,
        # so counting both
        else:
            for nd in self.tree.leaf_node_iter():
                # if not sampled ancestor
                if not (nd.edge_length < self.epsilon):
                    # counting extant terminal nodes
                    # 
                    # an extant terminal node is a leaf whose path
                    # to the origin/root has maximal length (equal
                    # to the age of the origin/root)
                    #
                    # nd.distance_from_root() returns distance to seed
                    if abs(self.seed_age - nd.distance_from_root()) \
                            <= self.epsilon:
                        if self.tree_read_as_newick:
                            self.n_extant_terminal_nodes += 1
                            extant_terminal_nd_labels_list. \
                                append(nd.taxon.label)

                        # if tree was created via simulation
                        # (and .alive has been set as an attribute)
                        elif nd.alive:
                            self.n_extant_terminal_nodes += 1
                            extant_terminal_nd_labels_list.append(nd.label)

                            if hasattr(nd, "sampled") and nd.sampled:
                                self.n_extant_sampled_terminal_nodes += 1
                                extant_sampled_terminal_nd_labels_list.\
                                    append(nd.label)

                    # counting extinct terminal nodes
                    #
                    # if terminal node's path is not maximal, it must
                    # have gone extinct
                    else:
                        if self.tree_read_as_newick:
                            self.n_extinct_terminal_nodes += 1
                            extinct_terminal_nd_labels_list \
                                .append(nd.taxon.label)

                        # if tree was created via simulation
                        # (and .alive a member of Node)
                        elif not nd.alive:
                            self.n_extinct_terminal_nodes += 1
                            extinct_terminal_nd_labels_list \
                                .append(nd.label)

                        else:
                            raise ec.AnnotatedTreeIncorrectAnnotationError(
                                ("Taxon had non-maximal age, but had "
                                 "\'.alive == True\'. This is not allowed. "
                                 "Exiting..."))

        if self.n_extant_terminal_nodes == 0 and not self.tree_died:
            raise ec.AnnotatedTreeIncorrectAnnotationError(
                "Zero extant terminal nodes were found, but tree was found to"
                " not be dead."
            )

        self.extant_terminal_nodes_labels = \
            tuple(extant_terminal_nd_labels_list)
        self.extant_sampled_terminal_nodes_labels = \
            tuple(extant_sampled_terminal_nd_labels_list)
        self.extinct_terminal_nodes_labels = \
            tuple(extinct_terminal_nd_labels_list)

    def _count_sampled_ancestors(self) -> None:
        """Count direct (sampled) ancestor nodes.

        This method counts direct (sampled) ancestor nodes, and updates
        the appropriate class members. Note that direct (sampled)
        ancestors are always sampled by definition, and thus always
        observed.
        
        This method depends on 'sa_lineage_dict' class member having
        been provided upon tree initialization.

        The initialization side-effects are on:
            (i)   self.n_sa_nodes
            (ii)  self.sa_obs_nodes_labels_list
        """

        sa_obs_nodes_labels_list: ty.List[str] = []

        if isinstance(self.sa_lineage_dict, dict):
            for _, sa_list in self.sa_lineage_dict.items():
                self.n_sa_nodes += len(sa_list)
                sa_obs_nodes_labels_list.extend([sa.label for sa in sa_list])

        self.sa_obs_nodes_labels = tuple(sa_obs_nodes_labels_list)

    def _count_node_states(self) -> None:
        """Count how many nodes are in each state.

        This method only visits terminal nodes. Importantly, direct
        (sampled) ancestors are represented as terminal nodes, but
        are not really terminal (nor extinct depending on how one sees)
        it.

        The initialization side-effects are on:
            (i)   self.state_count_dict
            (ii)  self.extant_terminal_state_count_dict
            (iii) self.extant_terminal_sampled_state_count_dict
            (iv)  self.sa_state_count_dict
        """

        # iterate over all terminal nodes, including nodes
        # that are not really terminal (direct ancestors),
        # but that are nonetheless implemented as such
        for nd in self.tree.leaf_node_iter():            
            # now we count
            self.state_count_dict[nd.state] += 1

            if nd.alive:
                self.extant_terminal_state_count_dict[nd.state] += 1

                if nd.sampled:
                    self.extant_terminal_sampled_state_count_dict[nd.state] += 1

            else:
                if nd.is_sa:
                    self.sa_state_count_dict[nd.state] += 1
                
                else:
                    self.extinct_terminal_state_count_dict[nd.state] += 1    

    def _populate_node_age_height_dicts(
            self,
            unit_branch_lengths: bool = False) -> None:
        """Recursively populate class members with nodes' ages.
        
        This method has the side-effect of populating two class members,
        node_ages_dict and node_heights_dict, whose keys are node names
        (str) and values (float) are their ages and heights,
        respectively.

        The age of the origin is maximal, and the height of the extant
        terminal nodes is maximal (and the same as the origin's age).
        """

        # side-effect:
        # populates self.node_heights_dict { node label (str): node age (float, ... }
        # or
        # self.node_heights_dict AND self.node_ages_dict { node label (str): node age (float, ... }
        def recur_node_ages_height(nd,
                                   remaining_height,
                                   tree_lvl) -> None:

            for ch_nd in nd.child_nodes():
                recur_node_ages_height(ch_nd,
                                       remaining_height - ch_nd.edge_length,
                                       tree_lvl + 1)

            nd_name: str = get_node_name(nd)

            if unit_branch_lengths:
                self.node_heights_dict[nd_name] = tree_lvl
            
            else:
                self.node_ages_dict[nd_name] = remaining_height
                self.node_heights_dict[nd_name] = \
                    self.seed_age - remaining_height

        # tree_lvl starts at 0
        recur_node_ages_height(self.tree.seed_node, self.seed_age, 0)

    def _prepare_taxon_namespace_for_nexus_printing(self) -> None:
        """Clean up taxon namespace for printing.
        
        Prepares taxon_namespace member of self.tree so Nexus produced
        by DendroPy is reasonable. Removing internal nodes from
        taxon_namespace does not affect a tree's Newick representation.
        """

        if not self.tree_read_as_newick:
            for nd in self.tree.nodes():
                if not nd.is_leaf():
                    self.tree.taxon_namespace.remove_taxon(nd.taxon)

    def is_extant_or_sa_on_both_sides_complete_tr_root(self, a_node: dp.Node) -> bool:
        """Verify one or more sampled nodes exist on both root sides.

        This method is called by extract_reconstructed_tree(), and
        verifies that there is at least one sampled node (direct
        sampled ancestor, or sampled extant) on both sides of the
        root of the complete tree. A root node labeled 'root' should
        be guaranteed to exist by this function's caller.

        Args:
            a_node (dendropy.Node): Node from which to start recurring
                so as to find if there are sampled taxa on both sides
                of root. E.g., the origin node.

        Returns:
            (bool): Whether there is at least one or more sampled
            (observed) taxon on both sides of the complete tree's
            root node.
        """

        def _recur_find_extant_or_sa(
                a_node: dp.Node,
                has_seen_root: bool,
                found_sampled_count: int,
                found_sampled: bool) -> int:

            if a_node.label == "root" or a_node.taxon == "root" \
                    or a_node.label == self.alternative_root_label \
                    or a_node.taxon == self.alternative_root_label:
                has_seen_root = True

            # can only return if this node is not the seed node
            # (origin or root) itself, i.e., we must have recurred
            # once
            #
            # then to save time, we return if we have already found
            # a sampled node on this side of the tree, provided we
            # are not the very first node on this side of the tree!
            # (this is why we have the a_node.parent_node.lavel != 'root')
            if a_node.parent_node and \
                    (a_node.parent_node.label not in ("root", self.alternative_root_label) \
                     and found_sampled):
                return found_sampled_count, found_sampled

            # if we have not returned above, then we need to recur;
            # a_node will have two sides; will do one and then
            # another
            for ch_node in a_node.child_node_iter():
                # if ch_node.label == "root" or ch_node.taxon == "root":
                #     has_seen_root = True

                # if early on, at the root, we found sampled nodes,
                # we just stop recurring 
                if (a_node.label not in ("root", self.alternative_root_label)) \
                        and found_sampled:
                    break

                if has_seen_root:
                    # as soon as a sampled descendant of the root is found,
                    # we are good, we break out of the loop (and return
                    # after it!)
                    if ch_node.sampled or ch_node.is_sa:
                        # debugging
                        # print("      ... this node was sampled! Breaking!")

                        found_sampled_count += 1
                        found_sampled = True

                        # found extant, stop recursion
                        # on current side of a_node if not root
                        if a_node.label not in ("root", self.alternative_root_label):
                            break

                    # still have not found a sampled node,
                    # so stay on this side, keep digging
                    else:
                        # if we're looking at the children of the root,
                        # we must recur no matter what
                        #
                        # the important thing is that we must pass 'False'
                        # as found_sampled because we cannot by definition
                        # have found sampled nodes yet
                        if a_node.label in ("root", self.alternative_root_lable):
                            found_sampled_count, found_sampled = \
                                _recur_find_extant_or_sa(
                                    ch_node,
                                    has_seen_root,
                                    found_sampled_count,
                                    False)

                        # we have seen the root in the past, and the node over
                        # whose children we are iterating is NOT the root, so
                        # this is a standard recursion
                        else:
                            found_sampled_count, \
                                found_sampled = \
                                _recur_find_extant_or_sa(
                                    ch_node,
                                    has_seen_root,
                                    found_sampled_count,
                                    found_sampled)

                else:
                    # no root yet, keep digging until
                    # we find root
                    found_sampled_count, \
                        found_sampled = \
                        _recur_find_extant_or_sa(
                            ch_node,
                            has_seen_root,
                            found_sampled_count,
                            found_sampled)

            # debugging
            # print("end of recurring " + a_node.label \
            #       + " \'found_extant_or_sa_count\' = " \
            #       + str(found_extant_or_sa_count))

            return found_sampled_count, found_sampled

        n_extant_or_sa, _ = _recur_find_extant_or_sa(a_node, False, 0, False)

        # debugging
        # print("\nWrapping up tree, n_extant_or_sa = " + str(n_extant_or_sa))

        if n_extant_or_sa == 2:
            return True

        else:
            return False

    def extract_reconstructed_tree(
            self,
            require_obs_both_sides: ty.Optional[bool] = None) -> dp.Tree:
        """Extract reconstructed tree from complete tree.

        This method was designed to be called outside of
        AnnotatedTree's initialization, and instead be prompted by
        rejection sampling (when conditioning on both sides of the
        root), and summarization and printing functions should the
        user ask for it. This is to save running time.
        
        The method deep-copies self.tree, which is of type
        dendropy.Tree, populates it appropriately and then returns it.

        When populating the reconstructed tree, this method effectively
        prunes extinct taxa, and re-roots the resulting tree at the 
        MRCA of the sampled (observed) taxa. This includes all taxa
        sampled at the present and direct ancestors.

        Args:
            require_obs_both_sides (bool, optional): Flag specifying
                if sampled (observed) taxa are required on both sides
                of the root.

        Returns:
            (dendropy.Tree): A tree instance containing the
                reconstructed tree extracted from AnnotatedTree's
                instance owning the method call.
        """

        # if method called by someone other than Tree obj,
        # then require_obs_both_sides won't be None
        require_obs_both_sides_root = True
        if require_obs_both_sides is not None:
            require_obs_both_sides_root = require_obs_both_sides

        # otherwise, we use the member inside Tree
        else:
            require_obs_both_sides_root = \
                self.condition_on_obs_both_sides_root

        if self.tree_reconstructed:
            return self.tree_reconstructed

        # taxon_name space of tree_reconstructed will not
        # have labels and taxon labels for internal nodes!
        self.tree_reconstructed = copy.deepcopy(self.tree)

        # if tree went extinct, return empty new tree
        if (self.n_extant_terminal_nodes + self.n_sa_nodes) <= 1:
            return dp.Tree()

        # filter_fn = lambda nd: \
        #     nd.is_sa or (nd.is_leaf() and nd.alive and nd.sampled)

        def filter_fn(nd):
            return nd.is_sa or (nd.is_leaf() and nd.alive and nd.sampled)

        self.tree_reconstructed.filter_leaf_nodes(
            filter_fn, suppress_unifurcations=True)

        # getting all root information #
        root_node: dp.Node = dp.Node()
        root_node_distance_from_seed = 0.0
        smallest_distance = 0.0

        root_node = self.tree_reconstructed.find_node_with_label("root")

        # could not get node with label, try with taxon label
        if root_node is None:
            root_node = \
                self.tree_reconstructed.find_node_with_taxon_label("root")
            
        if root_node:
            root_node_distance_from_seed = \
                root_node.distance_from_root()
            smallest_distance = root_node_distance_from_seed

        int_node_deeper_than_root: dp.Node = dp.Node()

        # we will now grab the internal node that is deeper
        # than the root, and also the furthest away from it
        # (internal_nd_distance will be more and more negative)
        for internal_nd in \
            self.tree_reconstructed.preorder_internal_node_iter(
                exclude_seed_node=False):

            internal_nd_distance = internal_nd.distance_from_root()

            # double-checking we're looking at just internal
            # nodes that are not the root
            if internal_nd.label not in ("root", "origin"):
                # we set int_node_deeper_than_root if
                # (i) we haven't looked at any internal nodes yet,
                # i.e., smallest_distance == 0.0,
                #
                # OR
                #
                # (ii) they are deeper than the last one we checked,
                # i.e., internal_nd_distance < smallest_distance
                #
                # AND
                #
                # we also only either
                # (iii) not have a root at all, possibly because tree
                # was read in as newick and nodes were unnamed (or maybe
                # there was no root indeed), and by definition all
                # internal nodes are deeper than the root, i.e., 
                # not root_node
                #
                # OR
                #
                # (iv) have a root, and it is some distance away from
                # the origin
                if (smallest_distance == 0.0 or
                        internal_nd_distance < smallest_distance) and \
                        (root_node_distance_from_seed != 0.0 or
                            not root_node):
                    int_node_deeper_than_root = internal_nd
                    smallest_distance = internal_nd_distance

        ##############
        # Re-rooting #
        ##############

        if not require_obs_both_sides_root:
            # no speciation happened, so complete tree has no root,
            # or tree was read in as a Newick string with unnamed
            # nodes; either way, we re-root at the node we found
            # to be the deepest in the tree except for the origin
            if root_node is None:
                self.tree_reconstructed \
                    .reroot_at_node(int_node_deeper_than_root)
                int_node_deeper_than_root.edge_length = 0.0
            # so there must be a
            # sampled ancestor before the root, so we need to re-root
            # above complete tree's root

            # there is a root in the complete tree, and it may or not
            # be the MRCA of all sampled taxa
            else:
                # if internal nodes in tree do not have a label because
                # they were read in as Newick strings instead of simulated
                # by PJ
                if int_node_deeper_than_root.label == "None":
                    self.tree_reconstructed.reroot_at_node(
                        int_node_deeper_than_root)

                    # getting origin and removing it if possible
                    origin_node_rec: dp.Node = dp.Node()
                    if self.origin_node:
                        origin_node_rec = self.tree_reconstructed. \
                            find_node_with_label("origin")

                        # could not get origin node with label,
                        # try with taxon label
                        if not origin_node_rec:
                            origin_node_rec = self.tree_reconstructed. \
                                find_node_with_taxon_label("origin")

                        # origin should be child because we re-rooted!
                        int_node_deeper_than_root.remove_child(origin_node_rec)

                    int_node_deeper_than_root.edge_length = 0.0

                # nodes in tree do have labels, so they were simulated in PJ;
                # we now find the MRCA of sampled (observed) taxa (it may or
                # not be the complete tree's root)
                else:
                    # DendroPy's mrca seems to have side-effect... ugh!
                    # will do things on deep copy of tree instead, and 
                    # look at labels
                    rec_tree_mrca_label = \
                        pj_get_name_mrca_obs_terminals(
                            self.tree_reconstructed.seed_node,
                            [leaf.label for leaf
                             in self.tree_reconstructed.leaf_node_iter()])

                    # re-root (re-seed in DendroPy) if necessary
                    # if so, it is because there must be a sampled ancestor
                    # before the root, so we need to re-root above complete
                    # tree's root
                    if rec_tree_mrca_label != "root":
                        rec_mrca_node = \
                            self.tree_reconstructed \
                                .find_node_with_label(rec_tree_mrca_label)
                        self.tree_reconstructed.reroot_at_node(rec_mrca_node)

                        rec_mrca_node.edge_length = 0.0

                    # complete tree root is rec tree root
                    else:
                        root_node.edge_length = 0.0

        # conditioning on sampled (observed) taxa on both sides of root
        else:
            # if we condition on sampled nodes on both sides of root,
            # then by definition the reconstructed tree is empty!
            if not root_node:
                return dp.Tree()

            else:
                # if there is a root, but sampled nodes were not found
                # on both sides, reconstructed tree is empty!
                if not self.is_extant_or_sa_on_both_sides_complete_tr_root(
                        root_node):
                    return dp.Tree()

                # there are sampled nodes on both sides of the root!
                #
                # but now we need to see if re-rooting is necessary;
                # if there is an internal node between the origin and
                # the root, then a direct ancestor sampling event took
                # place -- we then root at the dummy node who is the
                # parent of the direct sampled ancestor node (but the 
                # root in reality is the sampled ancestor itself!)
                if self.origin_node and \
                        int_node_deeper_than_root.taxon is not None:
                    self.tree_reconstructed \
                        .reroot_at_node(int_node_deeper_than_root)

                    # we remove the branch subtending the dummy node
                    int_node_deeper_than_root.edge_length = 0.0

        return self.tree_reconstructed
   
    def populate_nd_attr_dict(self,
                              attrs_of_interest_list: ty.List[str],
                              attr_dict_added_separately_from_tree: bool = False) \
                                -> None:
        """Populate member nested dictionary with node attributes.

        This method is not called upon initialization of class,
        but rather when information about nodes' states is required,
        e.g., when get_taxon_states_str() is called. The method
        takes the attribute values (for one or more attributes) stored
        in a DendroPy.Tree, and copies them into the member dictionary
        node_attr_dict.

        There is no return and only a side-effect.

        Args:
            attrs_of_interest_list (str): List of attribute names to
                store in member dictionary.
            attr_dict_added_separately_from_tree (bool): Flag specifying
                whether method is being called outside tree
                initialization (i.e., to a tree that already exists),
                or during tree initialization.
        """

        for nd in self.tree:
            if not nd.label:
                exit(("Must name all nodes before populating nd_attr_dict. "
                      "Exiting"))

            for att, att_val in nd.__dict__.items():
                if attr_dict_added_separately_from_tree:
                    if att == "_annotations":
                        for att_str, att_v in att_val.values_as_dict().items():
                            if att_str in attrs_of_interest_list:
                                try:
                                    self.node_attr_dict[nd.label][att_str] = \
                                        float(att_v)
                                
                                except:
                                    raise ec.AnnotatedTreeNodeMissingAttrError(
                                        nd,
                                        att_str,
                                        "Attribute value could not be converted to float."
                                    )

                else:
                    if att in attrs_of_interest_list:
                        self.node_attr_dict[nd.label][att] = att_val

    def __str__(self) -> str:
        return self.tree.as_string(
            schema="newick",
            suppress_annotations=False,
            suppress_internal_taxon_labels=True)

    def plot_node(self,
                  axes: plt.Axes,
                  node_attr: str = "state",
                  **kwargs) -> None:
        """Draw tree on provided Axes instance.

        This method is required whenever the user asks for a DAG node
        to be drawn.

        Args:
            axes (matplotlib.pyplot.Axes): Axes object where we are
                drawing the tree.
            node_attr (str): Name of the attribute according to which
                one wants to color the AnnotatedTree's branches with.
                Defaults to 'state'.
        """

        if not node_attr:
            plot_ann_tree(self, axes)

        else:
            self.populate_nd_attr_dict([node_attr])

            plot_ann_tree(self, axes, attr_of_interest=node_attr)

    def get_stats_dict(self) -> ty.Dict[str, ty.Union[int, float]]:
        """Get dictionary with AnnotatedTree's summary stats.

        This method is required whenever the user asks for a DAG node
        to be summarized.

        Returns:
            (dict): Dictionary with summary statistic names as keys
            and their values as values.
        """

        ks = ["Origin age",
              "Root age",
              "Total taxon count",
              "Extant taxon count",
              "Extinct taxon count",
              "Direct ancestor count"]
        vs = [self.origin_age,
              self.root_age,
              (self.n_extant_terminal_nodes + self.n_extinct_terminal_nodes),
              self.n_extant_terminal_nodes,
              self.n_extinct_terminal_nodes,
              self.n_sa_nodes]

        return dict((ks[i], str(vs[i])) for i in range(len(ks)))

    def _get_taxon_states_dict(self) -> ty.Tuple[ty.Dict[str, int], ...]:
        """Collect and return non-extinct node states into dicts.

        This method gets the value of attribute 'state' for all
        non-extinct nodes in the tree, and returns terminal and
        internal node states as two separate dictionaries.

        Return:
            (dict): A tuple of dictionaries containing the values of
            the state attribute for terminal (first dictionary) and
            internal nodes (second dictionary).
        """

        self.populate_nd_attr_dict(["state"])

        terminal_node_states_dict: ty.Dict[str, int] = dict()
        int_node_states_dict: ty.Dict[str, int] = dict()
        for taxon_name, attr_val_dict in self.node_attr_dict.items():
            a_node = self.tree.find_node_with_label(taxon_name)

            # internal nodes
            if a_node.is_internal():
                int_node_states_dict[taxon_name] = attr_val_dict["state"]

            # ignoring extinct nodes (after internal nodes have been
            # considered!)
            if not a_node.alive:
                continue

            terminal_node_states_dict[taxon_name] = attr_val_dict["state"]

        return terminal_node_states_dict, int_node_states_dict

    def get_taxon_states_str(self, nexus: bool = False) -> str:
        """Get states for all nodes in tree as single string.
        
        Args:
            nexus (bool): Flag specifying whether states are being
                collected for Nexus printing ('True' if so). Defaults
                to 'False'.

        Returns:
            (str): String containing the states of all nodes in tree.
        """

        terminal_node_state_dict, int_node_states_dict = \
            self._get_taxon_states_dict()

        terminal_node_states_str = ""

        for taxon_name, taxon_state in terminal_node_state_dict.items():
            terminal_node_states_str += \
                taxon_name + "\t" + str(taxon_state) + "\n"

        # only care about internal nodes in the reconstructed tree
        int_node_states_str = ""

        for taxon_name, taxon_state in int_node_states_dict.items():
            if self.tree_reconstructed.__str__() != "" and \
                    self.tree_reconstructed is not None:
                # note how maybe 'root' may not figure among the internal
                # nodes, because it may not be in the reconstructed tree
                int_node = \
                    self.tree_reconstructed.find_node_with_label(taxon_name)

                if int_node is not None:
                    int_node_states_str += \
                        taxon_name + "\t" + str(taxon_state) + "\n"

        # not sure if I should be adding self.n_sa_nodes here, because
        # I ignore dead taxa in _get_taxon_states_dict()
        if nexus:
            nexus_header = \
                "#Nexus\n\nBegin data;\nDimensions ntax=" + \
                str(self.n_extant_terminal_nodes + self.n_sa_nodes) + \
                " nchar=1;\nFormat datatype=Standard symbols=\"" + \
                "".join(str(i) for i in range(self.state_count)) + \
                "\" missing=? gap=-;\nMatrix\n"

            nexus_str = nexus_header + terminal_node_states_str + ";\nEnd;\n"

            return nexus_str

        return terminal_node_states_str, int_node_states_str


###########################
# Plotting tree functions #
###########################

def get_node_name(nd: dp.Node) -> str:
    """Return node name.

    Get name of node that may be stored in a DendroPy.Node as either
    its label, taxon, or taxon.label members. It tries them all to
    make sure a node name is found and returned. If none is found,
    an exception is raised.

    Args:
        nd (dendropy.Node): Node whose name one wants to collect.

    Returns:
        (str): Name of the node.
    """
    
    if nd.label:
        return nd.label

    elif nd.taxon and isinstance(nd.taxon, str):
        return nd.taxon

    elif nd.taxon and nd.taxon.label:
        return nd.taxon.label
    
    else:
        raise ec.AnnotatedTreeMissingNodeName()


def get_x_coord_from_nd_heights(ann_tr: AnnotatedTree,
                                use_age: bool = False,
                                unit_branch_lengths: bool = False) \
                                    -> ty.Dict[str, float]:
    """Return x-coordinates for all nodes in tree.
    
    This method returns a dictionary of node labels as keys,
    node x_coords (time) as values.

    Args:
        ann_tr (AnnotatedTree): Instance of AnnotatedTree that we
            are drawing.
        use_age (bool): Flag specifying if to use node age or not.
            Defaults to 'False'.
        unit_branch_lengths (bool): If branch lengths are all 1.0
            (currently not used). Defaults to 'False'.
    """

    if use_age and not unit_branch_lengths:
        return ann_tr.node_ages_dict

    return ann_tr.node_heights_dict


def get_y_coord_from_n_obs_nodes(ann_tr: AnnotatedTree,
                                 start_at_origin: bool = False,
                                 sa_along_branches: bool = True) \
                                    -> ty.Dict[str, float]:
    """Return y-coordinates for all nodes in tree.
    
    This method returns a dictionary of node labels as keys,
    y-coords as values. Y-coords here are integers that go from 1 to
    the total number of observable nodes. Every observable node will be
    1 y-unit away from each other.

    Args:
        ann_tr (AnnotatedTree): Instance of AnnotatedTree that we
            are drawing.
        start_at_origin (bool): Flag specifying if tree starts at
            the origin node. Defaults to 'False'.
        sa_along_branches (bool): Flag specifying if direct
            (sampled) ancestors should be placed along a branch.
            Defaults to 'True'.
    """

    # we define an inner recursive function
    # for obtaining y-coords at internal nodes
    def recursively_calculate_height(nd: dp.Node) -> None:
        children: ty.List[dp.Node] = list()

        if sa_along_branches:
            children = [ch for ch in nd.child_nodes() if not ch.is_sa]

        else:
            children = [ch for ch in nd.child_nodes()]

        if len(children) > 0:
            for ch_nd in children:
                if ch_nd not in y_coords:
                    recursively_calculate_height(ch_nd)

            # origin (and dummy nodes if user wants to place SA nodes along branches)
            if (start_at_origin and nd == ann_tr.origin_node) \
                    or (sa_along_branches and nd.is_sa_dummy_parent):
                y_coords[get_node_name(nd)] = \
                    y_coords[get_node_name(children[0])]

            # root and all others (assume bifurcation)
            else:
                y_coords[get_node_name(nd)] = \
                    (y_coords[get_node_name(children[0])]
                     + y_coords[get_node_name(children[1])]) / 2.0

    # now we will get the y-coords
    #
    # we start by obtaining the max height of graph,
    # which is given by number of observable nodes
    maxheight: int = 0
    if sa_along_branches:
        maxheight = ann_tr.n_extant_terminal_nodes \
            + ann_tr.n_extinct_terminal_nodes

    else:
        maxheight = ann_tr.n_extant_terminal_nodes \
            + ann_tr.n_extinct_terminal_nodes + ann_tr.n_sa_nodes

    # grab all terminal node names
    leaf_names = list()
    for nd in ann_tr.tree.leaf_node_iter():
        # if user wants to place SA nodes on branches,
        # we must ignore SA nodes here
        if sa_along_branches and nd.is_sa:
            continue

        leaf_names.append(get_node_name(nd))

    # initialize the y-coords with the values for terminal nodes
    y_coords: ty.Dict[str, float] = \
        {leaf_name: float(maxheight - i + 1)
         for i, leaf_name
         in enumerate(reversed(leaf_names))}
    
    # now take care of all internal nodes recursively
    recursively_calculate_height(ann_tr.tree.seed_node)

    return y_coords


def plot_ann_tree(ann_tr: AnnotatedTree,
                  axes: plt.Axes,
                  use_age: bool = False,
                  start_at_origin: bool = False,
                  attr_of_interest: str = "state",
                  sa_along_branches: bool = True) -> None:
    """Plot instance of AnnotatedTree on provided Axes instance.

    Plotting is a side-effect.

    Args:
        ann_tr (AnnotatedTree): Instance of AnnotatedTree that we
            are drawing.
        axes (matplotlib.pyplot.Axes): Axes object where we are
            drawing the tree.
        use_age (bool): Flag specifying if age or time is being
            used.
        start_at_origin (bool): Flag specifying if drawing starts
            at the origin node.
        attr_of_interest (str): Name of the attribute according to
            which states we are coloring tree branches. Defaults to
            'state'.
        sa_along_branches (bool): Flag specifying if direct
            (sampled) ancestors should be placed along a branch.
            Defaults to 'True'.
    """

    color_map: ty.Dict[int, str]
    attr_found: bool = True
    if attr_of_interest not in \
        ann_tr.node_attr_dict[ann_tr.tree.seed_node.label]:
        color_map = get_color_map(1)
        attr_found = False

    else:
        color_map = get_color_map(ann_tr.state_count)

    ####################################################
    # Setting up flags and checking tree content is OK #
    ####################################################
    if ann_tr.origin_node:
        start_at_origin = True  # scoped to draw_tree()

    if sa_along_branches and not ann_tr.sa_lineage_dict:
        # throw plotting exception, cannot plot sa along branches
        # without dict
        pass

    # grabbing key coordinates for drawing tree
    x_coords = get_x_coord_from_nd_heights(ann_tr,
                                           use_age=use_age)
    y_coords = \
        get_y_coord_from_n_obs_nodes(ann_tr,
                                     start_at_origin=start_at_origin,
                                     sa_along_branches=sa_along_branches)

    # arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # make sure we reset things
    axes.cla()

    # side-effect:
    # populates: horizontal_linecollections
    #            vertical_linecollections
    def _draw_clade_lines(x_end: float = 0.0,
                          use_linecollection: bool = False,
                          orientation: str = "horizontal",
                          y_here: int = 0,
                          x_start: int = 0,
                          x_here: int = 0,
                          y_bot: int = 0,
                          y_top: int = 0,
                          color: str = "black",
                          lw: float = ".1") -> None:
        """Draw a hor. or vert. line or add it to collection.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """

        if not use_linecollection:
            if orientation == "horizontal":
                axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
            
            elif orientation == "vertical":
                axes.vlines(x_here, y_bot, y_top, color=color)

        elif use_linecollection:
            if orientation == "horizontal":
                if not x_end:
                    horizontal_linecollections.append(
                        mpcollections.LineCollection(
                            [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                        )
                    )
                else:
                    horizontal_linecollections.append(
                        mpcollections.LineCollection(
                            [[(x_end, y_here), (x_here, y_here)]], color=color, lw=lw
                        )
                    )

            elif use_linecollection and orientation == "vertical":
                vertical_linecollections.append(
                    mpcollections.LineCollection(
                        [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                    )
                )

    def _draw_clade(nd: dp.Node,
                    x_start: int,
                    color: str,
                    lw: float,
                    use_age: bool,
                    start_at_origin: bool,
                    keep_it_black: bool = False) -> None:
        """Recursively draw a tree.

        This is the main drawing method.
        
        Args:
            nd (dendropy.Node): Current node being recurred and whose
                lines are being drawn.
            x_start (int): Position (integer) of node 'nd' on x-axis.
            color (str): Name of color associated to the state of the
                node parent (!) to focal node 'nd'. We pass the color
                of the parent node because branches can undergo
                attribute state changes (color changes), and the 'old'-
                end of the branch should have the color of the parent
                node.
            lw (float): Line width.
            use_age (bool): Flag specifying if age or time is being
                used.
            start_at_origin (bool): Flag specifying if drawing starts
                at the origin node.
        """

        nd_name = get_node_name(nd)
        attr_idx = 0
        
        # color of node parent to 'nd'!
        segment_colors = [color]

        # we override it if only one state to attribute
        if keep_it_black:
            segment_colors = ["black"]

        # in case nodes have attributes, there may be attribute state
        # changes, so the color of the branch will not necessarily
        # be the same of the node parent to 'nd' anymore
        #
        # sowe we update segment_colors with have the color matching
        # the state of nd (the node subtending the branch)
        if ann_tr.node_attr_dict and attr_of_interest is not None and \
            attr_of_interest in ann_tr.node_attr_dict[nd_name] and \
                not keep_it_black:
                    attr_idx = ann_tr.node_attr_dict[nd_name][attr_of_interest]
                    segment_colors = [color_map[attr_idx]]

        x_starts = [x_start]
        x_heres = []
        x_here_int_nodes = x_coords[nd_name]

        if ann_tr.at_dict is not None and not keep_it_black:
            # if branch subtending nd (i.e., branch starting at nd and going
            # into the past) underwent an attribute transition
            if nd_name in ann_tr.at_dict:
                # we get the list of attribute transitions (which goes from
                # old to young, i.e., from low time to high time)
                attr_trs = ann_tr.at_dict[nd_name]

                # iterating over attr transitions, from old to young
                for idx, at in enumerate(attr_trs):
                    # x_starts will go from low numbers to high numbers,
                    # i.e., old to young
                    x_starts.append(at.global_time)
                    x_heres.append(at.global_time)

                    # x_starts.append(attr_trs[-(idx+1)].global_time)
                    # x_heres.append(attr_trs[-(idx+1)].global_time)
                    
                    # when segment_colors is initialized, its first entry
                    # is the color of nd (which is the node subtending (!)
                    # the branch along which attribute transitions have
                    # happened
                    #
                    # this means that the first entry of segment_colors
                    # corresponds to the youngest state (actually, its color)
                    # 
                    # but because segment_colors has to match x_starts, we
                    # need to orient segment_colros so that it has older
                    # states (color) first, and then younget states
                    #
                    # this is why we reverse the indexing and add 1
                    segment_colors.insert(0, color_map[attr_trs[-(idx+1)].from_state])

            else:
                segment_colors.append(color)

            # age of focal node is always
            # the latest time (x_here)
            x_heres.append(x_coords[nd_name])

        # if no attribute transition dict, we have a single line with a single
        # color to draw (default behavior)
        else:
            x_heres = [x_coords[nd_name]]
            segment_colors = ["black"]

        y_here = y_coords[nd_name]

        from_x_here_to_this_x: float = 0.0  # default (present moment)
        # draws from origin/root_node to root_node x
        if use_age and nd in (ann_tr.origin_node, ann_tr.root_node):
            if isinstance(ann_tr.root_age, float):
                from_x_here_to_this_x = ann_tr.root_age

        #################################
        # Draw horizontal line (branch) #
        #################################

        # debugging
        # print("\nabout to draw hor. line for node " + nd_name)
        # print("where x_starts:")
        # print(x_starts)
        # print("where x_heres:")
        # print(x_heres)

        # drawing clade lines on the left (old) first,
        # then on the right (young)
        for idx in range(len(x_starts)):
            _draw_clade_lines(
                x_end=from_x_here_to_this_x,
                use_linecollection=True,
                orientation="horizontal",
                y_here=y_here,
                x_start=x_starts[idx],
                x_here=x_heres[idx],
                color=segment_colors[idx],
                lw=lw
            )

        #########################
        # Add node/taxon labels #
        #########################
        if nd.is_leaf():
            axes.text(x_heres[-1],
                      y_here,
                      f" {nd_name}",
                      verticalalignment="center")

        #################################
        # Draw sampled ancestors if any #
        #################################
        if sa_along_branches and \
                (not nd.is_sa_dummy_parent and nd.is_sa_lineage):

            sas: ty.List[pjsa.SampledAncestor] = []
            if isinstance(ann_tr.sa_lineage_dict, dict):
                sas = ann_tr.sa_lineage_dict[nd_name]

            for sa in sas:
                sa_x: float = sa.global_time
                sa_color: str = "k"

                if not keep_it_black:
                    color_map[sa.state]
                # sa_x: float = x_here - sa.time_to_lineage_node # works fine

                if use_age:
                    if start_at_origin and \
                            isinstance(ann_tr.origin_age, float):
                        sa_x = ann_tr.origin_age - sa_x
                    elif isinstance(ann_tr.root_age, float):
                        sa_x = ann_tr.root_age - sa_x

                sa_y = y_here

                _draw_sa_along_branch(sa.label,
                                      sa.global_time,
                                      sa_y,
                                      axes,
                                      color=sa_color)

        # recur if not leaf
        if len(nd.child_nodes()) > 0:
            children = nd.child_nodes()

            ################################################
            # Draw a vertical line connecting two children #
            ################################################
            # not origin
            if not start_at_origin or \
                (ann_tr.origin_node and nd != ann_tr.origin_node):

                # must be either (1) placing SA along branch, and then NOT a dummy, or
                # or             (2) regular plotting, AND have two children
                if (sa_along_branches and not nd.is_sa_dummy_parent) or \
                        (len(children) == 2 and not sa_along_branches):

                    y_top = y_coords[get_node_name(children[0])]
                    y_bot = y_coords[get_node_name(children[1])]

                    # last color in segment_colors will
                    # match the state of the node whose
                    # subtending branch we are drawing
                    _draw_clade_lines(
                        use_linecollection=True,
                        orientation="vertical",
                        x_here=x_here_int_nodes,
                        y_bot=y_bot,
                        y_top=y_top,
                        # color=segment_colors[0],
                        color=segment_colors[-1],
                        lw=lw,
                    )

            ####################
            # Draw descendents #
            ####################
            for child_nd in children:
                # if user wants to place SA nodes along branches, we must
                # ignore them here
                if sa_along_branches and child_nd.is_sa:
                    continue

                # segment_colors[-1] is the color of the (parent!) node whose
                # children we are visiting
                _draw_clade(
                    child_nd,
                    x_heres[-1],
                    segment_colors[-1],
                    lw,
                    False,
                    start_at_origin,
                    keep_it_black=keep_it_black)

    def _draw_time_slices(
            ann_tr: AnnotatedTree,
            axes: plt.Axes,
            use_age: bool) -> None:
        """Draw vertical lines representing boundaries of time slices.
        
        Args:
            ann_tr: Instance of AnnotatedTree that we
                are drawing.
            axes (matplotlib.pyplot.Axes): Axes object where we are
                drawing the tree.
            use_age (bool): Flag specifying if age or time is being
                used.
        """

        xs: ty.List[float] = []
        if use_age and ann_tr.slice_age_ends:
            xs = ann_tr.slice_age_ends[1:]  # ignore present

        elif ann_tr.slice_t_ends is not None:
            xs = [t_end for t_end in ann_tr.slice_t_ends[:-1]
                  if isinstance(t_end, float)]  # so mypy won't complain

        for x in xs:
            axes.axvline(x=x, color="gray", dashes=[2.0, 1.5])

    def _draw_sa_along_branch(
            sa_name: str,
            x: float,
            y: float,
            axes: plt.Axes,
            color: str = "k") -> None:
        """Draw direct (sampled) ancestor nodes on tree.
        
        This method places a dot on the given coordinate.

        Args:
            sa_name (str): Name of the direct (sampled) ancestor.
            x (float): X-axis coordinate where to draw node.
            y (float): Y-axis coordinate where to draw node.
            axes (matplotlib.pyplot.Axes): Axes object where we are
                drawing the tree.
            color (str): Color of the node being drawn. Defaults to
                'k'.
        """

        # zorder places dot on top of axes
        # k = black color, o = means marker
        axes.plot(x,
                  y,
                  marker="o",
                  color=color,
                  markersize=4,
                  zorder=10)
        axes.text(x,
                  y - 0.05,
                  f" {sa_name}",
                  rotation=45,
                  horizontalalignment="center",
                  fontsize=10)

    ########
    # call #
    ########
    color = "black"

    # if the attribute of interest is 'state' and
    # there is only one state, we keep it all black
    keep_it_black = False
    if attr_of_interest == "state" and ann_tr.state_count == 1:
        keep_it_black = True

    # grabbing the color at the start of the process
    if attr_found:
        attr_idx = ann_tr.node_attr_dict[ann_tr.tree.seed_node.label] \
            [attr_of_interest]
        color = color_map[attr_idx]

    # note that in this first call, the color being passed
    # is the color of the node in the first argument
    #
    # in the recurring calls inside, the color passed to
    # _draw_clade is the color of the parent node to the
    # one passed in the first argument; see comments inside
    # _draw_clade and its docstring
    _draw_clade(ann_tr.tree.seed_node,
                0,
                color,
                plt.rcParams["lines.linewidth"],
                use_age,
                start_at_origin,
                keep_it_black=keep_it_black)

    _draw_time_slices(ann_tr, axes, use_age)

    # draw lines
    for i in horizontal_linecollections:
        axes.add_collection(i)

    for i in vertical_linecollections:
        axes.add_collection(i)

    if use_age:
        axes.set_xlabel("Age")

    else:
        axes.set_xlabel("Time")

    #####################
    # final adjustments #
    #####################

    # add back bottom axis
    axes.spines['bottom'].set_visible(True)

    # add margins around the tree to prevent overlapping the axes
    xmax = max(x_coords.values())
    axes.set_xlim(-0.05 * xmax, 1.05 * xmax)

    # also invert the y-axis (origin at the top)
    # add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_coords.values()) + 0.8, 0.2)

    if use_age:
        axes.invert_xaxis()

    axes.yaxis.set_ticks([])
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)


def get_color_map(n_states: int) -> ty.Dict[int, str]:
    """Create and return a map from discrete state to color.

    This method uses a palette of very contrasting colors if the number
    of states is lesser or equal to 20. If greater than 20 and lesser
    or equal to 120, the method switches to another palette and
    carries out some truncation to avoid colors that are almost white.
    If the number of states is greater than 120, it switches to yet
    another palette, truncating it again. The choice of palette and
    truncation is arbitrary, and was informed by an experiment in
    plotting.pj_seeing_colors.py.

    Args:
        n_states (int): Number of states.

    Returns:
        (dict): A dictionary with an integer representing a
            discrete state, and a string containing the hex
            code of a color.
    """

    def truncate_colormap(cmap,
                          minval: float = 0.0,
                          maxval: float = 1.0,
                          n: int = 100) -> colors.LinearSegmentedColormap:
        """
        Truncate color map (palette) to remove white-ish colors.
        """
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n))
        )

        return new_cmap

    color_map: ty.Dict[int, str] = dict()

    if n_states <= 20:
        qual_cmap = matplotlib.colormaps.get_cmap('tab20')
        color_map = dict(
            (i, matplotlib.colors.rgb2hex(qual_cmap(i)[:3]))
            for i in range(qual_cmap.N))

    elif n_states <= 120:
        # cmap = matplotlib.colormaps.get_cmap('terrain', n_states)
        cmap = plt.get_cmap('terrain', n_states)
        new_cmap = truncate_colormap(cmap, minval=0.0, maxval=2.08, n=n_states)
        color_map = dict(
            (i, matplotlib.colors.rgb2hex(new_cmap(i)[:3]))
            for i in range(new_cmap.N)
        )

    # if n_states >120 up to 250
    else:
        # cmap = matplotlib.colormaps.get_cmap('terrain', n_states)
        cmap = plt.get_cmap('terrain', n_states)
        maxv = 120 * 2.08 / n_states
        new_cmap = truncate_colormap(cmap, minval=0.0, maxval=maxv, n=n_states)
        color_map = dict(
            (i, matplotlib.colors.rgb2hex(new_cmap(i)[:3]))
            for i in range(new_cmap.N)
        )

    return color_map

##########################
# General tree functions #
##########################

def pj_get_name_mrca_obs_terminals(nd: dp.Node,
                                   nd_label_list: ty.List[str]) -> str:
    """Get name of the MRCA of specified observed terminal nodes.

    This method recursively finds the most recent common ancestor
    of the terminal nodes whose names are specified as input.
    These nodes must be observed, i.e., be either direct (sampled)
    ancestors, or sampled extant terminal nodes.

    Args:
        nd (dendropy.Node): Node to recur and grab the name of.
        nd_label_list (str): List of node names whose MRCA's name
            is being searched.

    Returns:
        (str): Name of MRCA node
    """

    # side-effect recursion
    def recur_node(nd, nd_label_list: ty.List[str], mrca_node_label: str):
        """Populate visited_obs_terminals (side-effect)"""

        # not done: if tip, get label
        if nd.is_leaf() and (nd.is_sa or nd.sampled):
            return [nd.label], ""

        # not done: if internal node, we recur
        else:
            visited_obs_terminals: ty.List[str] = []

            for ch_node in nd.child_node_iter():
                # recur
                vot, mrca_node_label = recur_node(ch_node, nd_label_list, mrca_node_label)
                visited_obs_terminals += vot

            # we are actually done
            if set(visited_obs_terminals) == set(nd_label_list) and not mrca_node_label:
                mrca_node_label = nd.label

        return visited_obs_terminals, mrca_node_label

    mrca_node_label: str = ""

    # mrca node
    _, mrca_node_label = recur_node(nd, nd_label_list, mrca_node_label)

    return mrca_node_label


if __name__ == "__main__":
    # Assuming you opened the PhyloJunction/ (repo root) folder
    # on VSCode and that you want to run this as a standalone script,
    # i.e., "Run Without Debugging", you will need to configure your
    # launch.json to have:
    #
    # "env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
    #
    # and your settings.json to have:
    #
    # "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
    #
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable, e.g., [...]/Phylojunction/src/ should be
    # in PYTHONPATH), or to set it if it does not exist -- don't forget to export it!
    #
    # Then you can do:
    # $ python3 src/phylojunction/data/tree.py

    ##############################################
    # Below we have an environment for debugging #
    ##############################################
    origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
    origin_node.state = 0
    origin_node.alive = False
    origin_node.sampled = False
    origin_node.is_sa = False
    origin_node.is_sa_dummy_parent = False
    origin_node.is_sa_lineage = False

    dummy_node = dp.Node(taxon=dp.Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
    dummy_node.state = 0
    dummy_node.alive = False
    dummy_node.sampled = False
    dummy_node.is_sa = False
    dummy_node.is_sa_dummy_parent = True
    dummy_node.is_sa_lineage = False

    origin_node.add_child(dummy_node)

    # right child of dummy_node
    sa_node = dp.Node(taxon=dp.Taxon(label="sa1"), label="sa1", edge_length=0.0)
    sa_node.state = 0
    sa_node.alive = False
    sa_node.sampled = True
    sa_node.is_sa = True
    sa_node.is_sa_dummy_parent = False
    sa_node.is_sa_lineage = False

    # left child of dummy node
    root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.5)
    root_node.state = 1
    root_node.alive = False
    root_node.sampled = False
    root_node.is_sa = False
    root_node.is_sa_dummy_parent = False
    root_node.is_sa_lineage = True

    dummy_node.add_child(sa_node)
    dummy_node.add_child(root_node)

    # left child of root node
    extant_sp1 = dp.Node(taxon=dp.Taxon(label="sp1"), label="sp1", edge_length=0.25)
    extant_sp1.state = 2
    extant_sp1.alive = False
    extant_sp1.sampled = False
    extant_sp1.is_sa = False
    extant_sp1.is_sa_dummy_parent = False
    extant_sp1.is_sa_lineage = False

    # right child of root node
    extant_sp2 = dp.Node(taxon=dp.Taxon(label="sp2"), label="sp2", edge_length=0.5)
    extant_sp2.state = 3
    extant_sp2.alive = True
    extant_sp2.sampled = True
    extant_sp2.is_sa = False
    extant_sp2.is_sa_dummy_parent = False
    extant_sp2.is_sa_lineage = False

    root_node.add_child(extant_sp1)
    root_node.add_child(extant_sp2)

    tr_sa_with_root_survives = dp.Tree(seed_node=origin_node)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(origin_node.taxon)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(sa_node.taxon)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(dummy_node.taxon)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(root_node.taxon)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp1.taxon)
    tr_sa_with_root_survives.taxon_namespace.add_taxon(extant_sp2.taxon)

    # debugging
    # print("tr_sa_with_root_survives.seed_age = " + str(tr_sa_with_root_survives.max_distance_from_root()))
    # print(tr_sa_with_root_survives.as_string(schema="newick"))

    total_state_count = 4

    sa_global_time = 1.0
    time_to_sa_lineage_node = 0.5
    sa = pjsa.SampledAncestor(
        "sa1",
        "root",
        sa_global_time,
        state=0,
        time_to_lineage_node=time_to_sa_lineage_node)
    sa_lineage_dict = {"root": [sa]}

    at1 = pjat.AttributeTransition("state", "root", 1.25, 0, 1)
    at2 = pjat.AttributeTransition("state", "sp1", 1.6, 1, 2)
    at3 = pjat.AttributeTransition("state", "sp2", 1.8, 1, 3)
    at_dict = {
        "root": [at1],
        "sp1": [at2],
        "sp2": [at3]
    }

    max_age = 2.0

    ann_tr_sa_with_root_survives_max_age = \
        AnnotatedTree(
            tr_sa_with_root_survives,
            total_state_count,
            start_at_origin=True,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            at_dict=at_dict,
            epsilon=1e-12
        )

    print(ann_tr_sa_with_root_survives_max_age.tree.as_string(schema="newick"))  
    
    # this print statement invokes populate_nd_attr_dict() 
    # that then populates a member of AnnotatedTree, allowing
    # for branches to be painted according to the right attribute
    print(ann_tr_sa_with_root_survives_max_age.get_taxon_states_str(nexus=True))

    ######################
    # Preparing plotting #
    ######################
    # fig = Figure(figsize=(11,4.5))

    # note that pjgui uses matplotlib.figure.Figure
    # (which is part of Matplotlib's OOP class library)
    #
    # here, we instead use pyplot's figure, which is the
    # Matlab-like state-machine API
    fig = matplotlib.pyplot.figure()

    ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    ax.patch.set_alpha(0.0)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plot_ann_tree(ann_tr_sa_with_root_survives_max_age,
                  ax,
                  use_age=False,
                  start_at_origin=True,
                  sa_along_branches=True,
                  attr_of_interest="state")
    # sa_along_branches=False

    plt.show()
