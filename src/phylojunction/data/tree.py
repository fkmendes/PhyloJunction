from cgitb import small
import typing as ty
import dendropy as dp # type: ignore
import matplotlib # type: ignore
import copy # type: ignore
import numpy as np
import collections
# from dendropy import Node, Tree, Taxon 

# plotting tree
import matplotlib.pyplot as plt # type: ignore
import matplotlib.collections as mpcollections # type: ignore
from matplotlib.figure import Figure # type: ignore

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
import phylojunction.data.sampled_ancestor as pjsa

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class AnnotatedTree(dp.Tree):

    with_origin: bool
    origin_node: ty.Optional[dp.Node] # can be None
    origin_age: ty.Optional[float] # can be None
    root_node: ty.Optional[dp.Node] # can be None
    brosc_node: ty.Optional[dp.Node] # can be None
    root_age: ty.Optional[float] # can be None
    tree_read_as_newick_by_dendropy: bool
    tree: dp.Tree
    state_count: int
    seed_age: float
    max_age: ty.Optional[float]
    epsilon: float
    tree_died: ty.Optional[bool]
    tree_invalid: ty.Optional[bool]
    no_event: bool
    state_count_dict: ty.Dict[int, int]
    node_heights_dict: ty.Dict[str, float]
    node_ages_dict: ty.Dict[str, float]
    node_attr_dict: ty.Dict[str, ty.Dict[str, ty.Any]] # { node_name: { attr: value }}
    slice_t_ends: ty.List[ty.Optional[float]]
    slice_age_ends: ty.Optional[ty.List[float]] # can be None
    sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]] # can be None
    n_extant_obs_nodes: int
    n_extinct_obs_nodes: int
    n_sa_obs_nodes: int
    extant_obs_nodes_labels: ty.Tuple[str, ...]
    extinct_obs_nodes_labels: ty.Tuple[str, ...]
    sa_obs_nodes_labels: ty.Tuple[str, ...]

    # rec tree
    tree_reconstructed: dp.Tree
    condition_on_obs_both_sides_root: bool

    def __init__(self,
                a_tree: dp.Tree,
                total_state_count: int,
                start_at_origin: bool=False,
                condition_on_obs_both_sides_root: bool=True,
                max_age: ty.Optional[float]=None,
                slice_t_ends: ty.List[ty.Optional[float]]=[],
                slice_age_ends: ty.Optional[ty.List[float]]=None,
                sa_lineage_dict: ty.Optional[ty.Dict[str, ty.List[pjsa.SampledAncestor]]]=None,
                tree_died: ty.Optional[bool]=None,
                tree_invalid: ty.Optional[bool]=None,
                epsilon: float=1e-12):

        self.origin_node = None
        self.root_node = None
        self.brosc_node = None
        self.tree_read_as_newick_by_dendropy = False
        self.no_event = True

        # some initial tree specs
        self.tree = a_tree
        self.tree_reconstructed = None
        self.with_origin = start_at_origin
        self.tree_died = tree_died

        # state-related
        self.state_count = total_state_count
        self.state_count_dict = dict((int(s), 0) for s in range(self.state_count))

        # age related
        self.seed_age = self.tree.max_distance_from_root()
        self.max_age = max_age
        self.node_heights_dict: ty.Dict[str, float] = dict()
        self.node_ages_dict: ty.Dict[str, float] = dict()
        # self.node_attr_dict = pjh.autovivify(2) # populated on demand by populate_nd_attr_dict, cannot be pickled for some reason...
        self.node_attr_dict = collections.defaultdict(pjh.create_str_defaultdict) # this one can be pickled
        self.slice_t_ends = slice_t_ends
        self.slice_age_ends = slice_age_ends
        
        # reconstructed tree
        self.condition_on_obs_both_sides_root = condition_on_obs_both_sides_root

        # rejection sampling when stop condition
        # is age
        self.tree_invalid = tree_invalid
        if not isinstance(self.tree_invalid, bool):
            # (if flag not passed, we assume tree is valid)
            self.tree_invalid = False

        # other
        self.epsilon = epsilon
        self.sa_lineage_dict = sa_lineage_dict

        try:
            self.tree.seed_node.alive
        except:
            self.tree_read_as_newick_by_dendropy = True

        ########### Initializing important class members ############
        # (1) origin_age                                            #
        # (2) root_age (if there is one)                            #
        # (3) origin_edge (the edge that has the origin on one end) #
        # (4) tree_died (flag signaling all tips went extinct)      #
        #############################################################
        if self.with_origin:
            # making sure this is indeed origin...
            if len(self.tree.seed_node.child_nodes()) > 1:
                raise ec.AnnotatedTreeMisspec("Tree was specified as starting from origin, but it started from root. Exiting...")

            self.origin_node = self.tree.seed_node
            self.root_age = 0.0

            origin_children = [nd for nd in self.tree.nodes() if nd.label != "origin"]
            origin_children_labels = tuple([nd.label for nd in origin_children])
            
            # if there is a brosc node, we get it
            try: self.brosc_node = [nd for nd in origin_children if nd.label == "brosc"][0]
            except: pass # self.brosc_node will remain None

            # Case (a): no events took place, tree may or may not have gone extinct
            #
            # [origin] --- [brosc] ... [max age] (died)
            # [origin] ------ [brosc at max age] (survived)
            if self.origin_node.num_child_nodes() == 1 and origin_children_labels[0] == "brosc":
                self.origin_age = self.origin_edge_length = origin_children[0].edge_length # will have been set at stop_condition check in simulate()
                # print("Without root:")
                # print("  origin_edge_length = " + str(self.origin_edge_length))
                # print("  root_age = " + str(self.root_age))

            # Case (b): at least one event took place (root may or not have been born)
            # and tree may or not have died
            elif len(origin_children) > 1:
                self.no_event = False
                
                # Case (b.1) There is a root
                #
                # [origin] --- [root + children] ............. [stop condition] (died)
                # [origin] ---------------- [root + children at stop condition] (survived)
                if "root" in origin_children_labels:
                    self.origin_age = self.seed_age
                    self.root_node = [nd for nd in origin_children if nd.label == "root"][0] # there might be sampled ancestors before root, making sure...
                    self.origin_edge_length = self.recursively_find_node_age(self.root_node, 0.0)
                    self.root_age = self.origin_age - self.origin_edge_length # returns 0.0 if tree dies (TODO: revisit this...)
                    # print("With root:")
                    # print("  origin_edge_length = " + str(self.origin_edge_length))
                    # print("  root_age = " + str(self.root_age))
                
                # Case (b.2) There is no root, so event(s) must have been ancestor sampling
                #
                # [origin] --- [ancestor sampling(s)] --- [brosc] ... [max age] (died)
                # [origin] --- [ancestor sampling(s)] ------ [brosc at max age] (survived)
                else:
                    for nd in self.tree.seed_node.leaf_nodes():
                        if not nd.is_sa and abs(nd.edge_length) > self.epsilon:
                            self.origin_age = self.recursively_find_node_age(nd, 0.0)
                            
                    self.origin_edge_length = self.recursively_find_node_age(self.brosc_node, 0.0)
                    # print("Without root:")
                    # print("  origin_edge_length = " + str(self.origin_edge_length))
                    # print("  root_age = " + str(self.root_age))

            else:
                print("\n\nShould not be here\n\n")
                print(self.tree.as_string(schema="newick"))
                # print("origin children = ")
                # print(origin_children)

            # debugging
            # print("\nFinished initializing some class members:")
            # print("origin_age = " + str(self.origin_age))
            # print("origin_edge_length = " + str(self.origin_edge_length))
            # print("root_age = " + str(self.root_age))

            ##################################################
            # Figuring out if tree died, when tree_died flag #
            # was not passed upon initialization             #
            ##################################################
            if not isinstance(self.tree_died, bool):
                self.tree_died = True

                # with max_age
                if isinstance(max_age, float) and isinstance(self.origin_age, float) and \
                    (max_age - self.origin_age) <= self.epsilon:
                    self.tree_died = False
            
                # no max_age
                else:
                    # but there is a brosc_node
                    if not self.brosc_node == None:
                        # who is alive
                        if isinstance(self.brosc_node, dp.Node) and self.brosc_node.alive:
                            self.tree_died = False
                        # who is dead
                        else:
                            self.tree_died = True

                    # no brosc_node
                    elif self.brosc_node == None:
                        try:
                            self.tree_died = self.has_tree_died(self.origin_node)
                        # will get here if node doesn't have .alive
                        except:
                            pass
                    
                    # no way to tell, we assume tree died
                    # else:
                    #     self.tree_died = True

            # a bit of cleaning for printing the tree
            if self.tree_read_as_newick_by_dendropy and len(self.tree.leaf_nodes()) == 1 and\
                isinstance(self.root_node, dp.Node):
                self.root_node.label = str(self.root_node.taxon).replace("\'", "")

        # starting at root
        else:
            if len(self.tree.seed_node.child_nodes()) == 1:
                raise ec.AnnotatedTreeMisspec("Tree was specified as not starting from origin, but it starts from the origin. Exiting...")
            # if tree is read as a Newick string, the user forgot to specify "start_at_origin=True".
            # and there is a root edge, we need to add an origin node and update everything accordingly
            if len(self.tree.seed_node.child_nodes()) == 2 and self.tree.seed_node.edge_length:
                self.root_node = self.tree.seed_node
                self.origin_edge_length = self.root_node.edge_length
                self.origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
                self.origin_node.alive = False
                self.origin_node.add_child(self.root_node)
                self.tree.seed_node = self.origin_node
                self.seed_age = self.tree.max_distance_from_root()
                self.origin_age = self.seed_age
                self.root_age = self.origin_age - self.origin_edge_length

            # tree was built by simulator
            else:
                self.origin_age = None
                self.origin_edge_length = 0.0
                
                self.root_node = self.tree.seed_node
                self.root_age = self.tree.max_distance_from_root()
                root_children = [nd for nd in self.tree.nodes() if nd.label != "root"]

                ##################################################
                # Figuring out if tree died, when tree_died flag #
                # was not passed upon initialization             #
                ##################################################
                if not isinstance(self.tree_died, bool):
                    self.tree_died = True

                    # with max_age
                    if isinstance(max_age, float) and isinstance(self.root_age, float) and \
                        max_age and (max_age - self.root_age) <= self.epsilon:
                        self.tree_died = False
            
                    # no max_age
                    else:
                        try:
                            self.tree_died = self.has_tree_died(self.root_node)
                        # will get here if node doesn't have .alive
                        except:
                            pass

                # fixes root label if tree has only origin + root
                if self.tree_read_as_newick_by_dendropy and len(self.tree.leaf_nodes()) == 1:
                    self.root_node.label = self.root_node.taxon


        ###############
        # Count nodes #
        ###############
        self.n_extant_terminal_nodes = 0
        self.n_extinct_terminal_nodes = 0
        self.n_sa = 0

        # initializes
        # (i)  self.extant_obs_nodes_labels
        # (ii) self.extinct_obs_nodes_labels
        self._count_terminal_nodes() 

        # initializes
        # (i)  self.sa_obs_nodes_labels
        # (ii) self.n_sa
        self._count_sampled_ancestors()

        # initializes self.state_count_dict
        self._count_terminal_node_states() 

        # initializes
        # (i)  self.node_heights_dict
        # (ii) self.node_ages_dict
        self.populate_node_age_height_dicts()

        # prepare for dendropy's Nexus printing
        self.prepare_taxon_namespace_for_nexus_printing()


    def has_tree_died(self, origin_or_root_node) -> bool:
        """Scan all terminal nodes and returns boolean to store in self.tree_died
        
        Depends on nodes having .alive. If this member does not exist,
        tree is assumed dead.
        """
        
        # we scan all tips
        for nd in origin_or_root_node.leaf_iter():
            # if tree was read in instead of built, dendropy's Node instance
            # might not have .alive attribute
            try:
                if not nd.is_sa and nd.alive:
                    return False
            except:
                pass

        # if none of the non-SA tips is alive, or if there was something wrong
        # we assume the tree died
        return True

    # side-effect:
    # populates: self.n_sa_obs_nodes
    #            self.sa_obs_nodes_labels
    def _count_sampled_ancestors(self) -> None:
        """Count sampled ancestor nodes, store count and node labels into class members (side-effect)"""
        sa_obs_nodes_labels_list: ty.List[str] = []
        
        if isinstance(self.sa_lineage_dict, dict):
            for _, sa_list in self.sa_lineage_dict.items():
                self.n_sa += len(sa_list)
                sa_obs_nodes_labels_list.extend([sa.label for sa in sa_list])

        self.sa_obs_nodes_labels = tuple(sa_obs_nodes_labels_list)


    # side-effect:
    # populates: self.extant_obs_nodes_labels
    #            self.extinct_obs_nodes_labels
    def _count_terminal_nodes(self) -> None:
        """Count extant and extinct nodes, store counts and node labels into class members (side-effect)"""

        # nd.distance_from_root() gives distance to seed!
        extant_obs_nodes_labels_list: ty.List[str] = []
        extinct_obs_nodes_labels_list: ty.List[str] = []

        #################################################################
        # Annoying case: Tree died before stop condition was met;       #
        # We do not count SAs, and do count the only terminal node left #
        # that underwent extinction                                     #
        #################################################################
        if self.tree_died:
            # no root, single lineage from origin died
            if self.brosc_node:
                self.n_extinct_terminal_nodes = 1 # brosc node
                # if self.brosc_node:
                extinct_obs_nodes_labels_list.append(self.brosc_node.label)
            
            # with root, every lineage died
            else:
                for nd in self.tree.leaf_node_iter():
                    if not nd.is_sa and not nd.alive:
                        self.n_extinct_terminal_nodes += 1
                        extinct_obs_nodes_labels_list.append(nd.label)

        ###################################################
        # Tree did not die before stop condition was met, #
        # no "brosc" node                                 #
        ###################################################
        else:
            for nd in self.tree.leaf_node_iter():
                # an extant terminal node is a leaf whose path to the origin/root
                # has maximal length (equal to the age of the origin/root)

                # if not sampled ancestor
                if not (nd.edge_length < self.epsilon):

                    # nd.distance_from_root() returns distance to seed
                    if abs(self.seed_age - nd.distance_from_root()) <= self.epsilon:
                        if self.tree_read_as_newick_by_dendropy:
                            self.n_extant_terminal_nodes += 1
                            extant_obs_nodes_labels_list.append(nd.taxon.label)
                        # if tree was created via simulation (and .alive a member of Node)
                        elif nd.alive:
                            self.n_extant_terminal_nodes += 1
                            extant_obs_nodes_labels_list.append(nd.label)

                    # if terminal node's path is not maximal
                    else:
                        if self.tree_read_as_newick_by_dendropy:
                            self.n_extinct_terminal_nodes += 1
                            extinct_obs_nodes_labels_list.append(nd.taxon.label)
                        # if tree was created via simulation (and .alive a member of Node)
                        elif not nd.alive:
                            self.n_extinct_terminal_nodes += 1
                            extinct_obs_nodes_labels_list.append(nd.label)
                        else:
                            raise ec.AnnotatedTreeLineageMissannotation("Taxon had non-maximal age, but had \'.alive == True\'. This is not allowed. Exiting...")

        # just in case
        if self.n_extant_terminal_nodes == 0:
            self.tree_died = True

        self.extant_obs_nodes_labels = tuple(extant_obs_nodes_labels_list)
        self.extinct_obs_nodes_labels = tuple(extinct_obs_nodes_labels_list)


    # side-effect:
    # populates self.state_count_dict, {state (int): count (int), ... }
    def _count_terminal_node_states(self) -> None:
        # NOTE: we are counting ALL leaves, extinct and extant!
        for nd in self.tree.leaf_node_iter():
            # if nd.label in self.extant_obs_nodes_labels:
            try:
                self.state_count_dict[nd.state] += 1
            # if tree was read as newick, it might not have states defined
            # or the tree might have been created by hand and no states were
            # defined (e.g., Yule or birth-death processes)
            except: pass


    # TODO: add to .pyi
    def find_if_extant_or_sa_on_both_sides_complete_tr_root(self, a_node: dp.Node) -> bool:
        """
        Return True if there is at least one extant taxon or sampled ancestor
        on both sides of the root of the complete tree
        
        Note that if tree has an origin and sampled ancestors before the root,
        the reconstructed tree will be re-rooted at the dummy node that serves
        as the internal (parent) node of the sampled ancestor.
        Still, what this method checks is relative to the complete tree root node
        since a sampled ancestor does not have two sides as it is not a speciation
        event
        """
        def _recur_find_extant_or_sa(a_node: dp.Node, has_seen_root: bool, found_extant_or_sa_count: int, found_extant_or_sa: bool) -> int:
            # debugging
            # print("just started recurring at node " + a_node.label)
            # if a_node.parent_node: print("  my parent is " + a_node.parent_node.label)

            # if origin, no parent_node
            if a_node.parent_node and \
                (a_node.parent_node.label != "root" and found_extant_or_sa):
                # print("... but my parent has a child that is alive or sa, coming back!")
                return found_extant_or_sa_count, found_extant_or_sa

            has_seen_root = has_seen_root
            found_extant_or_sa = found_extant_or_sa
            
            if a_node.label == "root" or a_node.taxon == "root":
                has_seen_root = True
            
            # a_node will have two sides; will do one and then
            # another
            for ch_node in a_node.child_node_iter():
                # print("    will look at child " + ch_node.label)
                if ch_node.label == "root" or ch_node.taxon == "root":
                    has_seen_root = True

                if a_node.label != "root" and found_extant_or_sa:
                    break

                try:
                    if has_seen_root:
                        if ch_node.alive or ch_node.is_sa:
                            # print("      ... this node was alive! Breaking!")
                            found_extant_or_sa_count += 1
                            found_extant_or_sa = True

                            # found extant, stop recursion
                            # on current side of a_node if not root
                            if a_node.label != "root":
                                break

                        # stay on this side, keep digging
                        else:
                            # if root, we don't break if the sister is
                            # extant or sa, so we use False
                            if a_node.label == "root":
                                found_extant_or_sa_count, found_extant_or_sa = _recur_find_extant_or_sa(ch_node, has_seen_root, found_extant_or_sa_count, False)
                            
                            else:
                                found_extant_or_sa_count, found_extant_or_sa = _recur_find_extant_or_sa(ch_node, has_seen_root, found_extant_or_sa_count, found_extant_or_sa)

                    else:
                        # no root yet, keep digging until
                        # we find root
                        found_extant_or_sa_count, found_extant_or_sa = _recur_find_extant_or_sa(ch_node, has_seen_root, found_extant_or_sa_count, found_extant_or_sa)

                # tree must have been read as newick string
                # instead of being simulated by PJ (it doesn't
                # have .alive and .is_sa members)
                except:
                    pass

            return found_extant_or_sa_count, found_extant_or_sa
        
        n_extant_or_sa, _ = _recur_find_extant_or_sa(a_node, False, 0, False)

        # print("\nWrapping up tree, n_extant_or_sa = " + str(n_extant_or_sa))

        if n_extant_or_sa == 2:
            return True

        else:
            return False


    def recursively_find_node_age(self, a_node: dp.Node, running_age_sum: float) -> float:
        """Travel a node's path all the way to the start of the process (root or origin)
        and return path length"""
        # stop recursion
        if a_node == self.tree.seed_node:
            return running_age_sum

        # recur
        return self.recursively_find_node_age(a_node.parent_node, running_age_sum + a_node.edge_length)


    def prepare_taxon_namespace_for_nexus_printing(self) -> None:
        """Prepare taxon_namespace member of self.tree so Nexus produced by dendropy is reasonable
        
        Remove internal nodes from taxon namespace, does not affect a tree's newick representation
        """
        
        if not self.tree_read_as_newick_by_dendropy:
            for nd in self.tree.nodes():
                if not nd.is_leaf():
                    self.tree.taxon_namespace.remove_taxon(nd.taxon)


    # side-effect:
    # populates self.tree_reconstructed
    # 
    # if rejection sampling happened, is_tr_ok will populate it;
    # otherwise, population happend upon writing
    def extract_reconstructed_tree(self, require_obs_both_sides: ty.Optional[bool]=None) -> dp.Tree:
        """Make deep copy of self.tree, then prune extinct taxa from copy"""

        # if method called by someone other than Tree obj,
        # then require_obs_both_sides won't be None
        require_obs_both_sides_root = True
        if not require_obs_both_sides == None:
            require_obs_both_sides_root = require_obs_both_sides 
        # otherwise, we use the member inside Tree
        else:
            require_obs_both_sides_root = self.condition_on_obs_both_sides_root

        if self.tree_reconstructed:
            return self.tree_reconstructed

        filter_fn = lambda nd: nd.is_sa or (nd.is_leaf() and nd.alive)
        self.tree_reconstructed = copy.deepcopy(self.tree)

        # if tree went extinct, return empty new tree
        if (self.n_extant_terminal_nodes + self.n_sa) <= 1:
            return dp.Tree()
        
        self.tree_reconstructed.filter_leaf_nodes(
            filter_fn, suppress_unifurcations=True)

        # getting all root information #
        # root_node will be None if it does not exist
        smallest_distance = 0.0
        root_node: dp.Node = dp.Node()
        root_node = self.tree_reconstructed.find_node_with_label("root")
        if not root_node:
            root_node = self.tree_reconstructed.find_node_with_taxon_label("root")
        root_node_distance_from_seed = 0.0
        if root_node:
            root_node_distance_from_seed = root_node.distance_from_root()
            smallest_distance = root_node_distance_from_seed
        
        int_node_deeper_than_root: dp.Node = dp.Node()
        
        for internal_nd in self.tree_reconstructed.internal_nodes():
            internal_nd_distance = internal_nd.distance_from_root()

            if not internal_nd.label in ("root", "origin"):
                # we set int_node_deeper_than_root if we haven't
                # looked at any internal nodes yet, or if they are
                # deeper than the last one we checked
                #
                # we also only have an int_node_deeper_than_root
                # if there is (i) a root in the first place, or
                # (ii) if the root is not the seed_node
                if (smallest_distance == 0.0 or \
                   internal_nd_distance < smallest_distance) and \
                   (root_node_distance_from_seed != 0.0 or \
                    not root_node):
                    int_node_deeper_than_root = internal_nd
                    smallest_distance = internal_nd_distance


        #################
        # Special cases #
        #################

        if not require_obs_both_sides_root:
            # no speciation happened, so
            # complete tree has no root
            if not root_node:
                self.tree_reconstructed.reroot_at_node(int_node_deeper_than_root)
                int_node_deeper_than_root.edge_length = 0.0

                # if suppress_unifurcations set to False
                # origin_node = self.tree_reconstructed.seed_node
                # int_node_deeper_than_root.remove_child(origin_node)
                

            # there is a root in the complete tree
            else:
                # mrca seems to have side-effect (this shouldn't be the case...)
                # so we do things on deep copy of tree, and use labels

                # there must be an SA before the root
                # so we need to reroot above complete
                # tree's root
                if int_node_deeper_than_root.label == "None":
                    self.tree_reconstructed.reroot_at_node(int_node_deeper_than_root)

                    # getting origin and removing it
                    # if possible
                    origin_node_rec: dp.Node = dp.Node()
                    if self.origin_node:
                        origin_node_rec = self.tree_reconstructed.find_node_with_label("origin")
                        
                        if not origin_node_rec:
                            origin_node_rec = self.tree_reconstructed.find_node_with_taxon_label("origin")
                        
                        int_node_deeper_than_root.remove_child(origin_node_rec)

                    int_node_deeper_than_root.edge_length = 0.0

                # root is the deepest internal node
                else:
                    rec_tree_mrca_label = pj_get_mrca_obs_terminals(
                        self.tree_reconstructed.seed_node, [l.label for l in self.tree_reconstructed.leaf_node_iter()]
                    )

                    # re-seed if necessary
                    if rec_tree_mrca_label != "root":
                        rec_mrca_node = self.tree_reconstructed.find_node_with_label(rec_tree_mrca_label)
                        self.tree_reconstructed.reroot_at_node(rec_mrca_node)

                        # if suppress_unifurcations is set to False
                        #
                        # if there is a root AND an origin
                        # it is the origin that will dangle, so we
                        # remove it
                        # origin_node_rec: dp.Tree = dp.Tree()
                        # if self.origin_node:
                        #     origin_node_rec = self.tree_reconstructed.find_node_with_label("origin")
                        
                        #     if not origin_node_rec:
                        #         origin_node_rec = self.tree_reconstructed.find_node_with_taxon_label("origin")

                        #     rec_mrca_node.remove_child(origin_node_rec)
                    
                        # no origin! will delete dangling root
                        # else:
                        #     rec_mrca_node.remove_child(root_node)
                        
                        rec_mrca_node.edge_length = 0.0

                    # complete tree root is rec tree root
                    else:
                        root_node.edge_length = 0.0

        # require observed taxa on both sides of root
        else:
            # can't plot if need root
            if not root_node:
                return dp.Tree()

            else:
                # can't plot if taxa not on both sides of root
                if not self.find_if_extant_or_sa_on_both_sides_complete_tr_root(root_node):
                    return dp.Tree()

                if self.origin_node and int_node_deeper_than_root.taxon != None:
                    self.tree_reconstructed.reroot_at_node(int_node_deeper_than_root)

                    # if suppress_unifurcations is set to False
                    #
                    # getting origin and removing it
                    # if possible
                    # origin_node_rec: dp.Node = dp.Node()
                    # if self.origin_node:
                    #     origin_node_rec = self.tree_reconstructed.find_node_with_label("origin")
                        
                    #     if not origin_node_rec:
                    #         origin_node_rec = self.tree_reconstructed.find_node_with_taxon_label("origin")

                    # int_node_deeper_than_root.remove_child(origin_node_rec)
                    
                    int_node_deeper_than_root.edge_length = 0.0

        return self.tree_reconstructed
            

    def name_internal_nodes(self):
        pass


    # side-effect:
    # populates self.node_heights dict { node label (str): node height (float), ... }
    def populate_node_age_height_dicts(self, unit_branch_lengths=False) -> None:
        """Populate self.node_heights dict { node label: node height, ...}"""

        # side-effect:
        # populates self.node_heights_dict { node label (str): node age (float, ... }
        # or
        # self.node_heights_dict AND self.node_ages_dict { node label (str): node age (float, ... }
        def recur_node_ages_height(nd, remaining_height, tree_lvl) -> None:
            for ch_nd in nd.child_nodes():
                recur_node_ages_height(ch_nd, remaining_height - ch_nd.edge_length, tree_lvl + 1)

            if nd.label:
                if unit_branch_lengths:
                    self.node_heights_dict[nd.label] = tree_lvl
                else:
                    self.node_ages_dict[nd.label] = remaining_height
                    self.node_heights_dict[nd.label] = self.seed_age - remaining_height

            elif nd.taxon and type(nd.taxon) == str:
                if unit_branch_lengths:
                    self.node_heights_dict[nd.taxon] = tree_lvl
                else:
                    self.node_ages_dict[nd.taxon] = remaining_height
                    self.node_heights_dict[nd.taxon] = self.seed_age - remaining_height

            elif nd.taxon and nd.taxon.label:
                if unit_branch_lengths:
                    self.node_heights_dict[nd.taxon.label] = tree_lvl
                else:
                    self.node_ages_dict[nd.taxon.label] = remaining_height
                    self.node_heights_dict[nd.taxon.label] = self.seed_age - remaining_height

        tree_lvl = 0
        recur_node_ages_height(self.tree.seed_node, self.seed_age, tree_lvl)


    # side-effect:
    # populates self.node_attr_dict { node label (str): { attribute (str): val (Any)... }, ... }
    def populate_nd_attr_dict(self, attrs_of_interest_list) -> None:
        for nd in self.tree:
            if not nd.label: exit("Must name all nodes before populating nd_attr_dict. Exiting")

            for att, att_val in nd.__dict__.items():
                if att in attrs_of_interest_list:
                    self.node_attr_dict[nd.label][att] = att_val

        # self.node_attr_dict


    def __str__(self) -> str:
        return self.tree.as_string(schema="newick", suppress_annotations=False, suppress_internal_taxon_labels=True)


    def plot_node(self, axes: plt.Axes, node_attr=None, **kwargs) -> None:
        if not node_attr:
            return plot_ann_tree(self, axes)

        else:
            self.populate_nd_attr_dict([node_attr])
            return plot_ann_tree(self, axes, attr_of_interest=node_attr)


    def get_stats_dict(self) -> ty.Dict[str, ty.Union[int, float]]:
        ks = [ "Origin age", "Root age", "Total taxon count", "Extant taxon count", "Extinct taxon count", "Direct ancestor count" ]
        vs = [ self.origin_age, self.root_age, (self.n_extant_terminal_nodes + self.n_extinct_terminal_nodes), self.n_extant_terminal_nodes, self.n_extinct_terminal_nodes, self.n_sa]
        
        return dict((ks[i], str(vs[i])) for i in range(len(ks)))

    # TODO: add to .pyi
    def _get_taxon_states_dict(self, exclude_internal: bool=True) -> ty.Dict[str, int]:
        self.populate_nd_attr_dict(["state"])
        
        taxon_states_dict: ty.Dict[str, int] = dict()
        for taxon_name, attr_val_dict in self.node_attr_dict.items():
            # removing internal taxa
            if (exclude_internal and \
               self.tree.find_node_with_label(taxon_name).is_internal()):
               continue

            # removing extinct taxa if tree was simulated with
            # PJ (i.e., has .alive member)
            try:
                if not self.tree.find_node_with_label(taxon_name).alive:
                    continue
            
            except:
                pass

            taxon_states_dict[taxon_name] = attr_val_dict["state"]

        return taxon_states_dict


    def get_taxon_states_str(self, nexus: bool=False) -> str:
        taxon_states_dict = self._get_taxon_states_dict()
        states_str = ""

        for taxon_name, taxon_state in taxon_states_dict.items():
            states_str += taxon_name + "\t" + str(taxon_state) + "\n"

        if nexus:
            nexus_header = \
            "#Nexus\n\nBegin data;\nDimensions ntax=" + \
            str(self.n_extant_terminal_nodes + self.n_sa) + \
            " nchar=1;\nFormat datatype=Standard symbols=\"" + \
            "".join(str(i) for i in range(self.state_count)) + \
            "\" missing=? gap=-;\nMatrix\n"
            
            nexus_str = nexus_header + states_str + ";\nEnd;\n"
            
            return nexus_header

        return states_str


###########################
# Plotting tree functions #
###########################

color_map = {0: "deepskyblue", 1: "magenta", 2: "darkorange", 3: "gold", 4: "lawngreen"}

def get_node_name(nd):
    if nd.label: return nd.label
    elif nd.taxon and type(nd.taxon) == str: return nd.taxon
    elif nd.taxon and nd.taxon.label: return nd.taxon.label


def get_x_coord_from_nd_heights(ann_tr, use_age=False, unit_branch_lengths=False):
    """Get dictionary of node labels as keys, node y_coords (time) as values

    Args:
        ann_tr (AnnotatedTree): Annotated dendropy tree
    """

    if use_age and not unit_branch_lengths:
        return ann_tr.node_ages_dict

    return ann_tr.node_heights_dict


def get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin=False, sa_along_branches=False):
    """Get dictionary of node labels as keys, y-coords as values

    y-coords here are integers that go from 1 to the total number of
    observable nodes. Every observable node will be 1 y-unit away from
    each other.

    Args:
        ann_tr (AnnotatedTree): Annotated dendropy tree
        start_at_origin (bool): 'True' if tree starts at origin
        sa_along_branches (bool): if 'True', SA leaves are ignored so that later they are placed along branches
    """

    # max height of graph is given by number of observable nodes
    maxheight: int = 0
    if sa_along_branches:
        maxheight = ann_tr.n_extant_terminal_nodes + ann_tr.n_extinct_terminal_nodes
    else:
        maxheight = ann_tr.n_extant_terminal_nodes + ann_tr.n_extinct_terminal_nodes + ann_tr.n_sa

    # do leaves
    leaf_names = list()
    for nd in ann_tr.tree.leaf_node_iter():
        # if user wants to place SA nodes on branches,
        # we must ignore SA nodes here
        if sa_along_branches and nd.is_sa:
            continue
        
        leaf_names.append(get_node_name(nd))

    y_coords = { leaf_name: maxheight - i + 1 for i, leaf_name in enumerate(reversed(leaf_names)) }

    # do internal nodes
    def recursively_calculate_height(nd):
        children: ty.List[dp.Node] = []
        if sa_along_branches:       
            children = [ch for ch in nd.child_nodes() if not ch.is_sa]
        else:
            children = [ch for ch in nd.child_nodes()]

        if len(children) > 0:
            for ch_nd in children:
                if ch_nd not in y_coords:
                    recursively_calculate_height(ch_nd)

            # origin (and dummy nodes if user wants to place SA nodes along branches)
            if (start_at_origin and nd == ann_tr.origin_node) or (sa_along_branches and nd.is_sa_dummy_parent):
                y_coords[get_node_name(nd)] = y_coords[get_node_name(children[0])]

            # root and all others (assume bifurcation)
            else:
                y_coords[get_node_name(nd)] = (
                    y_coords[get_node_name(children[0])] + y_coords[get_node_name(children[1])]
                    ) / 2.0

    recursively_calculate_height(ann_tr.tree.seed_node)

    return y_coords

# side-effect:
# prints tree on axes (of class matplotlib.pyplot.Axes)
def plot_ann_tree(ann_tr: AnnotatedTree,
                    axes: plt.Axes,
                    use_age: bool=False,
                    start_at_origin: bool=False,
                    attr_of_interest: str=None,
                    sa_along_branches: bool=True) -> None:

    ####################################################
    # Setting up flags and checking tree content is OK # 
    ####################################################
    if ann_tr.origin_node:
        start_at_origin = True # scoped to draw_tree()
    if sa_along_branches and not ann_tr.sa_lineage_dict:
        # throw plotting exception, cannot plot sa along branches
        # without dict
        pass


    # grabbing key coordinates for drawing tree
    x_coords = get_x_coord_from_nd_heights(ann_tr, use_age=use_age)
    y_coords = get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin=start_at_origin, sa_along_branches=sa_along_branches)

    # arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # make sure we reset things
    axes.cla()

    # axes = None
    # if axes is None:
        # original
        # fig = plt.figure()
        # axes = fig.add_subplot(1, 1, 1)

        #print("oi")

        # fig = plt.figure(figsize=(11,4.5))
        # axes = fig.add_axes([0.1, 0.15, 0.8, 0.7]) # left, bottom, width, height

    # elif not isinstance(axes, plt.matplotlib.axes.Axes):
    #    raise ValueError(f"Invalid argument for axes: {axes}")

    # side-effect:
    # populates: horizontal_linecollections
    #            vertical_linecollections
    def _draw_clade_lines(
                        x_end=0.0,
                        use_linecollection=False,
                        orientation="horizontal",
                        y_here=0,
                        x_start=0,
                        x_here=0,
                        y_bot=0,
                        y_top=0,
                        color="black",
                        lw=".1"
                        ) -> None:

        """Create a line with or without a line collection object.
        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if not use_linecollection and orientation == "horizontal":
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif use_linecollection and orientation == "horizontal":
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

        elif not use_linecollection and orientation == "vertical":
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif use_linecollection and orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )


    def _draw_clade(nd, x_start, color, lw, use_age, start_at_origin) -> None:
        """Recursively draw a tree, down from the given node"""
        nd_name = get_node_name(nd)

        x_here = x_coords[nd_name]
        y_here = y_coords[nd_name]

        from_x_here_to_this_x: float=0.0 # default (present moment)
        # draws from origin/root_node to root_node x
        if use_age and nd in (ann_tr.origin_node, ann_tr.root_node):
            if isinstance(ann_tr.root_age, float):
                from_x_here_to_this_x = ann_tr.root_age

        color = "black"
        if attr_of_interest:
            attr_idx = ann_tr.node_attr_dict[nd_name][attr_of_interest]
            color = color_map[attr_idx]

        #################################
        # Draw horizontal line (branch) #
        #################################
        _draw_clade_lines(
            x_end=from_x_here_to_this_x,
            use_linecollection=True,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw
        )

        #########################
        # Add node/taxon labels #
        #########################
        if nd.is_leaf():
            axes.text(
                x_here,
                y_here,
                f" {nd_name}",
                verticalalignment="center",
            )

        #################################
        # Draw sampled ancestors if any #
        #################################
        if sa_along_branches and (not nd.is_sa_dummy_parent and nd.is_sa_lineage):
                sas: ty.List[pjsa.SampledAncestor]=[]
                if isinstance(ann_tr.sa_lineage_dict, dict):
                    sas = ann_tr.sa_lineage_dict[nd_name]
                
                for sa in sas:
                    sa_x: float = sa.global_time
                    # sa_x: float = x_here - sa.time_to_lineage_node # works fine
                    
                    if use_age:
                        if start_at_origin and isinstance(ann_tr.origin_age, float):
                            sa_x = ann_tr.origin_age - sa_x
                        elif isinstance(ann_tr.root_age, float):
                            sa_x = ann_tr.root_age - sa_x
                    
                    sa_y = y_here
                    
                    _draw_sa_along_branch(sa.label, sa.global_time, sa_y, axes)

        # recur if not leaf
        if len(nd.child_nodes()) > 0:           
            children = nd.child_nodes()   

            ################################################
            # Draw a vertical line connecting two children #
            ################################################
            # not origin
            if not start_at_origin or (ann_tr.origin_node and nd != ann_tr.origin_node):

                # must be either (1) placing SA along branch, and then NOT a dummy, or
                # or             (2) regular plotting, AND have two children
                if (sa_along_branches and not nd.is_sa_dummy_parent) or \
                    (len(children) == 2 and not sa_along_branches):
                    y_top = y_coords[get_node_name(children[0])]

                    # print("y_top for node " + children[0].label + " = " + str(y_top))

                    y_bot = y_coords[get_node_name(children[1])]

                    _draw_clade_lines(
                        use_linecollection=True,
                        orientation="vertical",
                        x_here=x_here,
                        y_bot=y_bot,
                        y_top=y_top,
                        color=color,
                        lw=lw,
                    )

            ####################
            # Draw descendents #
            ####################
            for child_nd in children:
                # if user wants to place SA nodes along branches, we must ignore them here
                # if sa_along_branches and (child_nd.is_sa or child_nd.is_sa_dummy_parent):
                if sa_along_branches and child_nd.is_sa:
                    continue

                # draws horizontal line (branch) inside 
                _draw_clade(child_nd, x_here, color, lw, False, start_at_origin)


    def _draw_time_slices(ann_tr: AnnotatedTree, axes: plt.Axes, use_age: bool) -> None:
        """Draw vertical lines representing boundaries of time slices"""
        xs: ty.List[float] = []
        if use_age and ann_tr.slice_age_ends:
            xs = ann_tr.slice_age_ends[1:] # ignore present
        elif ann_tr.slice_t_ends:
            xs = [t_end for t_end in ann_tr.slice_t_ends[:-1] if isinstance(t_end, float)] # so mypy won't complain

        for x in xs:
            axes.axvline(x=x, color="deeppink")

    def _draw_sa_along_branch(sa_name: str, x: float, y: float, axes: plt.Axes) -> None:
        """Draw SA nodes along internal branches"""
        # zorder places dot on top of axes
        # k = black color, o = means marker
        axes.plot(x, y, marker="o", color="k", markersize=4, zorder=10)               
        axes.text(x, y - 0.05, f" {sa_name}", rotation=45, horizontalalignment="center", fontsize=10)


    ########
    # call #
    ########
    _draw_clade(ann_tr.tree.seed_node, 0, "k", plt.rcParams["lines.linewidth"], use_age, start_at_origin)
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
    # axes.set_facecolor("gray") later maybe change colors
    # plt.show()
    # return plt.gcf()%

def pj_get_mrca_obs_terminals(a_node: dp.Node, nd_label_list: ty.List[str]) -> dp.Node:
    # side-effect recursion
    def recur_node(a_node, nd_label_list: ty.List[str], mrca_node_label: str):
        """Populate visited_obs_terminals (side-effect)"""        

        # not done: if tip, get label
        if a_node.is_leaf() and (a_node.is_sa or a_node.alive):
            return [a_node.label], ""

        # not done: if internal node, we recur
        else:
            visited_obs_terminals: ty.List[str] = []

            for ch_node in a_node.child_node_iter():
                # recur
                vot, mrca_node_label = recur_node(ch_node, nd_label_list, mrca_node_label)
                visited_obs_terminals += vot

            # we are actually done
            if set(visited_obs_terminals) == set(nd_label_list) and not mrca_node_label:
                mrca_node_label  = a_node.label

        return visited_obs_terminals, mrca_node_label

    mrca_node_label: str = ""
    
    # mrca node
    _, mrca_node_label = recur_node(a_node, nd_label_list, mrca_node_label)
    
    return mrca_node_label
 
    

##############################################################################

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
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 src/phylojunction/data/tree.py

    ##############################################
    # Below we have an environment for debugging #
    ##############################################
    origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
    origin_node.state = 0
    origin_node.alive = False
    origin_node.is_sa = False
    origin_node.is_sa_dummy_parent = False
    origin_node.is_sa_lineage = False

    dummy_node = dp.Node(taxon=dp.Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
    dummy_node.state = 0
    dummy_node.alive = False
    dummy_node.is_sa = False
    dummy_node.is_sa_dummy_parent = True
    dummy_node.is_sa_lineage = False
    
    origin_node.add_child(dummy_node)
    
    # right child of dummy_node
    sa_node = dp.Node(taxon=dp.Taxon(label="sa1"), label="sa1", edge_length=0.0)
    sa_node.state = 0
    sa_node.alive = False
    sa_node.is_sa = True
    sa_node.is_sa_dummy_parent = False
    sa_node.is_sa_lineage = False

    # left child of dummy node
    root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.5)
    root_node.state = 0
    root_node.alive = False
    root_node.is_sa = False
    root_node.is_sa_dummy_parent = False
    root_node.is_sa_lineage = True
    
    dummy_node.add_child(sa_node)
    dummy_node.add_child(root_node)

    # left child of root node
    extant_sp1 = dp.Node(taxon=dp.Taxon(label="sp1"), label="sp1", edge_length=0.25)
    extant_sp1.state = 0
    extant_sp1.alive = False
    extant_sp1.is_sa = False
    extant_sp1.is_sa_dummy_parent = False
    extant_sp1.is_sa_lineage = False

    # right child of root node
    extant_sp2 = dp.Node(taxon=dp.Taxon(label="sp2"), label="sp2", edge_length=0.5)
    extant_sp2.state = 1
    extant_sp2.alive = True
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

    total_state_count = 2

    sa_global_time = 1.0
    time_to_sa_lineage_node = 0.5
    sa = pjsa.SampledAncestor("sa1", "root", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
    sa_lineage_dict = { "root": [sa] }
    
    max_age = 2.0

    ann_tr_sa_with_root_survives_max_age = AnnotatedTree(
        tr_sa_with_root_survives,
        total_state_count,
        start_at_origin=True,
        max_age=max_age,
        sa_lineage_dict=sa_lineage_dict,
        epsilon=1e-12)

    print(ann_tr_sa_with_root_survives_max_age.tree.as_string(schema="newick"))
    print(ann_tr_sa_with_root_survives_max_age.get_taxon_states_nexus_str())

    ######################
    # Preparing plotting #
    ######################
    # fig = Figure(figsize=(11,4.5))
    fig = matplotlib.pyplot.figure() # note that pjgui uses matplotlib.figure.Figure (which is part of Matplotlib's OOP class library)
                                     # here, we instead use pyplot's figure, which is the Matlab-like state-machine API
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
                     sa_along_branches=True
                     # sa_along_branches=False
                     )

    plt.show()