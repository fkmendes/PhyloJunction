import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/data/)
import typing as ty
import dendropy as dp # type: ignore
# from dendropy import Node, Tree, Taxon 

# plotting tree
import matplotlib.pyplot as plt # type: ignore
import matplotlib.collections as mpcollections # type: ignore
from matplotlib.figure import Figure # type: ignore

# pj imports
import utility.exception_classes as ec
import utility.helper_functions as pjh

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class AnnotatedTree(dp.Tree):

    with_origin: bool
    origin_node: ty.Optional[dp.Node] # can be None
    origin_age: ty.Optional[float] # can be None
    root_node: ty.Optional[dp.Node] # can be None
    root_age: ty.Optional[float] # can be None
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
    slice_age_ends: ty.Optional[ty.List[float]] # can be None
    n_extant_obs_nodes: int
    n_extinct_obs_nodes: int
    n_sa_obs_nodes: int
    extant_obs_nodes_labels: ty.Tuple[str, ...]
    extinct_obs_nodes_labels: ty.Tuple[str, ...]

    def __init__(self,
                a_tree: dp.Tree,
                total_state_count: int,
                start_at_origin: bool=False,
                max_age: ty.Optional[float]=None,
                slice_t_ends: ty.List[ty.Optional[float]]=[],
                slice_age_ends: ty.Optional[ty.List[float]]=None,
                epsilon: float=1e-12):

        self.with_origin = start_at_origin
        self.origin_node = None
        self.root_node = None
        self.tree_read_as_newick_by_dendropy = False
        self.tree = a_tree
        self.state_count = total_state_count
        self.seed_age = self.tree.max_distance_from_root()
        self.epsilon = epsilon
        self.tree_died = False
        self.state_count_dict = dict((int(s), 0) for s in range(self.state_count))
        self.node_heights_dict: ty.Dict[str, float] = dict()
        self.node_ages_dict: ty.Dict[str, float] = dict()
        self.node_attr_dict = pjh.autovivify(2) # populated on demand by populate_nd_attr_dict
        self.slice_t_ends = slice_t_ends
        self.slice_age_ends = slice_age_ends

        try:
            self.tree.seed_node.alive
        except:
            self.tree_read_as_newick_by_dendropy = True

        if self.with_origin:
            # making sure this is indeed origin...
            if len(self.tree.seed_node.child_nodes()) > 1:
                raise ec.AnnotatedTreeMisspec("Tree was specified as starting from origin, but it started from root. Exiting...")

            self.origin_node = self.tree.seed_node
            self.origin_age = self.seed_age
            self.root_node = self.tree.seed_node.child_nodes()[0] # origin always have root child
            self.root_edge_length = self.root_node.edge_length
            self.root_age = self.origin_age - self.root_edge_length # returns 0.0 if tree dies (TODO: revisit this...)

            # fixes root label if tree has only origin + root
            if self.tree_read_as_newick_by_dendropy and len(self.tree.leaf_nodes()) == 1:
                self.root_node.label = str(self.root_node.taxon).replace("\'", "")

            # the tree might have died later, but if this is satisfied, it is certainly dead already
            if not self.root_node.child_nodes() and max_age and (abs(self.root_edge_length - max_age) > self.epsilon):
                self.tree_died = True

        else:
            if len(self.tree.seed_node.child_nodes()) == 1:
                raise ec.AnnotatedTreeMisspec("Tree was specified as not starting from origin, but it starts from the origin. Exiting...")
            # if tree is read as a Newick string, the user forgot to specify "start_at_origin=True".
            # and there is a root edge, we need to add an origin node and update everything accordingly
            if len(self.tree.seed_node.child_nodes()) == 2 and self.tree.seed_node.edge_length:
                self.root_node = self.tree.seed_node
                self.root_edge_length = self.root_node.edge_length
                self.origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
                self.origin_node.alive = False
                self.origin_node.add_child(self.root_node)
                self.tree.seed_node = self.origin_node
                self.seed_age = self.tree.max_distance_from_root()
                self.origin_age = self.seed_age
                self.root_age = self.origin_age - self.root_edge_length

            else:
                self.origin_age = None
                self.root_node = self.tree.seed_node
                self.root_age = self.tree.max_distance_from_root()
                self.root_edge_length = 0.0

                # fixes root label if tree has only origin + root
                if self.tree_read_as_newick_by_dendropy and len(self.tree.leaf_nodes()) == 1:
                    self.root_node.label = self.root_node.taxon

        ###############
        # Count nodes #
        ###############

        # side-effect initializes:
        #   self.n_extant_obs_nodes
        #   self.n_extinct_obs_nodes
        #   self.extant_obs_nodes_labels
        #   self.extinct_obs_nodes_labels
        self.n_extant_obs_nodes = 0
        self.n_extinct_obs_nodes = 0
        self.n_sa_obs_nodes = 0
        self.count_observable_nodes()
        # TODO: self.count_sampled_ancestors()

        self.count_observable_node_states() # initializes self.state_count_dict

        self.populate_node_age_height_dicts() # initializes self.node_heights_dict and self.node_ages_dict


    def count_sampled_ancestors(self):
        pass


    # side-effect:
    # populates: self.extant_obs_nodes_labels
    #            self.self.extinct_obs_nodes_labels
    def count_observable_nodes(self) -> None:
        # nd.distance_from_root() gives distance to seed!
        extant_obs_nodes_labels_list: ty.List[str] = []
        extinct_obs_nodes_labels_list: ty.List[str] = []

        # tree died before first speciation event
        if self.tree_died:
            self.n_extinct_obs_nodes = 1

            if self.root_node:
                extinct_obs_nodes_labels_list.append(self.root_node.label)

        # tree did not die before first speciation event
        else:
            for nd in self.tree:
                # root (un-speciated) survived all the way to max age, did not go extinct
                if nd == self.root_node and not self.root_node.child_nodes():
                    self.n_extant_obs_nodes = 1
                    extant_obs_nodes_labels_list.append(nd.label)
                    break

                elif nd.is_leaf():
                    # nd.distance_from_root() returns distance to seed
                    # checks that distance from seed to node is maximal
                    if (abs(self.seed_age - nd.distance_from_root()) <= self.epsilon):
                        if self.tree_read_as_newick_by_dendropy:
                            self.n_extant_obs_nodes += 1
                            extant_obs_nodes_labels_list.append(nd.taxon.label)
                        # if tree was created via simulation (and .alive a member of Node)
                        elif nd.alive:
                            self.n_extant_obs_nodes += 1
                            # nd.taxon = nd.label
                            extant_obs_nodes_labels_list.append(nd.label)

                    # if distance is not maximal
                    else:
                        if self.tree_read_as_newick_by_dendropy:
                            self.n_extinct_obs_nodes += 1
                            extinct_obs_nodes_labels_list.append(nd.taxon.label)
                        # if tree was created via simulation (and .alive a member of Node)
                        elif not nd.alive:
                            self.n_extinct_obs_nodes += 1
                            # nd.taxon = nd.label
                            extinct_obs_nodes_labels_list.append(nd.label)

                elif not nd.is_leaf() and not self.tree_read_as_newick_by_dendropy:
                    # remove internal node labels from taxon namespace, for Nexus printing
                    try:
                        self.tree.taxon_namespace.remove_taxon_label(nd.label, is_case_sensitive=True)
                    except:
                        print(self.tree.taxon_namespace)

        # used when stop condition is max age
        # if speciation events happened, but tree eventually died before max age
        if self.n_extant_obs_nodes == 0:
            self.tree_died = True

        self.extant_obs_nodes_labels = tuple(extant_obs_nodes_labels_list)
        self.extinct_obs_nodes_labels = tuple(extinct_obs_nodes_labels_list)


    # side-effect:
    # populates self.state_count_dict, {state (int): count (int), ... }
    def count_observable_node_states(self) -> None:
        for nd in self.tree:
            if nd.label in self.extant_obs_nodes_labels:
                try:
                    self.state_count_dict[nd.state] += 1
                # if tree was read as newick, it might not have states defined
                # or the tree might have been created by hand and no states were
                # defined (e.g., Yule or birth-death processes)
                except: pass


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


    def get_gcf(self, axes, node_attr=None, **kwargs) -> None:
        if not node_attr:
            return get_gcf_ann_tree(self, axes)

        else:
            self.populate_nd_attr_dict([node_attr])
            return get_gcf_ann_tree(self, axes, attr_of_interest=node_attr)


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


def get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin=False):
    """Get dictionary of node labels as keys, y-coords as values

    y-coords here are integers that go from 1 to the total number of
    observable nodes. Every observable node will be 1 y-unit away from
    each other.

    Args:
        ann_tr (AnnotatedTree): Annotated dendropy tree
    """

    # max height of graph is given by number of observable nodes
    maxheight = len(ann_tr.extant_obs_nodes_labels) + len(ann_tr.extinct_obs_nodes_labels)

    # do leaves
    leaf_names = list()
    for nd in ann_tr.tree.leaf_node_iter():
        leaf_names.append(get_node_name(nd))

    y_coords = { leaf_name: maxheight - i + 1 for i, leaf_name in enumerate(reversed(leaf_names)) }

    # do internal nodes
    def recursively_calculate_height(nd):
        children = nd.child_nodes()

        if len(children) > 0:
            for ch_nd in children:
                if ch_nd not in y_coords:
                    recursively_calculate_height(ch_nd)

            # origin
            if start_at_origin and nd == ann_tr.origin_node:
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
def get_gcf_ann_tree(ann_tr: AnnotatedTree, axes: plt.Axes, use_age: bool=False, start_at_origin: bool=False, attr_of_interest: str=None) -> None:
    if ann_tr.origin_node:
        start_at_origin = True # scoped to draw_tree()

    # grabbing key coordinates for drawing tree
    x_coords = get_x_coord_from_nd_heights(ann_tr, use_age=use_age)
    y_coords = get_y_coord_from_n_obs_nodes(ann_tr, start_at_origin=start_at_origin)

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

        from_x_here_to_this_x: ty.Optional[float]=0.0 # default
        if use_age and nd in (ann_tr.origin_node, ann_tr.root_node):
            from_x_here_to_this_x = ann_tr.root_age

        color = "black"
        if attr_of_interest:
            attr_idx = ann_tr.node_attr_dict[nd_name][attr_of_interest]
            color = color_map[attr_idx]

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

        # Add node/taxon labels
        if nd.is_leaf():
            axes.text(
                x_here,
                y_here,
                f" {nd_name}",
                verticalalignment="center"
            )

        # recur
        if len(nd.child_nodes()) > 0:
            children = nd.child_nodes()

            # not origin
            if not start_at_origin or (ann_tr.origin_node and nd != ann_tr.origin_node):
                # Draw a vertical line connecting all children
                y_top = y_coords[get_node_name(children[0])]
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

            # Draw descendents
            for child_nd in children:
                _draw_clade(child_nd, x_here, color, lw, False, start_at_origin)


    def _draw_time_slices(ann_tr: AnnotatedTree, axes: plt.Axes, use_age: bool) -> None:
        """Draw vertical lines representing boundaries of time slices"""
        xs: ty.List[float]
        if use_age and ann_tr.slice_age_ends:
            xs = ann_tr.slice_age_ends[1:] # ignore present
        elif ann_tr.slice_t_ends:
            xs = [t_end for t_end in ann_tr.slice_t_ends[:-1] if isinstance(t_end, float)] # so mypy won't complain

        for x in xs:
            axes.axvline(x=x, color="deeppink")

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

##############################################################################

if __name__ == "__main__":
    # can be called from data/
    # $ python3 tree.py
    # 
    # can also be called from phylojunction/
    # $ python3 data/tree.py
    # or
    # $ python3 -m data.tree
    #
    # can also be called from VS Code, if open folder is phylojuction/

    total_state_count = 1
    rootedge_tr_str = "(((nd6:3.0,nd7:1.0)nd3:1.0,(nd4:1.0,nd5:3.0)nd2:1.0)root:2.0)origin:0.0;"
    tr_origin = dp.Tree.get(data=rootedge_tr_str, schema="newick")
    tree_origin = AnnotatedTree(tr_origin, total_state_count, start_at_origin=True, epsilon=1e-12)
    print(tree_origin)