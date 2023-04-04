"""dn_discrete_sse.py: Class for distribution of discrete state-dependent speciation and extinction process"""

import typing as ty
import time
import numpy
# import statistics
import random
# from random import seed, choice
import dendropy as dp # type: ignore
# from dendropy import Tree, Node, Taxon, TaxonNamespace
# from matplotlib.pyplot import draw

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.helper_functions as pjh
import phylojunction.utility.exception_classes as ec
import phylojunction.distribution.dn_parametric as dnpar
from phylojunction.data.tree import AnnotatedTree
from phylojunction.data.sampled_ancestor import SampledAncestor # type: ignore
from phylojunction.data.attribute_transition import AttributeTransition # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class DnSSE(pgm.DistributionPGM):
    """Distribution class for State-dependent Speciation and Extinction (birth-death) process

    The sampling method in this class is a "rising tide" simulator,
    where all living lineages grow together.

    At all times we know: (i) how many lineages are alive and dead, and (ii) all the character states represented by living lineages

    This is not a recursive simulator (in which each lineage would take care of
    growing only itself, and recur upon birth events).

    Attributes:
        event_handler (MacroevolEventHandler): Object holding all rates we need to simulate.
        n (int): Number of simulations one wants.
        n_replicates (int): Number of tree replicates per simulation.
        stop (str): If "age", stops when origin or root age is equal to 'stop_condition_value'. If "size", stops when tree has 'stop_condition_value' observable nodes.
        stop_val (float): Either maximum tree age, or maximum count of observable nodes.
        start_at_origin (bool): Simulation starts at origin (i.e., seed == origin).
        start_states (int): List of integer representing the starting states of all n_sim simulations.
        state_count (int): Number of discrete states (obtained from 'events').
        n_time_slices (int): Number of time slices (obtained from 'events').
        slice_t_ends (int): List of end times from time slices (obtained from 'events').
        condition_on_speciation (bool): If true, rejects trees that go extinct before reaching stop_condition_value. Defaults to False.
        condition_on_survival (bool): If true, rejects trees that go extinct before reaching stop_condition_value. Defaults to True.
        condition_on_obs_both_sides_root (bool): If true, rejects trees that do not have observed taxa on both sides of complete tree's root. Defaults to False.
        min_rec_taxa (int): Required minimum number of observed taxa in reconstructed tree (at least or equal to). Defaults to zero.
        max_rec_taxa (int): Required maximum number of observed taxa in reconstructed tree (at most and equal to). Defaults to 1e12.
        abort_at_obs (int): If growing tree reaches this number of living taxa, we abort the sample (invalidate the tree). Defaults to 1e12.
        epsilon (float, optional): Any values difference below this value is set to 0.0. Defaults to 1e-12.
        runtime_limit (int, optional): Maximum number of minutes one is willing to wait until all 'n_sim' simulations are done. Defaults to 5.
        debug (bool, optional): Whether to print debugging messages

    Methods

    """

    DN_NAME = "DnSSE"

    n_sim: int
    n_repl: int
    with_origin: bool
    root_is_born: bool
    # tree_died: bool
    stop: str
    stop_val: ty.List[float]
    condition_on_speciation: bool
    condition_on_survival: bool
    condition_on_obs_both_sides_root: bool
    min_rec_taxa: int
    max_rec_taxa: int
    abort_at_obs: int
    sse_stash: sseobj.SSEStash
    events: sseobj.MacroevolEventHandler
    prob_handler: sseobj.DiscreteStateDependentProbabilityHandler
    start_states: ty.List[int]
    state_count: int
    n_time_slices: int
    slice_t_ends: ty.List[ty.Optional[float]]
    seed_age: ty.Optional[float]
    seeds: ty.Optional[ty.List[int]]
    epsilon: float    
    runtime_limit: int
    debug: bool

    # TODO: later make event_handler mandatory and update typing everywhere, as well as fix all tests
    def __init__(self,
                 sse_stash: sseobj.SSEStash,
                 stop_value: ty.List[float]=[],
                 n: int=1,
                 n_replicates: int=1,
                 stop: str="",
                 origin: bool=True,
                 start_states_list: ty.List[int]=[],
                 condition_on_speciation: bool=False,
                 condition_on_survival: bool=True,
                 condition_on_obs_both_sides_root: bool=False,
                 min_rec_taxa: int=0,
                 max_rec_taxa: int=int(1e12),
                 abort_at_obs: int=int(1e12),
                 seeds_list: ty.Optional[ty.List[int]]=None,
                 epsilon: float=1e-12,
                 runtime_limit: int=15,
                 debug: bool=False) -> None:

        # simulation parameters
        self.n_sim = int(n)
        self.n_repl = int(n_replicates) # number of replicate trees (in plate) for a given simulation
        self.with_origin = origin
        self.root_is_born = False
        # self.tree_died = False
        self.stop = stop
        self.stop_val = stop_value
        self.condition_on_speciation = condition_on_speciation
        self.condition_on_survival = condition_on_survival
        self.condition_on_obs_both_sides_root = condition_on_obs_both_sides_root

        # rejection sampling
        self.min_rec_taxa = min_rec_taxa
        self.max_rec_taxa = max_rec_taxa
        self.abort_at_obs = abort_at_obs

        # model parameters #
        self.start_states = start_states_list
        self.sse_stash = sse_stash
        # all rates, number of states and of slices
        self.events = sse_stash.get_meh()
        self.state_count = self.events.state_count
        self.n_time_slices = self.events.n_time_slices
        # sampling probabilities, if provided
        self.prob_handler = sse_stash.get_prob_handler()

        # for time-heterogeneous processes
        self.slice_t_ends = self.events.slice_t_ends
        self.seed_age = self.events.seed_age # used just for verifying inputs

        # other specs
        self.seeds = seeds_list # for randomization (this functionality is not working)
        self.epsilon = epsilon
        self.runtime_limit = runtime_limit
        self.debug = debug

        # checking number of provided values (vectorizing if necessary)
        self._check_sample_size()

        ##################
        # Input checking #
        ##################
        if self.n_sim < 1:
            raise ec.DnInitMisspec(self.DN_NAME, "Please specify a number of simulations >= 1. Exiting...")
        
        if not self.start_states:
            raise ec.DnInitMisspec(self.DN_NAME, "You must provide a list of " + str(self.n_sim) + " starting states. Exiting...")
        
        if self.stop != "age" and self.stop != "size":
            raise ec.DnInitMisspec(self.DN_NAME, "Stop condition must be \"age\" or \"size\". Exiting...")
        
        if self.n_time_slices > 1 and not self.stop == "age":
            raise ec.DnInitMisspec(self.DN_NAME, "If you specified more than a single time slice, you need to provide a maximum age as stopping condition to anchor the absolute slice ages. Exiting...")
        
        if self.n_time_slices > 1 and not self.seed_age:
            raise ec.DnInitMisspec(self.DN_NAME, "When providing time-slice age ends, you must specify a seed (origin or root) age to anchor the tree and ensure the tree spans those slices. Exiting...")
        
        if self.stop == "age" and self.seed_age:
            # TODO: make sure seed age in FIGManager is vectorized and that in script, the same variable is used for both DiscreteStateDependentParameterManager and stop condition value in the spec of this distribution
            for idx, a_stop_val in enumerate(self.stop_val):
                # TODO: will need to do self.seed_age[idx] later when vectorization is added
                if a_stop_val != self.seed_age:
                    raise ec.DnInitMisspec(self.DN_NAME, "The max. stopping age(s) for the simulation must match the specified seed (origin or root) age(s). Exiting...")
        
        # the checks below are also carried out at the grammar level, but we do it again
        # in case distribution is used without a script
        if self.stop == "size":
            for idx, a_stop_val in enumerate(self.stop_val):
                # must be a number
                if not isinstance(a_stop_val, (int, float)):
                    raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (number of terminal nodes) must be a number. Exiting...")
                # it is a number...
                else:
                    # can be a int-convertible float
                    if isinstance(a_stop_val, float):
                        if a_stop_val.is_integer():
                            self.stop_val[idx] = int(a_stop_val)
                        # float must be convertible to integer
                        else: 
                            raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (number of terminal nodes) could not be converted to integer. Exiting...")
                    # ideally should be int
                    elif isinstance(a_stop_val, int):
                        self.stop_val[idx] = a_stop_val
                    # once self.stop_val[idx] is initialized, it must be >= 0
                    if self.stop_val[idx] < 0:
                        raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (number of terminal nodes) cannot be negative. Exiting...")
                    # if it is 1, then the process cannot start at a root node, which by definition has two children already
                    if self.stop_val[idx] == 1 and not self.with_origin:
                        # TODO: add this check to the grammar level
                        raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (number of terminal nodes) cannot be 1 for a process starting at the root." + \
                            " This is only allowed if the process starts at the origin. Exiting...")
        
        if self.stop == "age":
            for idx, a_stop_val in enumerate(self.stop_val):
                # must be an integer
                if not isinstance(a_stop_val, float):
                    raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (tree height) must be an integer. Exiting...")
                else:
                    self.stop_val[idx] = a_stop_val
                    # must be >= 0
                    if self.stop_val[idx] < 0.0:
                        raise ec.DnInitMisspec(self.DN_NAME, "Stop condition value (tree height) cannot be negative. Exiting...")
                

    ##########################
    # Initialization methods #
    ##########################

    def _check_sample_size(self) -> None:
        # changing atomic_rate_params_matrix here
        # should affect state_dep_rate_manager
        #
        # 1D: time slices, 2D: list of atomic rate params
        atomic_rate_params_matrix = \
            self.events.state_dep_rate_manager.matrix_state_dep_params

        state_dep_probs_matrix = \
            self.prob_handler.state_dep_prob_manager.matrix_state_dep_params

        ############################
        # Checking rate parameters #
        ############################
        # k-th time slice
        for k, list_atomic_rate_params in enumerate(atomic_rate_params_matrix):
            for arp in list_atomic_rate_params:
                n_val = len(arp.value)

                # we multiply values if one is provided, but > 1 sims
                if n_val == 1 and (self.n_sim > 1):
                    arp.value = [arp.value[0] for i in range(self.n_sim)]
                
                # do not know how to multiply, error!
                elif n_val > 1 and n_val != self.n_sim:
                    raise ec.DimensionalityError(self.DN_NAME)

        ###################################
        # Checking sampling probabilities #
        ###################################
        for k, list_state_dep_probs in enumerate(state_dep_probs_matrix):
            for state_dep_prob in list_state_dep_probs:
                n_val = len(state_dep_prob.value)

                # we multiply values if one is provided, but > 1 sims
                if n_val == 1 and (self.n_sim > 1):
                    state_dep_prob.value = [state_dep_prob.value[0] for i in range(self.n_sim)]
                
                # do not know how to multiply, error!
                elif n_val > 1 and n_val != self.n_sim:
                    raise ec.DimensionalityError(self.DN_NAME)

        ############################
        # Checking starting states #
        ############################
        if len(self.start_states) > self.n_sim or (len(self.start_states) > 1 and len(self.start_states) < self.n_sim):
            raise ec.DimensionalityError(self.DN_NAME)

        # multiplying starting states if only one was passed and the number of simulations is > 1
        elif len(self.start_states) < self.n_sim and len(self.start_states) == 1:
            self.start_states = [int(self.start_states[0]) for i in range(self.n_sim)]

        ########################
        # Checking stop values #
        ########################
        if len(self.stop_val) > self.n_sim or (len(self.stop_val) > 1 and len(self.stop_val) < self.n_sim):
            raise ec.DimensionalityError(self.DN_NAME)

        # multiplying stop values if only one was passed and the number of simulations is > 1
        elif len(self.stop_val) < self.n_sim and len(self.stop_val) == 1:
            self.stop_val = [float(self.stop_val[0]) for i in range(self.n_sim)]

        return None # dummy return


    ######################
    # Simulation methods #
    ######################
    
    def get_next_event_time(self, total_rate: float, a_seed: ty.Optional[int]=None) -> float:
        """Draw next exponentially distributed event time

        Args:
            current_node_target_count (int): Number of nodes that can experience an event.
            total_rate (float): Sum of all rates currently represented by nodes that can experience an event.

        Returns:
            float: Time to next event.
        """

        if a_seed:
            numpy.random.seed(a_seed)

        next_time = float(dnpar.DnExponential.draw_exp(1, total_rate, rate_parameterization=True))

        return next_time


    def execute_birth(self,
                    tr_namespace: dp.TaxonNamespace,
                    chosen_node: dp.Node,
                    state_representation_dict: ty.Dict[int, ty.Set[str]],
                    sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                    state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                    untargetable_node_set: ty.Set[str],
                    cumulative_node_count: int,
                    macroevol_atomic_param: sseobj.DiscreteStateDependentRate,
                    event_t: float,
                    debug=False) -> ty.Tuple[dp.Node, int]:
        """Execute lineage birth (side-effect and return)

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo speciation.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            sa_lineage_dict (dict): Dictionary that keeps track of nodes that have sampled ancestor children..
            state_transition_dict (dict): Dictionary that keeps track of nodes subtending which state transitions happen (used for plotting)
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            macroevol_atomic_param (AtomicRateParameter): Atomic rate parameter containing departing/arriving state.
            event_t (float): Time of birth event taking place.
            debug (bool): If 'true', prints debugging messages. Defaults to False.

        Returns:
            (dendropy.Node, int): Tuple with last node to under go event and total (cumulative) node count.
        """

        left_arriving_state, right_arriving_state = \
            macroevol_atomic_param.arriving_states

        if debug:
            print("> SPECIATION of node " + chosen_node.label + " in state " + str(chosen_node.state) + \
                " into daughters with states " + str(left_arriving_state) + " and " + str(right_arriving_state))

        ############################################################
        # Special case: first speciation event (root must created) #
        ############################################################
        if chosen_node.label in ("origin", "brosc"):
            root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root")
            root_node.state = chosen_node.state # root takes origin or brosc state

            root_node.annotations.add_bound_attribute("state")
            root_node.alive = False # we will pick the root as the first event
            root_node.is_sa = False
            root_node.is_sa_dummy_parent = False
            root_node.is_sa_lineage = chosen_node.is_sa_lineage
            tr_namespace.add_taxon(root_node.taxon)
            self.root_is_born = True
            # state_representation_dict[root_node.state].add(root_node.label)

            # origin/brosc cannot be selected anymore
            chosen_node.alive = False # origin/brosc node is no longer alive
            untargetable_node_set.add(chosen_node.label) # cannot be targeted
            
            state_representation_dict[chosen_node.state].remove(chosen_node.label) # and does not represent state

            # updating origin
            if chosen_node.label == "origin":
                # at this point, the origin was chosen to undergo a birth event, but this node
                # will never be extended (the origin always has an origin_edge_length = 0.0);
                # the root does NOT exist yet here; for us to account for the evolution (branch length)
                # that has happened between the origin/brosc to this moment (when the root is born), we
                # must add this branch length (event_t) now when creating the root (unlike other nodes,
                # which are created with edge_length = 0.0)
                root_node.edge_length = event_t # if edge_length = event_t, fbd unit test passes
                chosen_node.add_child(root_node) # finally link root to origin

            # root replaces brosc (dummy_node might or not exist)
            # [ori] ---> (dummy_node) ---> [brosc] 
            # or
            # [ori] ---> [brosc]
            # becomes
            # [ori] ---> (dummy_node) ---> [root] 
            # and
            # [ori] ---> [brosc], respectively
            elif chosen_node.label == "brosc":
                # must add the evolution leading up to the brosc_node to the root node edge length
                # (note that the brosc_node edge length will always be 0.0 if it resulted from an
                # ancestor sampling event, but it could be > 0.0 if a state transition event happened)
                # root_node.edge_length += chosen_node.edge_length # I think this is wrong,
                
                root_node.edge_length = chosen_node.edge_length
                
                # because the extend_all_lineages() method does
                # extend brosc_node
                root_node.is_sa_lineage = chosen_node.is_sa_lineage

                # if true, means brosc's parent is a dummy node,
                # and brosc's sister taxon is a SA
                if root_node.is_sa_lineage:
                    sa_lineage_dict["root"] = sa_lineage_dict.pop("brosc")
                
                # could be origin, or dummy node
                # (in which case, brosc would have a SA sister taxon)
                brosc_parent_node = chosen_node.parent_node 
                
                # brosc is obliterated
                brosc_parent_node.remove_child(chosen_node)
                
                # we add the root node in its place
                brosc_parent_node.add_child(root_node) 

                if "brosc" in state_transition_dict:
                    state_transition_dict["root"] = state_transition_dict["brosc"]
                    del state_transition_dict["brosc"]
       
            # now make chosen_node become the root,
            # so the rest of the birth event can be executed
            chosen_node = root_node

        # normal birth event code (assuming only bifurcations)
        cumulative_node_count += 1
        left_label = "nd" + str(cumulative_node_count)
        left_child = dp.Node(
            taxon=dp.Taxon(label=left_label),
            label=left_label,
            edge_length=0.0
        )
        left_child.is_sa = False
        left_child.is_sa_lineage = False
        left_child.is_sa_dummy_parent = False
        left_child.alive = True
        left_child.state = left_arriving_state
        left_child.annotations.add_bound_attribute("state")
        state_representation_dict[left_child.state].add(
            left_child.label
        )

        cumulative_node_count += 1
        right_label = "nd" + str(cumulative_node_count)
        right_child = dp.Node(
            taxon=dp.Taxon(label=right_label),
            label=right_label,
            edge_length=0.0
        )
        right_child.is_sa = False
        right_child.is_sa_lineage = False
        right_child.is_sa_dummy_parent = False
        right_child.alive = True
        right_child.state = right_arriving_state
        right_child.annotations.add_bound_attribute("state")
        state_representation_dict[right_child.state].add(
            right_child.label
        )

        ####################################
        # Adding state transition instance #
        # to class member that bookkeeps   #
        # (for painting mid-branch)        #
        ####################################
        for idx, child_node in enumerate((left_child, right_child)):
            arriving_state: int=-1
            
            if idx == 0:
                arriving_state = left_arriving_state

            elif idx == 1:
                arriving_state = right_arriving_state

            if arriving_state != chosen_node.state:
                state_trans = AttributeTransition(
                    "state",
                    child_node.label,
                    event_t,
                    chosen_node.state,
                    arriving_state
                )

                try:
                    state_transition_dict[child_node.label].append(state_trans)
                
                except:
                    state_transition_dict[child_node.label] = [state_trans]

        # updating parent node
        # (1) adding both children
        chosen_node.add_child(left_child)
        tr_namespace.add_taxon(left_child.taxon)
        chosen_node.add_child(right_child)
        tr_namespace.add_taxon(right_child.taxon)
        
        # (2) parent node does not represent its state and cannot be targeted anymore
        try:
            state_representation_dict[chosen_node.state].remove(chosen_node.label) # state of parent is now irrelevant
        
        except:
            # if chosen_node is root, it has never been added to state_representation_dict, so we pass
            pass
        
        chosen_node.alive = False
        untargetable_node_set.add(chosen_node.label) # cannot pick parent anymore!
        
        # (3) mark parent node as the last node to undergo event
        #     will use this node to undo last event when number of species
        #     go beyond maximum
        last_node2speciate = chosen_node
        
        # (4) if chosen node was on a lineage with SAs, we update the SAs info
        if chosen_node.is_sa_lineage:
            self.update_sa_lineage_dict(event_t, sa_lineage_dict, [chosen_node.label], debug=debug)

        return last_node2speciate, cumulative_node_count


    def execute_death(self,
        tr_namespace: dp.TaxonNamespace,
        chosen_node: dp.Node,
        state_representation_dict: ty.Dict[int, ty.Set[str]],
        sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
        untargetable_node_set: ty.Set[dp.Node],
        event_t: float,
        debug=False) -> dp.Node:
        """Execute lineage death (side-effect and return)

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo extinction.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            sa_lineage_dict (dict): Dictionary that keeps track of nodes that have sampled ancestor children..
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            event_t (float): Time of death event taking place.
            debug (bool): If 'true', prints debugging messages. Defaults to False.
        """

        if debug:
            print("EXTINCTION of node " + chosen_node.label)

        # with this node gone, this state
        # is not represented by it anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        # we also cannot choose this extinct node
        # anymore to undergo an event
        chosen_node.alive = False
        untargetable_node_set.add(chosen_node.label)

        # if chosen node was on a lineage with SAs,
        # we update the SAs info
        if chosen_node.is_sa_lineage:
            self.update_sa_lineage_dict(event_t, sa_lineage_dict, [chosen_node.label])

        ######################################
        # Special case: origin went extinct, #
        # we slap a brosc node               #
        ######################################
        if chosen_node.label == "origin":
            # at this point, the origin was chosen to die, but this node will never be extended
            # (the origin always has an origin_edge_length = 0.0); for us to account for the
            # evolution (branch length) that has happened before this death -- between the origin
            # and the brosc_node being added -- we must add event_t as the brosc_node edge_length
            # (other nodes are always added with edge_length = 0.0, and have their edges extended
            # when a new event takes place)
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"), label="brosc", edge_length=event_t)
            brosc_node.alive = False
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False
            brosc_node.state = chosen_node.state
            untargetable_node_set.add(brosc_node.label)
            chosen_node.add_child(brosc_node)
            tr_namespace.add_taxon(brosc_node.taxon)

        last_node2die = chosen_node

        return last_node2die


    def execute_anatrans(self,
                        tr_namespace: dp.TaxonNamespace,
                        chosen_node: dp.Node,
                        state_representation_dict: ty.Dict[int, ty.Set[str]],
                        state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                        untargetable_node_set: ty.Set[str],
                        macroevol_rate_param: sseobj.DiscreteStateDependentRate,
                        event_t: float,
                        debug: bool=False) -> None:
        """
        Execute anagenetic trait-state transition on path to chosen
        node (side-effect)

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object
                recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo
                anagenetic trait-state transition.
            state_representation_dict (dict): Dictionary that keeps
                track of all states currently represented by living
                lineages.
            state_transition_dict (dict): Dictionary that keeps track
                of nodes subtending which state transitions happen
                (used for plotting)
            macroevol_rate_param (AtomicRateParameter): Atomic rate
                parameter containing departing/arriving state.
            event_t (float): Time of anagenetic transition event taking
                place.
            debug (bool): If 'True', prints debugging messages.
                Defaults to False.
        """

        if debug:
            print("TRANSITION of node " + chosen_node.label + " from state " + str(chosen_node.state) + " to state " + str(macroevol_rate_param.arriving_states[0]))

        # new state
        arriving_state = macroevol_rate_param.arriving_states[0]

        # old state is not represented anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        ########################################################
        # Special case: anagenetic trasition event before root #
        ########################################################
        if chosen_node.label == "origin":
            # at this point, the origin was chosen to undergo anagenetic state transition, but this
            # node will never become extended (the origin always has an origin_edge_length = 0.0);
            # for us to account for the evolution (branch length) that has happened before this death
            # -- between the origin and the brosc_node being added -- we must add event_t as the
            # brosc_node edge_length (other nodes are always added with edge_length = 0.0, and have
            # their edges extended when a new event takes place)
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"), label="brosc", edge_length=event_t)
            brosc_node.alive = True
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False
            brosc_node.state = chosen_node.state
            tr_namespace.add_taxon(brosc_node.taxon)

            # origin cannot be selected anymore
            chosen_node.alive = False # origin node is no longer alive
            untargetable_node_set.add(chosen_node.label) # cannot be targeted
            chosen_node.add_child(brosc_node) # brosc is now child of origin

            # now make chosen_node become brosc, so the rest of the anagenetic transition event can be executed
            chosen_node = brosc_node

        ####################################
        # Adding state transition instance #
        # to class member that bookkeeps   #
        # (for painting mid-branch)        #
        ####################################
        state_trans = AttributeTransition(
            "state",
            chosen_node.label,
            event_t,
            chosen_node.state,
            arriving_state)

        try:
            state_transition_dict[chosen_node.label].append(state_trans)
        
        except:
            state_transition_dict[chosen_node.label] = [state_trans]

        # last updates and bookkeping #
        # new state gets assigned
        chosen_node.state = arriving_state

        # new state is now represented
        state_representation_dict[chosen_node.state].add(chosen_node.label)


    def execute_sample_ancestor(self,
                                tr_namespace: dp.TaxonNamespace,
                                chosen_node: dp.Node,
                                state_representation_dict: ty.Dict[int, ty.Set[str]],
                                sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                                untargetable_node_set,
                                cumulative_sa_count: int,
                                event_t: float,
                                debug: bool=False) -> int:
        """Execute sampling of direct lineage ancestor (side-effect and return)

        Args:
            tr_namespace (dendropy.TaxonNamespace): 
            chosen_node (dendropy.Node): Node that will undergo event.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            sa_lineage_dict (dict): Dictionary that keeps track of nodes that have sampled ancestor children.
            cumulative_sa_count (int): Total number of sampled ancestors nodes in the tree (to be used in labeling).
            event_t (float): Time of ancestor sampling event taking place.
            debug (bool): If 'true', prints debugging messages. Defaults to False.
        """
        
        if debug:
            print("> ANCESTOR-SAMPLING of node " + chosen_node.label + " , keeping state " + str(chosen_node.state))
        
        # modeling sampled ancestor as bifurcating event, one child is the same as parent lineage,
        # second child is a terminal node that will never be elongated (branch length = 0.0);
        # the original parent is made into a dummy node

        if chosen_node.label == "origin":
            # ancestor sampling event took place before root was born so we create a replacement
            # for the root, the 'brosc_node' before-root-origin-single-child; this is the node that
            # will undergo the ancestor sampling event
            #
            # the origin was chosen to undergo the ancestor sampling, but this node will never become
            # extended (the origin always has an origin_edge_length = 0.0); for us to account for the
            # evolution (branch length) that has happened before this death -- between the origin and
            # the brosc node being added -- we must add event_t as the brosc node edge_length (other
            # nodes are always added with edge_length = 0.0, and have their edges extended when a new
            # event takes place)
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"), label="brosc", edge_length=event_t)
            brosc_node.alive = True
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = True
            brosc_node.state = chosen_node.state
            brosc_node.annotations.add_bound_attribute("state")
            state_representation_dict[chosen_node.state].add(brosc_node.label) # this node is alive and can be selected
            tr_namespace.add_taxon(brosc_node.taxon)
            
            # updating origin
            chosen_node.alive = False # origin is no longer alive
            untargetable_node_set.add("origin") # cannot be targeted
            state_representation_dict[chosen_node.state].remove(chosen_node.label) # and does not represent state
            chosen_node.add_child(brosc_node)

            chosen_node = brosc_node # we will execute the ancestor sampling on this new chosen_node

            # print("\n\nMust have a sampled ancestor before the root!")

        # doing lineage that remains alive
        parent_label = chosen_node.label # gets parent's label (parent's label will change below)
        left_child = dp.Node(taxon=dp.Taxon(label=parent_label), label=parent_label, edge_length=0.0)
        left_child.is_sa = False
        left_child.is_sa_lineage = True
        left_child.is_sa_dummy_parent = False
        left_child.alive = True
        left_child.state = chosen_node.state # gets parent's state
        left_child.annotations.add_bound_attribute("state")
        # we don't update state_representation_dict, because parent's label is already there

        # doing sampled ancestor
        cumulative_sa_count += 1
        right_label = "sa" + str(cumulative_sa_count)
        right_child = dp.Node(taxon=dp.Taxon(label=right_label), label=right_label, edge_length=0.0)
        right_child.is_sa = True
        right_child.is_sa_lineage = False
        right_child.is_sa_dummy_parent = False
        right_child.alive = False
        right_child.state = chosen_node.state # gets parent's state
        right_child.annotations.add_bound_attribute("state")
        # state_representation_dict[right_child.state].add(right_child.label) # I don't think sampled ancestor represents state
        # the way we're implementing sampled ancestors here, they cannot be picked for events;
        # it is their sister lineage that can be picked (their parent is becomes a dummy node)
        untargetable_node_set.add(right_child.label) # cannot pick sampled ancestor to undergo event the way its implemented

        ####################################################
        # Adding sampled ancestor and its lineage node     #
        # to class member that stashes it, or initializing #
        # class member                                     #
        ####################################################
        sa = SampledAncestor(right_label, parent_label, event_t)
        try:
            sa_lineage_dict[parent_label].append(sa)
        except:
            sa_lineage_dict[parent_label] = [sa]

        # updating parent node
        # (1) adding both children
        chosen_node.add_child(left_child)
        tr_namespace.add_taxon(left_child.taxon)
        chosen_node.add_child(right_child)
        tr_namespace.add_taxon(right_child.taxon)
        # (2) update parent's label
        chosen_node.label = "dummy" + str(cumulative_sa_count)
        # note that here, we do not remove parent's label from targeted nodes,
        # because left child will now have that label and should be targetable
        # (3) parent label cannot be targeted anymore
        chosen_node.alive = False
        chosen_node.is_sa_dummy_parent = True
        chosen_node.is_sa_lineage = False # left child will always be the sa_lineage
        untargetable_node_set.add(chosen_node.label) # cannot pick parent!

        return cumulative_sa_count


    def update_sa_lineage_dict(self,
                                a_time: float,
                                sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                                sa_lineage_node_labels: ty.List[str], 
                                debug: bool=False) -> None:
        """Update sa_lineage_dict (side-effect) when lineage node undergoes event and at tree stop condition
        
        This function is called every time a node (whose subtending branch has sampled ancestors, by asking if
        the node .is_sa_lineage == True) undergoes an event, and when the tree reaches its stop condition.
        
        When this function is called because a single node undergoes an event, then 'sa_lineage_node_label' will
        have a single element. When this function is called because the process reached its stop condition,
        then 'sa_lineage_node_label' will have the labels of all terminal nodes for which .is_sa_lineage == True.

        By updating the dictionary as the tree is built, we don't have to do tree traversals later to get
        sampled ancestor times for plotting (when initializing an AnnotatedTree).

        Args:
            a_time (float): Either the time of the last event an SA lineage node undergoes, or the simulation end time.
            sa_lineage_dict (dict): Dictionary that keeps track of nodes that have sampled ancestor children..
            sa_lineage_node_labels (str): List of node labels, whose subtending branch has SAs.
        """

        if debug:
            print("\n>> Updating SA lineage dictionary. Keys are is \'" + ", ".join(sa_label for sa_label in sa_lineage_node_labels) + \
                  "\' and time being added is " + str(a_time) + "\n   SA lineage dict:")
            print("      " + "\n      ".join(k + " has in its lineage " + ", ".join(sa.label for sa in sas) for k, sas in sa_lineage_dict.items()))

        for sa_lineage_node_label in sa_lineage_node_labels:
            sa_list = sa_lineage_dict[sa_lineage_node_label]
            
            for sa in sa_list:
                sa.time_to_lineage_node = a_time - sa.global_time
        # if isinstance(sa_lineage_node_label, str):
        #     sa_list = sa_lineage_dict[sa_lineage_node_label]
            
        #     for sa in sa_list:
        #         sa.time_to_lineage_node = a_time - sa.global_time

        # for _, sa_list in sa_lineage_dict.items():
        #     for sa in sa_list:
        #         sa.time_to_lineage_node = a_time - sa.global_time


    def execute_event(self,
                        tr_namespace,
                        macroevol_rate_param: sseobj.DiscreteStateDependentRate,
                        chosen_node: dp.Node,
                        state_representation_dict: ty.Dict[int, ty.Set[str]],
                        state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]],
                        sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]],
                        untargetable_node_set: ty.Set[str],
                        cumulative_node_count: int,
                        cumulative_sa_count: int,
                        last_chosen_node,
                        event_t: float,
                        debug: bool=False) -> ty.Tuple[dp.Node, int, int]:
        """Execute event on chosen node and bookkeep things

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object recording taxa in the tree.
            macroevol_rate_param (DiscreteStateDependentRate): Instance of DiscreteStateDependentRate carrying event information.
            chosen_node (dendropy.Node): Node that will undergo event.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            sa_lineage_dict (dict): Dictionary that keeps track of nodes that have sampled ancestor children..
            state_transition_dict (dict): Dictionary that keeps track of nodes below which state transitions happen (used for plotting).
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            cumulative_sa_count (int): Total number of sampled ancestor nodes in the tree (to be used in labeling).
            last_chosen_node (dendropy.Node): Last node undergoing event.
            event_t (float): Time of event taking place.
            debug (bool): If 'true', prints debugging messages. Defaults to False.

        Returns:
            (dendropy.Node, int): Tuple with last node to under go event and total (cumulative) node count.
        """
        macroevol_event = macroevol_rate_param.event

        if macroevol_event == sseobj.MacroevolEvent.W_SPECIATION or \
            macroevol_event == sseobj.MacroevolEvent.BW_SPECIATION or \
            macroevol_event == sseobj.MacroevolEvent.ASYM_SPECIATION:
            
            last_chosen_node, cumulative_node_count = \
                self.execute_birth(
                    tr_namespace,
                    chosen_node,
                    state_representation_dict,
                    sa_lineage_dict,
                    state_transition_dict,
                    untargetable_node_set,
                    cumulative_node_count,
                    macroevol_rate_param,
                    event_t,
                    debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.EXTINCTION:
            self.execute_death(
                tr_namespace,
                chosen_node,
                state_representation_dict,
                sa_lineage_dict,
                untargetable_node_set,
                event_t,
                debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.ANAGENETIC_TRANSITION:
            self.execute_anatrans(
                tr_namespace,
                chosen_node,
                state_representation_dict,
                state_transition_dict,
                untargetable_node_set,
                macroevol_rate_param,
                event_t,
                debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.ANCESTOR_SAMPLING:
            cumulative_sa_count = \
                self.execute_sample_ancestor(
                    tr_namespace,
                    chosen_node,
                    state_representation_dict,
                    sa_lineage_dict,
                    untargetable_node_set,
                    cumulative_sa_count,
                    event_t,
                    debug=debug)

        return last_chosen_node, cumulative_node_count, cumulative_sa_count

    ##########################
    # Main simulation method #
    ##########################

    def simulate(self, a_start_state: int, a_stop_value: ty.Union[int, float], value_idx: int=0, a_seed: ty.Optional[int]=None) -> AnnotatedTree:
        """Sample (simulate) tree with states at its terminal nodes]

        Args:
            a_start_state (int): State at seed node (origin or root)
            a_stop_value (float): Value to stop simulation with (number of tips or tree height)
            value_idx (int):

        Returns:
            (dendropy.Tree): One simulated tree.
        """

        def debug_print_state_dict(state_dict):
            for k, v in state_dict.items():
                print("  state " + str(k) + ": " + ", ".join(v))

        def extend_all_living_nodes(branch_length, end=False):
            for nd in tr:
                # the root node is extended (root_edge is made > 0.0) here too, when the root is picked
                if nd.label not in untargetable_node_set and nd.label != "origin" and nd.alive:
                    nd.edge_length += branch_length

        ### START initializing values pre-simulation ###

        # time variables
        # variables with "t" are times (0.0 is origin or root)
        t_stop = 0.0
        latest_t = 0.0 # will increase as simulation progresses

        max_obs_nodes = 0
        
        current_node_target_count = 1
        cumulative_node_count = 1
        cumulative_sa_count = 0
        
        untargetable_node_set = set() # not really used for anything, .alive takes care of everything
        reached_stop_condition = False
        start_state = a_start_state
        state_representation_dict: ty.Dict[int, ty.Set[str]] = dict((i, set()) for i in range(self.events.state_count)) # { 0: ["origin", "root"], 1: ["nd1", "nd2"], 2:["nd3",...], ... }
        last_chosen_node = None
        tr = dp.Tree()
        tree_invalid = False
        tree_died = False

        state_transition_dict: ty.Dict[str, ty.List[AttributeTransition]] = dict() # tree info for plotting
        sa_lineage_dict: ty.Dict[str, ty.List[SampledAncestor]] = dict() # tree info for plotting

        if self.stop == "age" and isinstance(a_stop_value, float):
            t_stop = a_stop_value

        elif self.stop == "size" and isinstance(a_stop_value, int):
            max_obs_nodes = a_stop_value

        # simulation starts at origin
        if self.with_origin:
            ############################################################################
            # brosc: before-root-origin-single-child                                   #
            # Node that might be necessary for tying process off if a speciation never #
            # happens before the process either dies, or reaches stop condition        #
            #                                                                          # 
            # Note: brosc can become an sa lineage node if ancestor sampling happens,  #
            # but a root is never born                                                 #
            ############################################################################
            
            # this particular brosc_node instance below will be used within the simulate() method,
            # but not inside "execute_sample_ancestor()", which has its own instance of brosc_node
            # being created 
            brosc_node = dp.Node(taxon=dp.Taxon(label="brosc"), label="brosc", edge_length=0.0)
            brosc_node.state = start_state
            brosc_node.annotations.add_bound_attribute("state")
            brosc_node.is_sa = False
            brosc_node.is_sa_dummy_parent = False
            brosc_node.is_sa_lineage = False

            # origin node will be the seed_node
            origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
            origin_node.state = start_state
            origin_node.annotations.add_bound_attribute("state")
            origin_node.alive = True
            origin_node.is_sa = False
            origin_node.is_sa_dummy_parent = False
            origin_node.is_sa_lineage = False
            state_representation_dict[origin_node.state].add(origin_node.label)

            # now make tree
            tr = dp.Tree(seed_node=origin_node) # will remain a "is_leaf() == True" until we add children
            tr.taxon_namespace.add_taxon(origin_node.taxon)

        # simulation starts at root
        else:
            # root node will be the seed_node, and will be untargetable from the get-go
            root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.0)
            root_node.state = start_state
            root_node.annotations.add_bound_attribute("state")
            root_node.alive = False
            root_node.is_sa = False
            root_node.is_sa_dummy_parent = False
            root_node.is_sa_lineage = False
            untargetable_node_set.add("root")
            self.root_is_born = True

            # we will actually start off with the left and right nodes being targetable
            left_node = dp.Node(taxon=dp.Taxon(label="nd1"), label="nd1", edge_length=0.0)
            left_node.is_sa = False
            left_node.is_sa_dummy_parent = False
            left_node.is_sa_lineage = False
            left_node.alive = True
            left_node.state = root_node.state # get state from parent (root)
            left_node.annotations.add_bound_attribute("state")
            root_node.add_child(left_node)
            state_representation_dict[left_node.state].add(left_node.label)
            cumulative_node_count += 1

            right_node = dp.Node(taxon=dp.Taxon(label="nd2"), label="nd2", edge_length=0.0)
            right_node.is_sa = False
            right_node.is_sa_dummy_parent = False
            right_node.is_sa_lineage = False
            right_node.alive = True
            right_node.state = root_node.state # get state from parent (root)
            right_node.annotations.add_bound_attribute("state")
            root_node.add_child(right_node)
            state_representation_dict[right_node.state].add(right_node.label)
            cumulative_node_count += 1

            # now make tree
            tr = dp.Tree(seed_node=root_node) # will remain a "is_leaf() == True" until we add children
            tr.taxon_namespace.add_taxon(root_node.taxon)
            tr.taxon_namespace.add_taxon(left_node.taxon)
            tr.taxon_namespace.add_taxon(right_node.taxon)

        # need these counts initialized for sampling first event
        living_nodes = [nd for nd in tr if nd.alive]
        current_node_target_count = len([nd for nd in tr if nd.alive])

        # debugging
        # if self.debug:
        #     print("\nLiving nodes: " + ", ".join(nd.label for nd in living_nodes))
        #     print("State representation:")
        #     debug_print_state_dict(state_representation_dict)

        ### END initializing values pre-simulation ###

        ### START simulation loop ###
        time_slice_idx = 0
        while (time_slice_idx < self.n_time_slices and not reached_stop_condition):

            # (1) Find out the overall total rate given the current targetable nodes' states
            rate_for_exponential_distn, state_total_rates = self.events.total_rate(latest_t,
                                                                    state_representation_dict,
                                                                    value_idx=value_idx,
                                                                    debug=self.debug) # [1] are all states with >= 1 representation

            # (2) get the time to the next event
            t_to_next_event = self.get_next_event_time(rate_for_exponential_distn)
            latest_t += t_to_next_event

            # (3) get end time of this slice (user provides it as end ages,
            # but FIGManager converts it to time ends)
            _next_max_t: float
            _t_end = self.slice_t_ends[time_slice_idx]
            if isinstance(_t_end, float):
                _next_max_t = _t_end

            # (4) check if new event time cut throw the end of a time slice, provided
            # there is more than 1 time slice, and that a max age was specified
            # through the "age stop condition"
            # 
            # next_max_t will be None if self.slice_t_ends is empty
            excess_t = 0.0
            if self.stop == "age" and self.n_time_slices > 1 and latest_t > _next_max_t:
                excess_t = latest_t - _next_max_t
                latest_t = _next_max_t
                time_slice_idx += 1

                # extend all lineages (could keep just the extend in the else-block
                # below) and keep everything up-to-date
                extend_all_living_nodes(t_to_next_event - excess_t)

                if _next_max_t == t_stop:
                    sa_lineage_node_labels = [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self.update_sa_lineage_dict(t_stop, sa_lineage_dict, sa_lineage_node_labels) # updates SA info for plotting
                    
                    # if origin is the only node (root always has children), we slap finish node at end of process
                    if self.with_origin and tr.seed_node.alive and len(tr.seed_node.child_nodes()) == 0:
                        brosc_node.edge_length = t_stop # we make origin edge length the max age of the tree
                        brosc_node.alive = True
                        tr.seed_node.add_child(brosc_node)
                        tr.taxon_namespace.add_taxon(brosc_node.taxon)

                    reached_stop_condition = True
                    
                    break

                # goes back to top of while-loop restarting at next time slice
                continue

            # (5) new event time did not go over a time slice, all good, we will now
            # draw an event and execute it
            else:
                # check if tree grew beyond max age (need this check here because the
                # max age check in the if-block above will not execute
                # with a single time slice)
                if self.stop == "age" and (latest_t > t_stop):
                    extend_all_living_nodes(t_stop - (latest_t - t_to_next_event), end=True)
                    
                    sa_lineage_node_labels = [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self.update_sa_lineage_dict(t_stop, sa_lineage_dict, sa_lineage_node_labels) # updates SA info for plotting

                    # if origin is the only node (root always has children), we slap brosc node at end of process
                    if self.with_origin and tr.seed_node.alive and len(tr.seed_node.child_nodes()) == 0:
                        brosc_node.edge_length = t_stop
                        brosc_node.alive = True
                        tr.seed_node.add_child(brosc_node)
                        tr.taxon_namespace.add_taxon(brosc_node.taxon)
                    
                    reached_stop_condition = True
                    
                    break

                # elongate all lineages
                extend_all_living_nodes(t_to_next_event) # original implmn w/o slices

                # (6) draw a node we'll apply the event to
                #
                # a lineage will be chosen in proportion to the total rate of its state
                lineage_weights = [state_total_rates[nd.state] for nd in living_nodes]
                chosen_node = random.choices(living_nodes, weights=lineage_weights)[0]

                # debugging
                if self.debug:
                    print("\nChosen node: " + chosen_node.label + " in state " + str(chosen_node.state) + "\n")

                # (7) draw an event
                macroevol_event_state = chosen_node.state
                state_conditioned_total_rate = state_total_rates[macroevol_event_state] # denominator for sampling event
                macroevol_atomic_param_in_list = self.events.sample_event_atomic_parameter(
                    state_conditioned_total_rate, latest_t, [macroevol_event_state], value_idx=value_idx, debug=self.debug) # total_rate and latest_t have been updated
                macroevol_atomic_param = macroevol_atomic_param_in_list[0]

                # (8) execute event
                last_chosen_node, cumulative_node_count, cumulative_sa_count = \
                    self.execute_event(
                        tr.taxon_namespace,
                        macroevol_atomic_param,
                        chosen_node,
                        state_representation_dict,
                        state_transition_dict,
                        sa_lineage_dict,
                        untargetable_node_set,
                        cumulative_node_count,
                        cumulative_sa_count,
                        last_chosen_node,
                        latest_t,
                        debug=self.debug
                    )

                # (9) update number of lineages after event
                #
                # what counts for targetable node here are living terminal nodes
                # TODO: at some point, make 'living_nodes' a set that is passed to all execute_X functions
                # so that we don't need to traverse the tree every event to list living nodes that can be 
                # targeted
                living_nodes = [nd for nd in tr if nd.alive]
                current_node_target_count = len(living_nodes)

                # (10) check for stop conditions
                #
                # check if all lineages went extinct, stop!
                if current_node_target_count == 0:
                    tree_died = True

                    sa_lineage_node_labels = [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self.update_sa_lineage_dict(latest_t, sa_lineage_dict, sa_lineage_node_labels) # updates SA info for plotting
                    
                    reached_stop_condition = True

                    break

                # [Condition is "age", but rejection sampling
                # was requested for taxon count range]
                # 
                # this helps stop trees growing out of control
                # within age limit
                if current_node_target_count == self.abort_at_obs:
                    tree_invalid = True
                    reached_stop_condition = True

                    break


                # check tree grew beyond its max taxon count, stop!
                #
                # it's '>' max_obs_nodes below, instead of '==' or '>='
                # because we want to have one extra species (that then gets removed),
                # i.e., we want the longest possible tree up to max_obs_nodes
                # by waiting until one-too-many speciation events
                if self.stop == "size" and current_node_target_count > max_obs_nodes:                    
                    last_chosen_node.clear_child_nodes() # delete last created children
                    last_chosen_node.alive = True # parent is back alive
                    untargetable_node_set.remove(last_chosen_node.label) # we must add parent back so it can be extended
 
                    if max_obs_nodes == 1:
                        # special case: birth happens on just-added root, and we want just a single lineage;
                        # (only makes sense if starting at origin, if starting at root, two lineages already exist)
                        # we will replace the root node with a brosc node
                        if self.with_origin:
                            # getting root attributes
                            brosc_node.state = last_chosen_node.state
                            brosc_node.edge_length = last_chosen_node.edge_length
                            brosc_node.alive = last_chosen_node.alive # should be True
                            last_chosen_node.parent_node.add_child(brosc_node) # connecting brosc to origin
                            tr.taxon_namespace.add_taxon(brosc_node.taxon)
                            last_chosen_node.parent_node.remove_child(last_chosen_node) # removing root
                            last_chosen_node = brosc_node
                    
                    state_representation_dict[last_chosen_node.state].add(last_chosen_node.label) # won't use it again, but to be safe
                    
                    sa_lineage_node_labels = [nd.label for nd in living_nodes if nd.is_sa_lineage]
                    self.update_sa_lineage_dict(t_stop, sa_lineage_dict, sa_lineage_node_labels) # updates SA info for plotting

                    reached_stop_condition = True

                    break

        # (11) got out of while loop because met stop condition
        # 'at' is scoped to simulate_a_tree() function
        
        # debugging
        # for nd_name, at_list in state_transition_dict.items():
        #     print("nd_name = " + nd_name)
        #     for at in at_list:
        #         print(str(at))

        if self.stop == "age":
            at = AnnotatedTree(tr,
                               self.events.state_count,
                               start_at_origin=self.with_origin,
                               max_age=a_stop_value,
                               slice_t_ends=self.slice_t_ends,
                               slice_age_ends=self.events.slice_age_ends,
                               sa_lineage_dict=sa_lineage_dict,
                               at_dict=state_transition_dict,
                               condition_on_obs_both_sides_root=self.condition_on_obs_both_sides_root,
                               tree_invalid=tree_invalid,
                               tree_died=tree_died)
        
        elif self.stop == "size":
            at = AnnotatedTree(tr,
                               self.events.state_count,
                               start_at_origin=self.with_origin,
                               sa_lineage_dict=sa_lineage_dict,
                               at_dict=state_transition_dict,
                               condition_on_obs_both_sides_root=self.condition_on_obs_both_sides_root,
                               tree_invalid=tree_invalid)

        if self.debug:
            # print(at.tree)
            # print(at.tree.as_string(schema="nexus", suppress_annotations=True, suppress_internal_taxon_labels=True))
            print(at.tree.as_string(schema="newick", suppress_annotations=False, suppress_internal_taxon_labels=True))

        return at
    ### END simulation loop ###


    def generate(self) -> ty.List[AnnotatedTree]:
        # output
        output: ty.List[AnnotatedTree] = []

        # do something with self.runtime_limit
        start_time = time.time()
        ith_sim = 0
        j = 0
        while len(output) < (self.n_sim * self.n_repl):
            ellapsed_time = pjh.get_ellapsed_time_in_minutes(start_time, time.time())

            if ellapsed_time >= self.runtime_limit:
                print("Aborted simulations, hit runtime limit")
                break

            # simulate!
            repl_size = 0
            while repl_size < self.n_repl:
                # print("\n\nsimulation #" + str(ith_sim))
                tr = self.simulate(self.start_states[ith_sim], self.stop_val[ith_sim], value_idx=ith_sim)

                # check if tr has right specs
                if self.is_tr_ok(tr, self.stop_val[ith_sim]):
                    output.append(tr)
                    repl_size += 1

                # tree not good, stay in while loop
                else:
                    j += 1
                    continue

            ith_sim += 1
            # print(ith_sim)

        return output

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        return rev_str_list

    #########################
    # Output health methods #
    #########################
    def is_tr_ok(self, ann_tr: AnnotatedTree, a_stop_value: ty.Union[int, float]) -> bool:
        """Check that simulated tree meets stop conditions, so that it can be returned

        Args:
            a_tree (AnnotatedTree): A tree with extra annotations
            a_stop_value [int or float]: Value specified for stopping simulation

        Returns:
            bool: If tree meets stop condition value.
        """
        # print("self.tree_died inside SSE dn: " + str(ann_tr.tree_died))
        
        if self.condition_on_survival and ann_tr.tree_died:
            return False

        if self.condition_on_speciation and not isinstance(ann_tr.root_node, dp.Node):
            return False

        # need to worry about rejection sampling
        # when user conditions reconstructed tree
        if self.condition_on_obs_both_sides_root:
            if self.condition_on_speciation and not isinstance(ann_tr.root_node, dp.Node):
                return False

            rec_tree = ann_tr.extract_reconstructed_tree()
            
            if rec_tree.__str__() == "":
                # debugging
                # print("\nRejected reconstructed tree")
                # print("    Complete tree:")
                # print("    " + ann_tr.tree.as_string(schema="newick", suppress_annotations=True, suppress_internal_node_labels=False, suppress_internal_taxon_labels=True))
                # for ch in ann_tr.tree.leaf_node_iter():
                #     print(ch.label + " alive = " + str(ch.alive))
                # print("\n")
                return False
            
            # debugging
            # print("\nSuccessful reconstructed tree:")
            # print(rec_tree.as_string(schema="newick", suppress_annotations=True, suppress_internal_node_labels=False, suppress_internal_taxon_labels=True))

        if self.stop == "size":
            if self.condition_on_speciation and isinstance(ann_tr.root_node, dp.Node) and len(ann_tr.root_node.child_nodes()) == 0:
                return False

            if self.condition_on_survival and ann_tr.n_extant_terminal_nodes != a_stop_value:
                return False

            # if conditions are not specified, returns all trees
            return True

        elif self.stop == "age":
            # grew out of control! #
            if ann_tr.tree_invalid:
                # tree was cut short of stop_condition
                # because it grew too much!
                print("Rejected: tree grew out of control!")
                #
                # starting from root
                if (not ann_tr.with_origin and isinstance(ann_tr.root_age, float)) and \
                    (a_stop_value - ann_tr.root_age) > self.epsilon:
                    return False
                # starting from origin

            # print("WARNING: don't forget to remove n-extant-tip bracketing in is_tr_ok() after FIG validation!!!")
            # TODO: remove these two ifs after done with FIG validation
            # reconstructed tree is too small
            if ann_tr.n_extant_terminal_nodes < self.min_rec_taxa:
                # print("Rejected: too small, " + str(ann_tr.n_extant_terminal_nodes))
                return False
            
            # reconstructed tree is too large
            if ann_tr.n_extant_terminal_nodes > self.max_rec_taxa:
                # print("Rejected: too large, " + str(ann_tr.n_extant_terminal_nodes))
                return False

            if ann_tr.with_origin and isinstance(ann_tr.origin_age, float):
                # tree is ok if either it reached the max age by epsilon, or if its age is smaller than max age
                if abs(a_stop_value - ann_tr.origin_age) <= self.epsilon or ann_tr.origin_age < a_stop_value:
                    # print("    tree was accepted!")
                    # print("    " + ann_tr.tree.as_string(schema="newick", suppress_annotations=True))
                    return True
                elif ann_tr.origin_age > a_stop_value:
                    # print("    tree was rejected, origin_age = " + str(ann_tr.origin_age))
                    # print("    " + ann_tr.tree.as_string(schema="newick", suppress_annotations=True))
                    return False

            elif not ann_tr.with_origin and isinstance(ann_tr.root_age, float):
                if (a_stop_value - ann_tr.root_age) <= self.epsilon:
                    return True
                else:
                    return False

        # stop condition is neither "size" nor "age", and no conditioning was violated
        # then something is off...

        return False