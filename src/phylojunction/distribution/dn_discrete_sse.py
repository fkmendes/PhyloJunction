"""dn_discrete_sse.py: Class for distribution of discrete state-dependent speciation and extinction process"""

import typing as ty
import time
import numpy
# import statistics
import random
# from random import seed, choice
import dendropy as dp
# from dendropy import Tree, Node, Taxon, TaxonNamespace
# from matplotlib.pyplot import draw

# pj imports
import pgm.pgm as pgm
import distribution.dn_parametric as dnpar
# from tpsimulator_utils import *
# from tpsimulator_math_lib import draw_exp
import calculation.discrete_sse as sseobj
import utility.helper_functions as pjh
import utility.exception_classes as ec
from data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class DnSSE(pgm.DistributionPGM):
    """Distribution class for State-dependent Speciation and Extinction (birth-death) process

    The sampling method in this class is a "rising tide" simulator,
    where all living lineages grow together.

    At all times we know:
        * How many lineages are alive and dead
        * All the character states represented by living lineages

    This is _not_ a recursive simulator (in which each lineage would take care of
    growing only itself, and recur upon birth events).

    Attributes
    ----------
        n_sim (int): Number of simulations one wants.
        n_repl (int): Number of (successful) trees per simulation.
        stop (str): If "age", stops when origin or root age is equal to 'stop_condition_value'. If "size", stops when tree has 'stop_condition_value' observable nodes.
        stop_val (float): Either maximum tree age, or maximum count of observable nodes.
        start_at_origin (bool): Simulation starts at origin (i.e., seed = origin).
        events (MacroEvolEventHandler): Object holding all parameter values we need to simulate.
        start_states (int): List of integer representing the starting states of all n_sim simulations.
        state_count (int): Number of discrete states (obtained from 'events').
        n_time_slices (int): Number of time slices (obtained from 'events').
        slice_t_ends (int): List of end times from time slices (obtained from 'events').
        condition_on_speciation (bool): If true, rejects trees that go extinct before reaching stop_condition_value. Defaults to False.
        condition_on_survival (bool): If true, rejects trees that go extinct before reaching stop_condition_value. Defaults to False.
        epsilon (float, optional): Any values difference below this value is set to 0.0. Defaults to 1e-12.
        runtime_limit (int, optional): Maximum number of minutes one is willing to wait until all 'n_sim' simulations are done. Defaults to 10.
        debug (bool, optional): Whether to print debugging messages

    Methods
    -------

    """

    DN_NAME = "DnSSE"

    def __init__(self, n=1, n_replicates=1, stop=None, stop_value=None, origin=True, event_handler=None,
                start_states_list=[], condition_on_speciation=False, condition_on_survival=False,
                seeds_list=None, epsilon=1e-12, runtime_limit=60, debug=False):

        # simulation parameters
        self.n_sim = int(n)
        self.n_repl = n_replicates # number of replicate trees (in plate) for a given simulation
        self.with_origin = origin
        self.stop = stop
        self.stop_val = stop_value
        self.condition_on_speciation = condition_on_speciation
        self.condition_on_survival = condition_on_survival

        # model parameters
        self.events = event_handler # carries all parameters, number of states and of slices
        self.start_states = start_states_list
        self.state_count = self.events.state_count
        self.n_time_slices = self.events.n_time_slices
        self.slice_t_ends = self.events.slice_t_ends
        self.seed_age = self.events.seed_age # used just for verifying inputs

        # other specs
        self.seeds = seeds_list
        self.epsilon = epsilon
        self.runtime_limit = runtime_limit
        self.debug = debug

        # dealing with vectorization
        self.check_sample_size()

        # making sure inputs are ok
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
        if self.stop == "age" and self.seed_age and self.seed_age != self.stop_val:
            raise ec.DnInitMisspec(self.DN_NAME, "The max. stopping age for the simulation must match the specified seed (origin or root) age. Exiting...")


    #########################
    # Vectorization methods #
    #########################
    def check_sample_size(self):
        # changing atomic_rate_params_matrix here should affect fig_rates_manager
        atomic_rate_params_matrix = self.events.fig_rates_manager.atomic_rate_params_matrix # 1D: time slices, 2D: list of atomic rate params

        # rate parameters
        # k-th time slice
        for k, list_atomic_rate_params in enumerate(atomic_rate_params_matrix):
            for arp in list_atomic_rate_params:
                n_val = len(arp.value)

                # we multiply values if one is provided, but > 1 sims
                if n_val == 1 and self.n_sim > 1:
                    arp.value = [arp.value[0] for i in range(self.n_sim)]
                # do not know how to multiply, error!
                elif n_val > 1 and n_val < self.n_sim:
                    raise ec.DimensionalityError(self.DN_NAME)

        # starting state
        if len(self.start_states) > self.n_sim or (len(self.start_states) > 1 and len(self.start_states) < self.n_sim):
            raise ec.DimensionalityError(self.DN_NAME)

        elif len(self.start_states) < self.n_sim and len(self.start_states) == 1:
            try:
                self.start_states = [int(self.start_states[0]) for i in range(self.n_sim)]
            except:
                self.start_states = [self.start_states for i in range(self.n_sim)]


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


    def execute_birth(self, tr_namespace: dp.TaxonNamespace, chosen_node: dp.Node, state_representation_dict,
                                    untargetable_node_set, cumulative_node_count: int, last_node2speciate: dp.Node,
                                    macroevol_atomic_param: sseobj.AtomicSSERateParameter, debug=False) -> ty.Tuple[dp.Node, int]:
        """Execute lineage birth

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object recording taxa in the tree.
            chosen_node (dendropy.Node): Node that will undergo speciation.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            last_node2speciate (dendropy.Node): Last node to undergo speciation, to be deleted when tree gets too large
            macroevol_atomic_param (AtomicRateParameter): Atomic rate parameter containing departing/arriving state.
            debug (bool): If 'true', prints debugging messages. Defaults to False.

        Returns:
            (dendropy.Node, int): Tuple with last node to under go event and total (cumulative) node count.
        """

        left_arriving_state, right_arriving_state = macroevol_atomic_param.arriving_states

        if debug:
            print("SPECIATION of node " + chosen_node.label + " in state " + str(chosen_node.state) + \
                " into daughters with states " + str(left_arriving_state) + " and " + str(right_arriving_state))

        # assume only bifurcations
        cumulative_node_count += 1
        left_label = "nd" + str(cumulative_node_count)
        left_child = dp.Node(taxon=dp.Taxon(label=left_label), label=left_label, edge_length=0.0)
        left_child.alive = True
        left_child.state = left_arriving_state
        left_child.annotations.add_bound_attribute("state")
        state_representation_dict[left_child.state].add(left_child.label)

        cumulative_node_count += 1
        right_label = "nd" + str(cumulative_node_count)
        right_child = dp.Node(taxon=dp.Taxon(label=right_label), label=right_label, edge_length=0.0)
        right_child.alive = True
        right_child.state = right_arriving_state
        right_child.annotations.add_bound_attribute("state")
        state_representation_dict[right_child.state].add(right_child.label)

        # updating parent node
        # (1) adding both children
        chosen_node.add_child(left_child)
        tr_namespace.add_taxon(left_child)
        chosen_node.add_child(right_child)
        tr_namespace.add_taxon(right_child)
        # (2) parent node does not represent its state and cannot be targeted anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label) # state of parent is now irrelevant
        chosen_node.alive = False
        untargetable_node_set.add(chosen_node.label) # cannot pick parent anymore!
        # (3) mark parent node as the last node to under go even
        #     will use this node to undo last event when number of species
        #     go beyond maximum
        last_node2speciate = chosen_node

        return last_node2speciate, cumulative_node_count

    def execute_death(self, chosen_node, state_representation_dict, untargetable_node_set, debug=False):
        """Execute lineage death (side-effect)

        Args:
            chosen_node (dendropy.Node): Node that will undergo extinction.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            debug (bool): If 'true', prints debugging messages. Defaults to False.
        """

        if debug:
            print("EXTINCTION of node " + chosen_node.label)

        # with this node gone, and this state is not represented by it anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        # we also cannot choose this extinct node anymore to undergo an event
        chosen_node.alive = False
        untargetable_node_set.add(chosen_node.label)

        last_node2die = chosen_node

        return last_node2die


    def execute_anatrans(self, chosen_node, state_representation_dict, macroevol_atomic_param, debug=False):
        """Execute anagenetic trait-state transition on path to chosen node (side-effect)

        Args:
            chosen_node (dendropy.Node): Node that will undergo anagenetic trait-state transition.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            macroevol_atomic_param (AtomicRateParameter): Atomic rate parameter containing departing/arriving state.
            debug (bool): If 'true', prints debugging messages. Defaults to False.
        """

        if debug:
            print("TRANSITION of node " + chosen_node.label + " from state " + str(chosen_node.state) + " to state " + str(macroevol_atomic_param.arriving_states[0]))

        # old state is not represented anymore
        state_representation_dict[chosen_node.state].remove(chosen_node.label)

        # new state gets assigned
        arriving_state = macroevol_atomic_param.arriving_states[0]
        chosen_node.state = arriving_state

        # new state is now represented
        state_representation_dict[chosen_node.state].add(chosen_node.label)


    def execute_event(self, tr_namespace, macroevol_atomic_param, chosen_node, state_representation_dict, untargetable_node_set, cumulative_node_count, last_chosen_node, debug=False):
        """Execute event on chosen node and bookkeep things

        Args:
            tr_namespace (dendropy.TaxonNamespace): Dendropy object recording taxa in the tree.
            macroevol_event (MacroevolEvent): Class specifying event that was drawn.
            chosen_node (dendropy.Node): Node that will undergo event.
            state_representation_dict (dict): Dictionary that keeps track of all states currently represented by living lineages.
            untargetable_node_set (set): Set of Node labels that cannot be targeted for events anymore (went extinct).
            cumulative_node_count (int): Total number of nodes in the tree (to be used in labeling).
            debug (bool): If 'true', prints debugging messages. Defaults to False.

        Returns:
            (dendropy.Node, int): Tuple with last node to under go event and total (cumulative) node count.
        """
        macroevol_event = macroevol_atomic_param.event

        # TODO: this should be come if macroevol_event == MacroevolEvent.SPECIATION or BW-SPECIATION or ASYM-SPECIATION
        #       and then execute_birth deals with them (avoiding code redundancy)
        if macroevol_event == sseobj.MacroevolEvent.W_SPECIATION or macroevol_event == sseobj.MacroevolEvent.BW_SPECIATION or macroevol_event == sseobj.MacroevolEvent.ASYM_SPECIATION:
            # print("node " + chosen_node.label + " split")
            last_chosen_node, cumulative_node_count = self.execute_birth(tr_namespace, chosen_node, state_representation_dict, untargetable_node_set, cumulative_node_count, last_chosen_node, macroevol_atomic_param, debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.EXTINCTION:
            # print("node " + chosen_node.label + " died")
            self.execute_death(chosen_node, state_representation_dict, untargetable_node_set, debug=debug)

        elif macroevol_event == sseobj.MacroevolEvent.ANAGENETIC_TRANSITION:
            self.execute_anatrans(chosen_node, state_representation_dict, macroevol_atomic_param, debug=debug)
            pass

        return last_chosen_node, cumulative_node_count


    def simulate(self, a_start_state, value_idx=0, a_seed=None):
        """Sample (simulate) tree with states at its terminal nodes]

        Args:
            a_start_state (int): State at seed node (origin or root).

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
        untargetable_node_set = set()
        reached_stop_condition = False
        start_state = a_start_state
        state_representation_dict = dict((i, set()) for i in range(self.events.state_count)) # { 0: ["origin", "root"], 1: ["nd1", "nd2"], 2:["nd3",...], ... }
        last_chosen_node = None
        tr = dp.Tree()

        if self.stop == "age":
            t_stop = float(self.stop_val)

        elif self.stop == "size":
            max_obs_nodes = int(self.stop_val)

        # simulation starts at origin
        if self.with_origin:
            # origin node will be the seed_node
            origin_node = dp.Node(taxon=dp.Taxon(label="origin"), label="origin", edge_length=0.0)
            origin_node.state = start_state
            origin_node.annotations.add_bound_attribute("state")
            origin_node.alive = False
            untargetable_node_set.add("origin")

            # creating root already upon initialization for practical purposes (at the moment root edge = 0.0)
            root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.0)
            root_node.alive = True # we will pick the root as the first event
            root_node.state = origin_node.state
            root_node.annotations.add_bound_attribute("state")
            state_representation_dict[root_node.state].add(root_node.label)

            # make root the child of origin
            origin_node.add_child(root_node)

            # now make tree
            tr = dp.Tree(seed_node=origin_node) # will remain a "is_leaf() == True" until we add children

        # simulation starts at root
        else:
            # root node will be the seed_node, and will be untargetable from the get-go
            root_node = dp.Node(taxon=dp.Taxon(label="root"), label="root", edge_length=0.0)
            root_node.state = start_state
            root_node.alive = False
            untargetable_node_set.add("root")

            # we will actually start off with the left and right nodes being targetable
            left_node = dp.Node(taxon=dp.Taxon(label="nd1"), label="nd1", edge_length=0.0)
            left_node.alive = True
            left_node.state = root_node.state # get state from parent (root)
            left_node.annotations.add_bound_attribute("state")
            root_node.add_child(left_node)
            state_representation_dict[left_node.state].add(left_node.label)
            cumulative_node_count += 1

            right_node = dp.Node(taxon=dp.Taxon(label="nd2"), label="nd2", edge_length=0.0)
            right_node.alive = True
            right_node.state = root_node.state # get state from parent (root)
            right_node.annotations.add_bound_attribute("state")
            root_node.add_child(right_node)
            state_representation_dict[right_node.state].add(right_node.label)
            cumulative_node_count += 1

            # now make tree
            tr = dp.Tree(seed_node=root_node) # will remain a "is_leaf() == True" until we add children
            tr.taxon_namespace.add_taxon(left_node)
            tr.taxon_namespace.add_taxon(right_node)

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
            # print("self.slice_t_ends = ")
            # print(self.slice_t_ends)
            # print("time_slice_idx = " + str(time_slice_idx))
            next_max_t = self.slice_t_ends[time_slice_idx]

            # (4) check if new event time cut throw the end of a time slice, provided
            # there is more than 1 time slice, and that a max age was specified
            # through the "age stop condition"
            excess_t = 0.0
            if self.stop == "age" and self.n_time_slices > 1 and latest_t > next_max_t:
                excess_t = latest_t - next_max_t
                latest_t = next_max_t
                time_slice_idx += 1

                # extend all lineages (could keep just the extend in the else-block
                # below) and keep everything up-to-date
                extend_all_living_nodes(t_to_next_event - excess_t)

                if next_max_t == t_stop:
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
                last_chosen_node, cumulative_node_count = self.execute_event(tr.taxon_namespace, macroevol_atomic_param, chosen_node,
                                    state_representation_dict, untargetable_node_set,
                                    cumulative_node_count, last_chosen_node, debug=self.debug)

                # (9) update number of lineages after birth/death event
                #
                # what counts for targetable node here are living terminal nodes
                #
                # for large trees, this scales badly, because we have lots of
                # nodes to visit
                if not macroevol_atomic_param.event == sseobj.MacroevolEvent.ANAGENETIC_TRANSITION:
                    living_nodes = [nd for nd in tr if nd.alive]
                    current_node_target_count = len(living_nodes)

                # (10) check for stop conditions
                #
                # check if all lineages went extinct, stop!
                if current_node_target_count == 0:
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
                    state_representation_dict[last_chosen_node.state].add(last_chosen_node.label) # won't use it again, but to be safe
                    reached_stop_condition = True
                    break

        # (11) got out of while loop because met stop condition
        # 'at' is scoped to simulate_a_tree() function
        if self.stop == "age":
            at = AnnotatedTree(tr, self.events.state_count, start_at_origin=self.with_origin, max_age=self.stop_val, slice_t_ends=self.slice_t_ends, slice_age_ends=self.events.slice_age_ends)
        elif self.stop == "size":
            at = AnnotatedTree(tr, self.events.state_count, start_at_origin=self.with_origin)

        if self.debug:
            # print(at.tree.as_string(schema="nexus", suppress_annotations=True, suppress_internal_taxon_labels=True))
            print(at.tree.as_string(schema="newick", suppress_annotations=False, suppress_internal_taxon_labels=True))
            # print(at.tree)

        return at
    ### END simulation loop ###


    def generate(self):
        # output
        output = list()

        # do something with self.runtime_limit
        start_time = time.time()
        ith_sim = 0
        j = 0
        while len(output) < (self.n_sim * self.n_repl):
            ellapsed_time = pjh.get_ellapsed_time_in_minutes(start_time, time.time())

            if ellapsed_time >= self.runtime_limit:
                break

            # simulate!
            repl_size = 0
            while repl_size < self.n_repl:
                tr = self.simulate(self.start_states[ith_sim], value_idx=ith_sim)

                # check if tr has right specs
                if self.is_tr_ok(tr):
                    output.append(tr)
                    repl_size += 1

                # tree not good, stay in while loop
                else:
                    j += 1
                    continue

            ith_sim += 1

        return output

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        rev_str_list: ty.List[str] = []

        return rev_str_list

    #########################
    # Output health methods #
    #########################
    def is_tr_ok(self, ann_tr):
        """Check that simulated tree meets stop conditions, so that it can be returned

        Args:
            a_tree (AnnotatedTree): A tree with extra annotations.

        Returns:
            bool: If tree meets stop condition value.
        """

        if self.condition_on_survival and ann_tr.tree_died:
            return False

        if self.stop == "size":
            if self.condition_on_speciation and len(ann_tr.root_node.child_nodes()) == 0:
                return False

            if self.condition_on_survival and (ann_tr.n_extant_obs_nodes + ann_tr.n_sa_obs_nodes) != self.stop_val:
                return False

            # if conditions are not specified, returns all trees
            return True

            # if (ann_tr.n_extant_obs_nodes + ann_tr.n_sa_obs_nodes) in (self.stop_val, 0):
            #     return True
            # else:
            #     return False

        elif self.stop == "age":
            if ann_tr.with_origin:
                # tree is ok if either it reached the max age by epsilon, or if its age is smaller than max age
                if abs(self.stop_val - ann_tr.origin_age) <= self.epsilon or ann_tr.origin_age < self.stop_val:
                    return True
                elif ann_tr.origin_age > self.stop_val:
                    return False

            else:
                if self.stop_val >= ann_tr.root_age:
                    return True
                else:
                    return False