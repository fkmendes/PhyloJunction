import os
import copy
import typing as ty
import matplotlib.pyplot as plt  # type: ignore

# pj imports
from phylojunction.data.attribute_transition import AttributeTransition
import phylojunction.data.tree as pjtr
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.evol_event as pjev
import phylojunction.functionality.biogeo as pjbio
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class StochMap(pjev.EvolRelevantEvent):
    """Parent class for a single stochastic map.
    
    Parameters:
        from_state (int): Integer representation of source state (i.e.,
            source biogeographic range).
        from_state_bit_patt (str): Bit pattern representation of
            source state (i.e., source biogeographic range).
        to_state (int): Integer representation of target state (i.e.,
            target biogeographic range).
        to_state_bit_pattern (str): Bit pattern representation of
            target state (i.e., target biogeographic range).
        to_state2 (int): Integer representation of second target state
            if stochastic map is cladogenetic (i.e., second target
            biogeographic range).
        to_state2_bit_pattern (str): Bit pattern representation of
            second target state if stochastic map is cladogenetic
            (i.e., second target biogeographic range).
        focal_node_name (str): Name of the node subtending the branch
            along which the stochastic map was placed.
        parent_node_name (str): Name of the node parent to the focal
            node.
        child_node_name (str): Name of the first child node of the
            focal node.
        child2_node_name (str): Name of the second child node of the
            focal node.
        map_type (str): Class of stochastic map (e.g., 'expansion',
        'contraction', 'range_split', 'no_change'). Defaults to None.
    """

    from_state: int
    from_state_bit_patt: str
    to_state: int
    to_state_bit_patt: str
    to_state2: int
    to_state2_bit_patt: str
    focal_node_name: str
    parent_node_name: str
    child_node_name: str
    child2_node_name: str
    map_type: str

    def __init__(self,
                 n_chars: int,
                 age: ty.Optional[float],
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 focal_node_name: str,
                 parent_node_name: str,
                 child_node_name: str,
                 map_type: ty.Optional[str] = None,
                 time: ty.Optional[int] = None,
                 to_state2: ty.Optional[int] = None,
                 to_state2_bit_patt: ty.Optional[str] = None,
                 child2_node_name: ty.Optional[str] = None) -> None:

        # initializing EvolRelevantEvent members
        super().__init__(n_chars, age, time=time)

        self.map_type = map_type
        self.from_state = from_state
        self.from_state_bit_patt = from_state_bit_patt
        self.to_state = to_state
        self.to_state_bit_patt = to_state_bit_patt
        self.focal_node_name = focal_node_name
        self.parent_node_name = parent_node_name
        self.child_node_name = child_node_name

        # if cladogenetic range split
        self.child2_node_name = child2_node_name
        self.to_state2 = to_state2
        self.to_state2_bit_patt = to_state2_bit_patt


class NoChange(StochMap):
    """Stochastic map class for no-change.

    When nothing happens along a branch, we want to still record
    this for updating nodes in an AnnotatedTree object.

    Parameters:
        str_representation (str): String representation of stochastic
            map for printing.
        state (int): Integer representation of node's state.
    """

    state: int
    size_of_unchanging_range: int
    str_representation: str

    def __init__(self,
               n_regions: int,
               state: int,
               state_bit_patt: str,
               focal_node_name: str,
               parent_node_name: str,
               child_node_name: str) -> None:

        super().__init__(n_regions,
                         None,
                         state,
                         state_bit_patt,
                         state,
                         state_bit_patt,
                         focal_node_name,
                         parent_node_name,
                         child_node_name,
                         map_type="no_change")

        self.size_of_unchanging_range = sum(int(i) for i in state_bit_patt)
        self._init_str_representation()

    def _init_str_representation(self) -> None:
        self.str_representation = \
            "No change" \
            + "\n  Node subtending no changes: " + self.focal_node_name \
            + "\n    Unchanging state " + str(self.from_state) + ", bits \'" \
            + self.from_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_unchanging_range)

    def __str__(self) -> str:
        return self.str_representation


class RangeExpansion(StochMap):
    """Stochastic map class for (anagenetic) dispersal.
    
    When a lineage disperses, its range is incremented by one region.

    Parameters:
        str_representation (str): String representation of stochastic
            map for printing.
        size_of_expanding_range (int): Size of starting range.
        size_of_final_range (int): Size of resulting range after
            dispersal.
        region_gained_idx (int): Index of region that was added to the
            starting range (index is the position of the region in the
            bit pattern).
        over_barrier (bool, optional): Flag specifying if dispersal was over a
            barrier. Note that even though the dispersal may have been
            over a barrier, it may not matter for a focal range-split
            (speciation) event.
        split_relevant (bool, optional): Flag specified if dispersal
            matters for range split happening in the future.
        erased_by_extinction (bool, optional): Flag specifying if
            the expanded range contracted in the reverse direction
            later on, due to extinction.
    """

    str_representation: str
    size_of_expanding_range: int
    size_of_final_range: int
    from_region_idx: ty.Optional[int]
    region_gained_idx: int

    # additional annotation for testing hypotheses
    over_barrier: ty.Optional[bool]
    split_relevant: ty.Optional[bool]
    erased_by_extinction: ty.Optional[bool]

    def __init__(self,
                 region_gained_idx: int,
                 n_regions: int,
                 age: int,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 focal_node_name: str,
                 parent_node_name: str,
                 child_node_name: str,
                 from_region_idx: ty.Optional[int] = None,
                 time: ty.Optional[float] = None):
        
        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         focal_node_name,
                         parent_node_name,
                         child_node_name,
                         map_type="expansion",
                         time=time)

        self.from_region_idx = from_region_idx
        self.region_gained_idx = region_gained_idx
        self.size_of_expanding_range = \
            sum(int(i) for i in from_state_bit_patt)
        self.size_of_final_range = \
            sum(int(i) for i in to_state_bit_patt)

        # additional annotation for testing hypotheses
        self.over_barrier = None
        self.split_relevant = None
        self.erased_by_extinction = None

    def _init_str_representation(self) -> None:
        self.str_representation = \
            "Range expansion / Dispersal (at age = " + str(self.age) + ")" \
            + "\n  Node subtending range expansion: " + self.focal_node_name \
            + "\n    From state " + str(self.from_state) + ", bits \'" \
            + self.from_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_expanding_range) \
            + "\n    To state " + str(self.to_state) + ", bits \'" \
            + self.to_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_final_range) \
            + "\n    Dispersal from region (index): " + str(self.from_region_idx) \
            + "\n    Region gained (index): " + str(self.region_gained_idx) \
            + "\n  Dispersal over barrier: " + str(self.over_barrier) \
            + ("\n  Dispersal relevant to split: ") + str(self.split_relevant)

    # setter
    def gained_region_is_lost_in_future(self):
        self.range_expansion_erased_by_extinction = True

    def __str__(self) -> str:
        self._init_str_representation()

        return self.str_representation


class RangeContraction(StochMap):
    """Stochastic map class for local (anagenetic) extinction.
    
    When a lineage goes locally extinct, its range contracts by
    one region.

    Attributes:
        str_representation (str): String representation of stochastic
            map for printing.
        size_of_contracting_range (int): Size of starting range.
        size_of_final_range (int): Size of resulting range after
            dispersal.
        region_lost_idx (int): Index of region that was added to the
            starting range (index is the position of the region in the
            bit pattern).
    """

    str_representation: str
    size_of_contracting_range: int
    size_of_final_range: int
    region_lost_idx: int

    def __init__(self,
                 region_lost_idx: int,
                 n_regions: int,
                 age: float,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 focal_node_name: str,
                 parent_node_name: str,
                 child_node_name: str,
                 time: ty.Optional[float] = None) -> None:
        
        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         focal_node_name,
                         parent_node_name,
                         child_node_name,
                         map_type="contraction",
                         time=time)

        self.region_lost_idx = region_lost_idx
        self.size_of_contracting_range = \
            sum(int(i) for i in from_state_bit_patt)
        self.size_of_final_range = \
            sum(int(i) for i in to_state_bit_patt)

        self._init_str_representation()

    def _init_str_representation(self) -> None:
        self.str_representation = \
            "Range contraction / Local extinction (at age = " + str(self.age) + ")" \
            + "\n  Node subtending range contraction: " + self.focal_node_name \
            + "\n    From state " + str(self.from_state) + ", bits \'" \
            + self.from_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_contracting_range) \
            + "\n    To state " + str(self.to_state) + ", bits \'" \
            + self.to_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_final_range) \
            + "\n    Region lost (index): " + str(self.region_lost_idx)

    def __str__(self) -> str:
        return self.str_representation


class RangeSplitOrBirth(StochMap):
    """Stochastic map class for cladogenetic event.
    
    Range either splits between children, e.g., ABC -> AB, C 
    (upon b/w-region speciation)
    
    or

    One-sized range is born and inherited by only one child,
    e.g., ABC -> ABC, A (upon within-region speciation).

    Attributes:
        str_representation (str): String representation of stochastic
            map for printing.
        range_split (bool): Flag specifying if range split.
        size_of_splitting_node_range (int): Size of parent's starting
            range, before split.
        size_of_child1_range (int): Size of first child range.
        size_of_child2_range (int): Size of second child range.
    """

    str_representation: str
    
    # this is when, under GeoSSE, at an internal node
    # we have speciation within a region, e.g., AB -> AB,A
    range_split: bool  # used for printing only
    
    size_of_splitting_node_range: int
    size_of_child1_range: int
    size_of_child2_range: int

    def __init__(self,
                 n_regions: int,
                 age: float,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 focal_node_name: str,
                 parent_node_name: str,
                 child_node_name: str,
                 child2_node_name: str,
                 to_state2: ty.Optional[int] = None,
                 to_state2_bit_patt: ty.Optional[str] = None,
                 time: ty.Optional[float] = None) -> None:

        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         focal_node_name,
                         parent_node_name,
                         child_node_name,
                         child2_node_name=child2_node_name,
                         to_state2=to_state2,
                         to_state2_bit_patt=to_state2_bit_patt,
                         map_type="range_split",
                         time=time)
        
        self.size_of_splitting_node_range = \
            sum(int(i) for i in from_state_bit_patt)
        self.size_of_child1_range = \
            sum(int(i) for i in to_state_bit_patt)
        self.range_split = False
        
        self._init_str_representation()

    def add_child2_info(self,
                        to_state2: int,
                        to_state2_bit_patt: str) -> None:

        self.to_state2 = to_state2
        self.to_state2_bit_patt = to_state2_bit_patt
        self.size_of_child2_range = \
            sum(int(i) for i in to_state2_bit_patt)
        self._update_str_representation()
        self.range_split = True

    def _init_str_representation(self) -> None:
        self.str_representation = \
            "Range split (at age = " + str(self.age) + ")" \
            + "\n  Splitting node: " + self.focal_node_name \
            + " (state " + str(self.from_state) + ", bits \'" \
            + self.from_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_splitting_node_range) + ")" \
            + "\n  Child 1 node: " + self.child_node_name \
            + " (state " + str(self.to_state) + ", bits \'" \
            + self.to_state_bit_patt + "\', " \
            + "range size " + str(self.size_of_child1_range) + ")"
    
    def _update_str_representation(self) -> None:
        self.str_representation += \
            "\n  Child 2 node: " + self.child2_node_name \
            + " (state " + str(self.to_state2) + ", bits \'" \
            + self.to_state2_bit_patt + "\', " \
            + "range size " + str(self.size_of_child2_range) + ")"
        
    def __str__(self) -> str:
        if self.range_split:
            return self.str_representation
        
        else:
            return self.str_representation \
                + "\n  Child 2 node: same as parent"


class StochMapsOnTree():
    """Collection of stochastic maps on a single tree.

    This collection of stochastic maps consists of all stochastic
    maps logged during a single iteration.

    Parameters:
        ann_tr (AnnotatedTree): AnnotatedTree object on which maps
            are placed. There should be one different tree for each
            StochMapsOnTree object.
        no_change_stoch_maps_dict (dict): Dictionary with node names
            as keys and a single StochMap object as values.
        anagenetic_stoch_maps_dict (dict): Dictionary with node names
            as keys and a list of StochMap objects as values.
        cladogenetic_stoch_maps_dict (dict): Dictionary with node_names
            as keys and a single StochMap object as values.
        identical_cladogenetic_stoch_maps_dict (dict):
        iteration_idx (int): 
        state2bit_lookup (State2BitLookup): Object of class for
            converting between integer and bit pattern state
            representations.
        n_regions (int): Total number of atomic regions in the system.
        maps_counts_dict (dict): Dictionary with two-sized string
            tuples as keys (string are bit patterns for source and
            target ranges), and counts as values. It keeps track of
            how many stochastic maps of different kinds there are.
        n_no_changes (int): Total count of no change maps.
        n_total_maps (int): Total count of stochastic maps.
        n_anagenetic_maps (int): Count of stochastic maps consisting
            of anagenetic range changes.
        n_range_exp_maps (int): Count of range expansion stochastic
            maps.
        n_range_cont_maps (int): Count of range contraction stochastic
            maps.
        pctg_anagenetic_maps (float): Percentage of anagenetic
            stochastic maps.
        forbidden_higher_order_maps_counts_dict (dict): Dictionary with
            two-sized string tuples as keys (string are bit patterns
            for source and target ranges), and counts as values. This
            dictionary keeps track of how many forbidden higher-order
            (more than one region is gained/lost simultaneously) maps
            there are.
        n_forbidden_higher_order_anagenetic_maps (int): Count of
            forbidden higher-order anagenetic stochastic maps. Higher-
            order anagenetic changes happen when one or more regions
            are gained or lost simultaneously.
        pctg_forbidden_higher_order_anagenetic_maps (float): Percentage
            of forbidden higher-order anagenetic stochastic maps.
        n_different_clado_maps (int): Count of different cladogenetic
            stochastic maps.
        n_identical_clado_maps (int): Count of identical cladogenetic
            stochastic maps.
        str_representation (str): String representation of collection
            of stochastic maps for printing.
    """

    # class members
    ann_tr: pjtr.AnnotatedTree
    no_change_stoch_maps_dict: \
        ty.Dict[str, StochMap]  # node name as key
    anag_stoch_maps_dict: \
        ty.Dict[str, ty.List[StochMap]]  # node name as key
    clado_stoch_maps_dict: \
        ty.Dict[str, StochMap]  # node name as key
    identical_clado_stoch_maps_dict: \
        ty.Dict[str, ty.List[str]]  # iteration idx as key, node idxs as values
    iteration_idx: int
    str_representation: str
    state2bit_lookup: pjbio.State2BitLookup
    n_chars: int

    # for summary and printing
    maps_counts_dict: ty.Dict[ty.Tuple[str], int]
    n_no_change_maps: int
    n_total_maps: int
    n_anag_maps: int
    n_range_exp_maps: int
    n_range_cont_maps: int
    pctg_anagenetic_maps: float
    forbidden_higher_order_maps_counts_dict: \
        ty.Dict[ty.Tuple[str], int]
    n_forbidden_higher_order_anag_maps: int
    pctg_forbidden_higher_order_anag_maps: float
    n_different_clado_maps: int
    n_identical_clado_maps: int

    # iteration
    iteration_idx: int
    
    def __init__(self,
                 iteration_idx: int,
                 ann_tr: pjtr.AnnotatedTree,
                 state2bit_lookup: pjbio.State2BitLookup,
                 node_attr_file_path: str,
                 stoch_map_attr_name: str) -> None:
        
        self.ann_tr = ann_tr
        self.iteration_idx = iteration_idx
        self.state2bit_lookup = state2bit_lookup
        self.n_regions = self.state2bit_lookup.n_char

        self.no_change_stoch_maps_dict = dict()
        self.anag_stoch_maps_dict = dict()
        self.clado_stoch_maps_dict = dict()
        self.identical_clado_stoch_maps_dict = dict()
        
        # summary / counts
        self.maps_counts_dict = dict()
        self.forbidden_higher_order_maps_counts_dict = dict()
        self.n_no_change_maps = 0
        self.n_total_maps = 0
        self.n_anag_maps = 0
        self.n_range_exp_maps = 0
        self.n_range_cont_maps = 0
        self.pctg_anagenetic_maps = 0.0
        self.n_forbidden_higher_order_anag_maps = 0
        self.pctg_forbidden_higher_order_anag_maps = 0.0
        self.n_different_clado_maps = 0
        self.pctg_diff_clado_maps = 0.0
        self.n_identical_clado_maps = 0
        self.pctg_ident_clado_maps = 0.0

        # reading terminal node attribute values and updating
        # AnnotatedTree object
        if node_attr_file_path and stoch_map_attr_name:
            pjr.read_node_attr_update_tree(node_attr_file_path,
                                           stoch_map_attr_name,
                                           int,
                                           self.ann_tr)

    def add_map(self, stoch_map_str_list: ty.List[str]) -> None:
        """Add new stochastic map to bookeeping class members.

        This method reads a list of strings carrying information about
        a stochastic map, and depending on what those strings are, the
        appropriate stochastic map objects are created and stored.

        Side-effects of this method populate:
            (i)    self.n_total_maps
            (ii)   self.n_no_change_maps
            (iii)  self.n_anag_maps
            (iv)   self.n_identical_clado_maps
            (v)    self.n_different_clado_maps
            (vi)   self.n_range_exp_maps
            (vii)  self.n_cont_exp_maps
            (viii) self.n_forbidden_higher_order_anag_maps
            (ix)   self.no_change_stoch_maps_dict
            (x)    self.anag_stoch_maps_dict
            (xi)   self.clado_stoch_maps_dict
            (xii)  self.identical_clado_stoch_maps_dict
            (xiii) self.forbidden_higher_order_maps_counts_dict
        """

        def _which_region_changed(from_state: int, to_state: int) \
                -> ty.Tuple[bool, int]:
            pass

        nd_name, parent_nd_name, ch1_nd_name, ch2_nd_name = \
            str(), str(), str(), str()

        # print("n_forbidden before add_map=", self.n_forbidden_higher_order_anag_maps)
        nd_idx, br_start_age, br_end_age, from_state, to_state, \
        event_age, event_type, parent_idx, ch1_idx, ch2_idx \
            = stoch_map_str_list

        # populating node names
        nd_name = self.ann_tr.tree.find_node(
            filter_fn=lambda nd: nd.index == nd_idx).label

        parent_nd = self.ann_tr.tree.find_node(
            filter_fn=lambda nd: nd.index == parent_idx)
        if parent_nd:
            parent_nd_name = parent_nd.label

        ch1_nd = self.ann_tr.tree.find_node(
            filter_fn=lambda nd: nd.index == ch1_idx)
        if ch1_nd:
            ch1_nd_name = ch1_nd.label

        ch2_nd = self.ann_tr.tree.find_node(
            filter_fn=lambda nd: nd.index == ch2_idx)
        if ch2_nd:
            ch2_nd_name = ch2_nd.label

        # debugging
        # print("parent_nd_name", parent_nd_name, "ch1_nd_name", ch1_nd_name, "ch2_nd_name", ch2_nd_name)

        # bit patterns
        from_state_bit_patt = self.state2bit_lookup.get_bit(int(from_state))
        to_state_bit_patt = self.state2bit_lookup.get_bit(int(to_state))

        try:
            event_time = self.ann_tr.seed_age - float(event_age)

        except ValueError:
            event_time = None

        if event_type == "no_change":
            self.n_no_change_maps += 1
            stoch_map = NoChange(self.n_regions,
                                 int(from_state),
                                 from_state_bit_patt,
                                 nd_name,
                                 parent_nd_name,
                                 ch1_nd_name)

            self.no_change_stoch_maps_dict[nd_name] = stoch_map

        elif event_type == "anagenetic":
            self.n_total_maps += 1
            seen_different = False
            is_higher_order = False
            is_expansion = True
            stoch_map: StochMap = None

            # parsing bit by bit
            for idx, b_from, in enumerate(from_state_bit_patt):
                b_to = to_state_bit_patt[idx]

                if seen_different and b_from != b_to:
                    is_higher_order = True
                    break
            
                # determining if expansion or contraction
                elif not seen_different and b_from != b_to:
                    seen_different = True
                    
                    if b_from == "0" and b_to == "1":
                        continue
                    
                    elif b_from == "1" and b_to == "0":
                        is_expansion = False

            #####################################
            # Getting counts for printing later #
            #####################################

            if is_higher_order:
                # populate forbidden maps counts dict
                if (from_state_bit_patt, to_state_bit_patt) \
                        not in self.forbidden_higher_order_maps_counts_dict:
                    self.forbidden_higher_order_maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] = 1
                
                else:
                    self.forbidden_higher_order_maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] += 1

                self.n_forbidden_higher_order_anag_maps += 1
            
            else:
                # populate maps counts dict
                if (from_state_bit_patt, to_state_bit_patt) \
                        not in self.maps_counts_dict:
                    self.maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] = 1
                
                else:
                    self.maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] += 1
                    
                self.n_anag_maps += 1

                #####################
                # Initializing maps #
                #####################
                if is_expansion:
                    # find out which region was gained
                    region_gained_idx = -1
                    count_gained = count_lost = 0
                    for i, b in enumerate(from_state_bit_patt):
                        if b == "0" and to_state_bit_patt[i] == "1":
                            if count_gained > 0:
                                exit("I found two bits that were different. This should not have happened.")
                            region_gained_idx = i
                            count_gained += 1

                    # now build stoch map!
                    stoch_map = RangeExpansion(region_gained_idx,
                                               self.n_regions,
                                               float(event_age),
                                               int(from_state),
                                               from_state_bit_patt,
                                               int(to_state),
                                               to_state_bit_patt,
                                               nd_name,
                                               parent_nd_name,
                                               ch1_nd_name,
                                               time=event_time)
                    
                    self.n_range_exp_maps += 1

                else:
                    # find out which region was lost
                    region_lost_idx = -1
                    count_lost = 0
                    for i, b in enumerate(from_state_bit_patt):
                        if b == "1" and to_state_bit_patt[i] == "0":
                            if count_lost > 0:
                                exit("I found two bits that were different. This should not have happened.")
                            region_lost_idx = i
                            count_lost += 1

                    # now build stoch map!
                    stoch_map = RangeContraction(region_lost_idx,
                                                 self.n_regions,
                                                 event_age,
                                                 int(from_state),
                                                 from_state_bit_patt,
                                                 int(to_state),
                                                 to_state_bit_patt,
                                                 nd_name,
                                                 parent_nd_name,
                                                 ch1_nd_name,
                                                 time=event_time)
                    
                    self.n_range_cont_maps += 1
                        
                # TODO: later, replace idx with node name
                if nd_name not in self.anag_stoch_maps_dict:
                    self.anag_stoch_maps_dict[nd_name] = [stoch_map]

                else:
                    self.anag_stoch_maps_dict[nd_name].append(stoch_map)

        elif event_type == "cladogenetic":
            # first time we add cladogenetic map
            if nd_name not in self.clado_stoch_maps_dict:
                self.n_total_maps += 1
                self.n_different_clado_maps += 1

                stoch_map = RangeSplitOrBirth(self.n_regions,
                                              event_age,
                                              int(from_state),
                                              from_state_bit_patt,
                                              int(to_state),
                                              to_state_bit_patt,
                                              nd_name,
                                              parent_nd_name,
                                              ch1_nd_name,
                                              ch2_nd_name,
                                              time=event_time)
                
                self.clado_stoch_maps_dict[nd_name] = stoch_map

            # second cladogenetic map entry (child 2!)
            else:
                # just counting here
                if self.clado_stoch_maps_dict[nd_name].from_state \
                    == int(from_state) and \
                    self.clado_stoch_maps_dict[nd_name].to_state \
                        == int(to_state):
                    
                    if self.iteration_idx \
                            not in self.identical_clado_stoch_maps_dict:
                        self.identical_clado_stoch_maps_dict \
                            [self.iteration_idx] = [nd_name]
                    
                    else:
                        self.identical_clado_stoch_maps_dict \
                            [self.iteration_idx].append(nd_name)
                    
                    self.n_identical_clado_maps += 1
                    self.n_different_clado_maps -= 1

                # updating the cladogenetic map with info from the second child
                self.clado_stoch_maps_dict[nd_name] \
                    .add_child2_info(int(to_state), to_state_bit_patt)
                
            # debugging
            # for nd_name, rsom in self.cladogenetic_stoch_maps_dict.items():
            #     print("Iteration ", str(self.iteration_idx) + ",", "node", nd_name, "\n", rsom)

    def get_clado_distribution(self) -> ty.Dict[int, int]:
        widespread_range_sizes_dict = dict()
        for nd_name, clado_stoch_map \
                in self.clado_stoch_maps_dict.items():
            
            # print(clado_stoch_map.from_state, clado_stoch_map.from_state_bit_patt)
            
            widespread_range_size = \
                sum(int(i) for i in list(clado_stoch_map.from_state_bit_patt))
            
            if not widespread_range_size in widespread_range_sizes_dict:
                widespread_range_sizes_dict[widespread_range_size] = 1

            else:
                widespread_range_sizes_dict[widespread_range_size] += 1

        return widespread_range_sizes_dict

    def update_tree_attributes(self,
                               stoch_map_attr_name: str,
                               debug: bool = False) -> None:
        """Go through stochastic maps and update tree.

        This function is called outside of StochMapsOnTree.
        It goes through the two stochastic maps dictionaries,
        and updates the AnnotatedTree member with respect to
        its (i) node's attributes, and (ii) node_attr_dict member

        Args:
            stoch_map_attr_name (str): Name of the node attribute
                the stochastic maps carry a value for (e.g., 'state')

        Returns:
            None, this function has a side-effect
        """

        # debugging
        # print("Iteration", self.iteration_idx, ": updating tree attributes!")
        # for nd_name, sm in self.cladogenetic_stoch_maps_dict.items():
        #     print("clado", nd_name, sm.from_state, sm.age, sm.time)
        # for nd_name, sm_list in self.anagenetic_stoch_maps_dict.items():
        #     for sm in sm_list:
        #         print("anag", nd_name, sm.from_state, sm.to_state, sm.age, sm.time)

        ########################################
        # Updating the states nodes subtending #
        # lineages with no changes             #
        ########################################

        no_changes_nd_idx_set: ty.Set[str] = set()
        for nd_name, smap in self.no_change_stoch_maps_dict.items():
            nd = self.ann_tr.tree.find_node(filter_fn = \
                                            lambda nd: nd.label == nd_name)

            setattr(nd, stoch_map_attr_name, smap.from_state)
            no_changes_nd_idx_set.add(nd_name)

        ################################################
        # Creating and populating clado_at_dict member #
        # of AnnotatedTree object                      #
        #
        # This step must be carried out before going   #
        # through anagenetic stochastic maps because   #
        # cladogenetic state transition events         #
        # contribute the first anagenetic transition   #
        # events                                       #
        ################################################

        for nd_name, smap in self.clado_stoch_maps_dict.items():
            nd = self.ann_tr.tree.find_node(filter_fn = \
                                                lambda nd: nd.label == nd_name)

            # we set the state of the node in the AnnotatedTree
            setattr(nd, stoch_map_attr_name, smap.from_state)
            if debug:
                print("nd_name", nd_name, stoch_map_attr_name, smap.from_state)

            clado_state_trans = AttributeTransition(
                "state",
                nd_name,
                smap.time,
                smap.from_state,
                smap.to_state,
                smap.to_state2,
                at_speciation=True
            )

            # adding cladogenetic state transition
            self.ann_tr.clado_at_dict[nd.label] = clado_state_trans

            # for child_nd in nd.child_node_iter():
            #     smap_to_state = -1
            #
            #     # debugging
            #     # print("smap.child_node_name", smap.child_node_name, "smap.child2_node_name", smap.child2_node_name)
            #     # print("child_nd.label", child_nd.label)
            #
            #     if child_nd.label == smap.child_node_name:
            #         smap_to_state = smap.to_state
            #
            #     elif child_nd.label == smap.child2_node_name:
            #         smap_to_state = smap.to_state2
            #
            #     # if child_nd.label == "nd4":
            #     #     print("Inside clado_stoch_maps:")
            #     #     print("smap.time", smap.time)
            #     #     print("smap.from_state", smap.from_state, "smap_to_state", smap.to_state)
            #
            #     ana_state_trans = AttributeTransition(
            #         "state",
            #         child_nd.label,
            #         smap.time,
            #         smap.from_state,
            #         smap_to_state,
            #         at_speciation=True
            #     )
            #
            #     # adding anagenetic state transition (used for plotting)
            #     self.ann_tr.at_dict[child_nd.label] = [ana_state_trans]

        ##############################################
        # Creating and populating at_dict member of  #
        # of AnnotatedTree object                    #
        ##############################################

        for nd_name, smap_list in self.anag_stoch_maps_dict.items():
            nd = self.ann_tr.tree.find_node(filter_fn = \
                                                lambda nd: nd.label == nd_name)

            for smap in smap_list:
                # we set the state of the node in the AnnotatedTree
                setattr(nd, stoch_map_attr_name, smap.to_state)
                if debug:
                    print("nd_name", nd_name, stoch_map_attr_name, smap.to_state)

                ana_state_trans = AttributeTransition(
                    "state",
                    nd_name,
                    smap.time,
                    smap.from_state,
                    smap.to_state
                )

                # if nd_name == "nd4":
                #     print("Inside anag_stoch_maps:")
                #     print("smap.time", smap.time)
                #     print("smap.from_state", smap.from_state, "smap_to_state", smap.to_state)

                # adding anagenetic state transition (used for plotting)
                if nd.label in self.ann_tr.at_dict:
                    self.ann_tr.at_dict[nd_name].append(ana_state_trans)

                else:
                    self.ann_tr.at_dict[nd_name] = [ana_state_trans]

        ############################################################
        # Now we update the node_attr_dict member of AnnotatedTree #
        ############################################################

        self.ann_tr.populate_nd_attr_dict([stoch_map_attr_name])

    # internal
    def _preorder_traverse_tree_annotating_maps(self, ann_tr):
        pass

    def _init_str_representation(self) -> None:
        self._stoch_maps_summary()

        self.str_representation = "Stochastic maps in MCMC iteration " \
            + str(self.iteration_idx) + ":"
        
        self.str_representation += \
            "\n    Number of anagenetic changes = " \
            + str(self.n_anag_maps) \
            + " (" + str(self.pctg_anagenetic_maps) + "%)" \
            + "\n        Range contractions = " + str(self.n_range_cont_maps) \
            + "\n        Range expansions = " + str(self.n_range_exp_maps) \
            + "\n    Number of higher-order anagenetic changes = " \
            + str(self.n_forbidden_higher_order_anag_maps) \
            + " (" + str(self.pctg_forbidden_higher_order_anag_maps) + "%)" \
            + "\n    Number of cladogenetic changes = " \
            + str(self.n_different_clado_maps) \
            + " (" + str(self.pctg_diff_clado_maps) + "%)" \
            + "\n    Number of identical cladogenetic changes = " \
            + str(self.n_identical_clado_maps) \
            + " (" + str(self.pctg_ident_clado_maps) + "%)" \

    # updates
    #     self.pctg_anagenetic_maps
    #     self.pctg_forbidden_higher_order_anagenetic_maps
    #     self.pctg_diff_clado_maps
    #     self.pctg_ident_clado_maps
    def _stoch_maps_summary(self) -> None:
        self.pctg_anagenetic_maps = round(100.0 * (self.n_anag_maps \
            / self.n_total_maps), 2)
        
        self.pctg_forbidden_higher_order_anag_maps = \
            round(100.0 * (self.n_forbidden_higher_order_anag_maps \
            / self.n_total_maps), 2)
        
        self.pctg_diff_clado_maps = round(100.0 * (self.n_different_clado_maps \
            / self.n_total_maps), 2)
        
        self.pctg_ident_clado_maps = round(100.0 * (self.n_identical_clado_maps \
            / self.n_total_maps), 2)

    def __str__(self) -> str:
        self._init_str_representation()

        return self.str_representation


class StochMapsOnTreeCollection():
    """Collection of StochMapsOnTree objects.

    This class reads through a table with (normally multiple)
    iterations of logged stochastic maps, and maps them to
    one or more annotated trees -- storing them as objects of
    StochMapsOnTree. These are in turn stored in class members.

    Attributes:
        n_char (int): Number of characters.
        stoch_maps_tree_dict (dict): Dictionary with iteration
            indices as keys, and StochMapsOnTree objects as values.
        n_stoch_map_iterations (int): Total number of MCMC iterations
            where stochastic maps were logged.
        sorted_it_idxs (int): This is the sorted indices of
            iterations for which stochastic maps were logged. Not
            all iterations have maps, which is why we keep track of
            those who do.
        str_representation (str): String representation of class
            instance for printing.
    """

    n_char: int
    stoch_maps_tree_dict: \
        ty.Dict[int, StochMapsOnTree]
    n_stoch_map_iterations: int
    sorted_it_idxs: ty.List[int]
    str_representation: str

    def __init__(self,
                 stoch_maps_file_path: str,
                 ann_trs: ty.List[pjtr.AnnotatedTree],
                 state2bit_lookup: pjbio.State2BitLookup,
                 stoch_map_attr_name: str = "",
                 node_states_file_path: str = "") -> None:

        self.stoch_maps_tree_dict = dict()
        self.sorted_it_idxs = list()

        # side-effects:
        # initializes
        # (i)   n_char
        # (i)   stoch_maps_tree_list
        # (ii)  sorted_iteration_idxs (also sorts it)
        # (iii) n_stoch_map_iterations
        #
        # updates members of AnnotatedTree objects passed
        # in 'ann_trs'
        self._read_stoch_maps_file(stoch_maps_file_path,
                                   node_states_file_path,
                                   stoch_map_attr_name,
                                   ann_trs,
                                   state2bit_lookup)
        
        # looks for errors with cladogenetic maps
        self._report_clado_map_issue_str()

        # self.str_representation = ""
        self._init_str_representation()


    # internal
    def _read_stoch_maps_file(self,
                              stoch_maps_file_path: str,
                              node_states_file_path: str,
                              stoch_map_attr_name: str,
                              ann_trs: ty.List[pjtr.AnnotatedTree],
                              state2bit_lookup: \
                              pjbio.State2BitLookup) -> None:
        """Read text file with stochastic maps and populate members.
        
        Side-effects of this method populate:
            (i)   self.n_char
            (ii)  self.stoch_maps_tree_dict
            (iii) self.n_stoch_map_iterations
            (iv)  self.sorted_it_idxs

        Args:
            stoch_maps_file_path (str): Path to the .tsv file
                containing stochastic maps, in the format outputted by
                RevBayes.
            node_states_file_path (str): Path to the .tsv file
                containing the states of all terminal nodes. Nodes
                must be named as 'nd' + index, where index is provided
                as metadata for each node.
            stoch_map_attr_name (str): Name of the attribute (e.g.,
                'state') whose transitions are being stochastically
                mapped.
            ann_trs (AnnotatedTree): List of AnnotatedTree objects. Can
                have a single tree inside, in which case it is deep-
                cloned should there be multiple iterations of
                stochastic maps.
            state2bit_lookup (State2BitLookup): State2BitLookup object
                for......
        """

        self.n_char = state2bit_lookup.n_char

        if not os.path.isfile(stoch_maps_file_path):
            exit("ERROR: Could not find " + stoch_maps_file_path + ". Exiting.")
        
        with open(stoch_maps_file_path, "r") as infile:
            infile.readline() # skip header
            # load all maps in memory
            stoch_map_on_tree_str_list = infile.read().splitlines()

        # how many MCMC iterations we have maps for
        self.n_stoch_map_iterations = \
            len(set(line.split("\t")[0] \
                    for line in stoch_map_on_tree_str_list))
            
        # either we have a single tree (n_samples = 1)
        # and multiple maps (stochastic maps from
        # a single analysis at different MCMC generations),
        # in which that tree is deep copied to match the number
        # of maps
        #
        # or we have multiple trees (n_samples > 1 and/or
        # n_repl > 1) and the number of MCMC stochastic maps
        # must match the number of trees
        n_trees: int = len(ann_trs)
        if n_trees > 1 and n_trees != self.n_stoch_map_iterations:
            raise ec.ObjInitIncorrectDimensionError(
                "StochMapCollection",
                "ann_trs",
                n_trees,
                self.n_stoch_map_iterations)

        # deep-cloning AnnotatedTree objects so that each can get
        # a different set of stochastic maps (one set per iteration)
        if n_trees == 1:
            n_tree_copies_to_make = self.n_stoch_map_iterations - 1

            for i in range(n_tree_copies_to_make, 0, -1):
                ann_trs.append(copy.deepcopy(ann_trs[0]))

        # debugging
        # print("n_iterations", self.n_stoch_map_iterations)

        # iterating over MCMC iterations when stochastic maps were logged
        n_different_iterations: int = 0
        for line in stoch_map_on_tree_str_list:
            stoch_map_on_tree_str_list = line.split("\t")
            it_idx = int(stoch_map_on_tree_str_list[0])

            if it_idx not in self.sorted_it_idxs:
                self.sorted_it_idxs.append(it_idx)

            if it_idx not in self.stoch_maps_tree_dict:
                self.stoch_maps_tree_dict[it_idx] \
                    = StochMapsOnTree(it_idx,
                                      ann_trs[n_different_iterations],
                                      state2bit_lookup,
                                      node_states_file_path,
                                      stoch_map_attr_name)

                n_different_iterations += 1

            # we append this map to the members of the
            # already initialized smot
            self.stoch_maps_tree_dict[it_idx] \
                .add_map(stoch_map_on_tree_str_list[1:])
                
                # print("just added map", iteration_idx, "n_forbidden =", self.stoch_maps_tree_dict[iteration_idx].n_forbidden_higher_order_anagenetic_maps)
        
        sorted(self.sorted_it_idxs)

        #############################################
        # Make every SMOT update its tree's members #
        # according to the stochastic maps of their #
        # assigned iteration                        #
        #############################################
        for it_idx, smot in self.stoch_maps_tree_dict.items():
            # sets critical properties in stoch maps
            # if iteration_idx == 11:
            #     smot.update_tree_attributes(stoch_map_attr_name, debug=True)
            # else:
            smot.update_tree_attributes(stoch_map_attr_name)


    def _report_clado_map_issue_str(self) -> None:
        clado_issue_str = str()
        identical_cladogenetic_maps_dict = dict()
        nonwidespread_cladogenetic_maps_dict = dict()
        
        for it_idx, stoch_maps_on_tree \
                in self.stoch_maps_tree_dict.items():
            for nd_name, range_split_map \
                in stoch_maps_on_tree \
                    .clado_stoch_maps_dict.items():

                # range that splits is not widespread!
                # (we look at a bit pattern and if the sum of '1's is
                # <= 1, that means the range is an atomic region, i.e.,
                # not widespread)
                if sum(int(i) for i in \
                       range_split_map.from_state_bit_patt) <= 1:
                    if it_idx not in nonwidespread_cladogenetic_maps_dict:
                        nonwidespread_cladogenetic_maps_dict[it_idx] = [nd_name]

                    else:
                        nonwidespread_cladogenetic_maps_dict[it_idx].append(nd_name)

            # now we interrogate the StochMapsOnTree instance
            # for (repeated) identical RangeSplit maps
            if it_idx in stoch_maps_on_tree.identical_clado_stoch_maps_dict:
                nd_idx_list = stoch_maps_on_tree \
                    .identical_clado_stoch_maps_dict[it_idx]
                identical_cladogenetic_maps_dict[it_idx] = nd_idx_list
 
        if len(nonwidespread_cladogenetic_maps_dict) > 0 \
            or len(identical_cladogenetic_maps_dict) > 0:
                clado_issue_str = ("\nWARNING: Cladogenetic stochastic maps "
                                   "had issues:\n\n")
                
                clado_issue_str += ("\nIterations with repeated entries"
                                    " (and the affected nodes)\n")
                
                for it_idx, nd_list \
                        in identical_cladogenetic_maps_dict.items():
                    clado_issue_str += str(it_idx) \
                        + "\t" + ", ".join(i for i in nd_list) + "\n"
                    
                clado_issue_str += ("\nIterations with non-widespread parent ranges"
                                    " (and the affected nodes)\n")
                
                for it_idx, nd_list \
                        in nonwidespread_cladogenetic_maps_dict.items():
                    clado_issue_str += str(it_idx) \
                        + "\t" + ", ".join(i for i in nd_list) + "\n"
        
        if clado_issue_str != "":
            print(clado_issue_str)

    def _init_str_representation(self) -> None:
        self.str_representation = "Stochastic mapping collection"

        for it_idx, smot in self.stoch_maps_tree_dict.items():
            anag_node_names_str = ", ".join(
                smot.anag_stoch_maps_dict.keys())
            clado_node_names_str = ", ".join(
                smot.clado_stoch_maps_dict.keys())
            
            self.str_representation += "\n  Iteration " + str(it_idx) \
                + "\n    Nodes with anagenetic maps: " \
                    + anag_node_names_str \

            # preparing string for anagenetic events
            if len(smot.anag_stoch_maps_dict) > 0:
                range_expansion_count = 0
                range_contraction_count = 0

                for node_name, sm_list in \
                        smot.anag_stoch_maps_dict.items():
                    for sm in sm_list:
                        if isinstance(sm, RangeExpansion):
                            range_expansion_count += 1

                        if isinstance(sm, RangeContraction):
                            range_contraction_count += 1
                
                range_expansion_str = \
                    "\n      Number of range expansion (dispersal) events: " \
                    + str(range_expansion_count)
                    
                range_contraction_str = \
                    "\n      Number of range contraction (extinction) events: " \
                    + str(range_contraction_count)
                
                self.str_representation += \
                    range_expansion_str + range_contraction_str

            self.str_representation += "\n    Nodes with cladogenetic maps: " \
                + clado_node_names_str
            
            # preparing string for cladogenetic events
            if len(smot.clado_stoch_maps_dict) > 0:
                self.str_representation += \
                        "\n      Range size (of nodes above): "
                range_split_str = ""
                range_birth_str = ""
                first_element = True

                for splitting_node_name, rsb in \
                        smot.clado_stoch_maps_dict.items():
                    
                    # getting node names for nodes with cladogenetic maps
                    if first_element:
                        self.str_representation += \
                            str(rsb.size_of_splitting_node_range)

                        first_element = False

                    else:
                        self.str_representation += ", " + \
                            str(rsb.size_of_splitting_node_range)

                    # getting node names for b/w or within-speciation
                    if range_split_str == "" and rsb.range_split:
                        range_split_str = \
                            "\n      Between-region speciation nodes: " \
                            + splitting_node_name
                    
                    elif range_split_str != "" and rsb.range_split:
                        range_split_str += ", " + splitting_node_name

                    elif range_birth_str == "" and not rsb.range_split:
                        range_birth_str = \
                            "\n      Within-region speciation nodes: " \
                            + splitting_node_name
                    
                    elif range_birth_str != "" and not rsb.range_split:
                        range_birth_str += ", " + splitting_node_name

                self.str_representation += \
                    range_split_str + range_birth_str

    def __str__(self) -> str:
        return self.str_representation
    
    def iter_mcmc_generations(self) -> ty.Iterator[StochMapsOnTree]:
        for idx in self.sorted_it_idxs:
            yield self.stoch_maps_tree_dict[idx]


if __name__ == "__main__":
    
    # "examples/trees_maps_files/turtle.tre",
    trs = [pjr.read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree1.tre",
                                "read_tree",
                                 node_names_attribute="index",
                                 n_states=3,
                                 in_file=True)]
    
    # bit_patts = list()
    # with open("examples/trees_maps_files/kadua_range_label.csv", "r") as infile:
    #     infile.readline()

    #     for line in infile:
    #         tokens = line.rstrip().split(",")
    #         bit_patts.append(tokens[1])

    # unique_dict = dict()
    # with open("examples/trees_maps_files/geosse_maps.tsv", "r") as infile:
    #     infile.readline()  # skip first
    #     for line in infile:
    #         tokens = line.rstrip().split("\t")    
    #         if tokens[7] == "anagenetic":
    #             start_state = int(tokens[4])
    #             end_state = int(tokens[5])
    #             if start_state == 1 and end_state == 2:
    #                 print(tokens)
    #             if start_state == 2 and end_state == 1:
    #                 print(tokens)
    #             tup = (start_state, end_state)
    #             if tup not in unique_dict:
    #                 unique_dict[tup] = 1
    #             else:
    #                 unique_dict[tup] += 1
    # for k, v in unique_dict.items():
    #     print(k, v)

    # state2bit_lookup = pjbio.State2BitLookup(8, 2, geosse=True)
    n_chars = 2
    state2bit_lookup = pjbio.State2BitLookup(n_chars, 2, geosse=True)

    smap_coll = \
        StochMapsOnTreeCollection("examples/trees_maps_files/geosse_dummy_tree1_maps.tsv",
                           trs,
                           state2bit_lookup,
                           node_states_file_path="examples/trees_maps_files/geosse_dummy_tree1_tip_states.tsv",
                           stoch_map_attr_name="state")

    fig = plt.figure()

    ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    ax.patch.set_alpha(0.0)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # trees have been multiplied inside StochMapsOnTreeCollection
    i = 1 # 0
    # print("i =", i)
    # for nd in trs[i].tree.preorder_node_iter():
    #     print(nd.label, nd.state)
    #
    # if "nd4" in trs[i].at_dict:
    #     for at in trs[i].at_dict["nd4"]:
    #         print(at)
    # if "nd3" in trs[i].at_dict:
    #     for at in trs[i].at_dict["nd3"]:
    #         print(at)
    # if "nd2" in trs[i].at_dict:
    #     for at in trs[i].at_dict["nd2"]:
    #         print(at)
    # if "nd1" in trs[i].at_dict:
    #     for at in trs[i].at_dict["nd1"]:
    #         print(at)

    pjtr.plot_ann_tree(trs[i],
                       ax,
                       use_age=False,
                       start_at_origin=False,
                       sa_along_branches=False,
                       attr_of_interest="state")

    plt.show()

    # choose one:
    # print(stoch_mapcoll)
    #
    # or
    #
    # for smot in stoch_mapcoll.iter_mcmc_generations():s
    #     print(smot)
        # k: node index, v: stoch map
        # print(smot.anagenetic_stoch_maps_dict)
        # print(smot.cladogenetic_stoch_maps_dict)

    # for it_idx, smot in stoch_mapcoll.stoch_maps_tree_dict.items():
    #     print(smot)
    #     print(smot)
    #     widespread_range_size_dict = smot.get_clado_distribution()

    #     print("    Size of range being split, count:")
    #     for k, v in widespread_range_size_dict.items():
    #         print("    ", k, ", ", v, sep="")
