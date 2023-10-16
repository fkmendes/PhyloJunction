import sys
import os
import enum
import math
import numpy as np
import typing as ty

# pj imports
from phylojunction.data.tree import AnnotatedTree
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.biogeo as pjbio

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class StochMap():
    """Parent class for stochastic map"""

    n_regions: int = 0
    age: float = None
    time: float = None
    from_state: int = None
    from_state_bit_patt: str = None
    to_state: int = None
    to_state_bit_patt: str = None
    to_state2: int = None  # optional
    to_state2_bit_patt: str = None
    parent_node_name: str = ""
    child_node_name: str = ""
    child2_node_name: str = ""  # optional

    def __init__(self,
                 n_regions: int,
                 age: int,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 parent_node_name: str,
                 child_node_name: str,
                 time: ty.Optional[int] = None,
                 to_state2: ty.Optional[int] = None,
                 to_state2_bit_patt: ty.Optional[str] = None,
                 child2_node_name: ty.Optional[str] = None) -> None:
            
        self.age = age
        self.time = time
        self.from_state = from_state
        self.from_state_bit = from_state_bit_patt
        self.to_state = to_state
        self.to_state_bit = to_state_bit_patt
        self.parent_node_name = parent_node_name
        self.child_node_name = child_node_name
        self.n_regions = n_regions

        # if cladogenetic range split
        self.child2_node_name = child2_node_name
        self.to_state2 = to_state2
        self.to_state_bit2 = to_state2_bit_patt


class RangeExpansion(StochMap):
    """Stochastic map class for dispersal"""

    region_gained_idx: int = 0
    ancestral_over_barrier: bool = False
    range_expansion_erased_by_extinction: bool = False

    def __init__(self,
                 region_gained_idx: int,
                 n_regions: int,
                 age: int,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 parent_node_name: str,
                 child_node_name: str,
                 time: ty.Optional[int] = None):
        
        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         parent_node_name,
                         child_node_name,
                         time=time)

        self.region_gained_idx = region_gained_idx


    # setter
    def gained_region_is_lost_in_future(self):
        self.range_expansion_erased_by_extinction = True

    def dispersal_is_over_barrier(self):
        self.ancestral_over_barrier = True


class RangeContraction(StochMap):
    """Stochastic map class for local extinction"""

    region_lost_idx: int = 0

    def __init__(self,
                 region_lost_idx: int,
                 n_regions: int,
                 age: int,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 parent_node_name: str,
                 child_node_name: str,
                 time: int = None) -> None:
        
        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         parent_node_name,
                         child_node_name,
                         time=time)

        self.region_lost_idx = region_lost_idx


class RangeSplit(StochMap):
    """Stochastic class for cladogenetic event with range split"""

    def __init__(self,
                 n_regions: int,
                 age: int,
                 from_state: int,
                 from_state_bit_patt: str,
                 to_state: int,
                 to_state_bit_patt: str,
                 parent_node_name: str,
                 child_node_name: str,
                 child2_node_name: str = None,
                 to_state2: int = None,
                 to_state2_bit_patt: str = None,
                 time: int = None) -> None:
        
        super().__init__(n_regions,
                         age,
                         from_state,
                         from_state_bit_patt,
                         to_state,
                         to_state_bit_patt,
                         parent_node_name,
                         child_node_name,
                         child2_node_name=child2_node_name,
                         to_state2=to_state2,
                         to_state2_bit_patt=to_state2_bit_patt,
                         time=time)


class StochMapsOnTree():

    anagenetic_stoch_maps_on_tree_dict: \
        ty.Dict[str, ty.List[StochMap]] = dict()  # node name as key
    cladogenetic_stoch_maps_on_tree_dict: \
        ty.Dict[str, StochMap] = dict()  # node name as key
    iteration_idx: int = None
    str_representation: str = ""
    state2bit_lookup: pjbio.State2BitLookup = None
    n_regions: int

    def __init__(self,
                 stoch_map_str_list: ty.List[str],
                 iteration_idx: int,
                 ann_tr: AnnotatedTree,
                 state2bit_lookup: pjbio.State2BitLookup) -> None:
        
        self.iteration_idx = iteration_idx
        self.state2bit_lookup = state2bit_lookup
        self.n_regions = self.state2bit_lookup.n_char
        
        # initializes
        #     stoch_maps_on_tree_dict
        self._read_stoch_map_str_list(stoch_map_str_list)

        # sets critical properties in stoch maps
        # self._preorder_traverse_tree_annotating_maps(ann_tr)

    # internal
    def _read_stoch_map_str_list(self, stoch_map_str_list) -> None:
        
        def _which_region_changed(from_state: int, to_state: int) \
                -> ty.Tuple[bool, int]:
            pass

        nd_idx, br_start_age, br_end_age, from_state, to_state, \
        event_age, event_type, parent_idx, ch1_idx, ch2_idx \
            = stoch_map_str_list
        from_state_bit_patt = self.state2bit_lookup.get_bit(int(from_state))
        to_state_bit_patt = self.state2bit_lookup.get_bit(int(to_state))
        
        if event_type == "anagenetic":
            seen_different = False
            stoch_map: StochMap = None

            print("from, to")
            print(from_state_bit_patt)
            print(to_state_bit_patt)
            print(stoch_map_str_list)

            for idx, b_from, in enumerate(from_state_bit_patt):
                b_to = to_state_bit_patt[idx]

                if seen_different and b_from != b_to:
                    exit("ERROR: Found more than one different bit in bit patterns" +
                          b_from + " " + b_to + ". Exiting.")
                
                elif not seen_different and b_from != b_to:
                    seen_different = True
                    
                    if b_from == "0" and b_to == "1":
                        # TODO: later replace idx with node name
                        stoch_map = RangeExpansion(idx,
                                                   self.n_regions,
                                                   event_age,
                                                   int(from_state),
                                                   b_from,
                                                   int(to_state),
                                                   b_to,
                                                   parent_idx,
                                                   ch1_idx)
                    
                    elif b_from == "1" and b_to == "0":
                        # TODO: later replace idx with node name
                        stoch_map = RangeContraction(self.n_regions,
                                                    event_age,
                                                    int(from_state),
                                                    b_from,
                                                    int(to_state),
                                                    b_to,
                                                    parent_idx,
                                                    ch1_idx,
                                                    ch2_idx)
                        
            # TODO: later, replace idx with node name
            if nd_idx not in self.anagenetic_stoch_maps_on_tree_dict:
                self.anagenetic_stoch_maps_on_tree_dict[nd_idx] = [stoch_map]

            else:
                self.anagenetic_stoch_maps_on_tree_dict[nd_idx].append(stoch_map)

        elif event_type == "cladogenetic":
            if nd_idx not in self.cladogenetic_stoch_maps_on_tree_dict:
                stoch_map = RangeSplit(self.n_regions,
                                       event_age,
                                       int(from_state),
                                       b_from,
                                       int(to_state),
                                       b_to,
                                       parent_idx,
                                       ch1_idx,
                                       ch2_idx)
                
                self.anagenetic_stoch_maps_on_tree_dict[nd_idx] = stoch_map

            else:
                self.anagenetic_stoch_maps_on_tree_dict[nd_idx] \
                    .to_state2 = to_state
                self.anagenetic_stoch_maps_on_tree_dict[nd_idx] \
                    .to_state2_bit_patt = to_state_bit_patt

                # TODO: set child2 stuff
                # self.anagenetic_stoch_maps_on_tree_dict[nd_idx].update()

    def _preorder_traverse_tree_annotating_maps(self, ann_tr):
        pass

    def _init_str_representation(self) -> None:
        pass

    def stoch_maps_summary(self) -> str:
        pass


class StochMapCollection():

    stoch_maps_tree_list: ty.Dict[str, ty.List[str]] = dict()

    def __init__(self,
                 stoch_maps_file_path: str,
                 ann_tr: AnnotatedTree,
                 state2bit_lookup: pjbio.State2BitLookup) -> None:

        # initializes stoch_maps_tree_list
        self._read_stoch_maps_file(stoch_maps_file_path,
                                   ann_tr,
                                   state2bit_lookup)


    # internal
    def _read_stoch_maps_file(self,
                              stoch_maps_file_path: str,
                              ann_tr: AnnotatedTree,
                              state2bit_lookup: pjbio.State2BitLookup) -> None:

        if not os.path.isfile(stoch_maps_file_path):
            exit("ERROR: Could not find " + stoch_maps_file_path + ". Exiting.")
        
        with open(stoch_maps_file_path, "r") as infile:
            infile.readline() # skip header

            for line in infile:
                stoch_map_on_tree_str_list = line.rstrip().split("\t")
                iteration_idx = stoch_map_on_tree_str_list[0]

                stoch_maps_on_tree = \
                    StochMapsOnTree(stoch_map_on_tree_str_list[1:],
                                    iteration_idx,
                                    ann_tr,
                                    state2bit_lookup)
                
                # print(stoch_maps_on_tree.summary())


if __name__ == "__main__":
    
    tr = pjr.read_nwk_tree_str(
        "examples/trees_maps_files/turtle.tre",
        "read_tree",
        node_names_attribute="index")
    
    state2bit_lookup = pjbio.State2BitLookup(8, 2)
    
    stoch_mapcoll = \
        StochMapCollection("examples/trees_maps_files/turtle_maps.tsv",
                           tr,
                           state2bit_lookup)