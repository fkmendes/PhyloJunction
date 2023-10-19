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
    age: float
    time: float
    from_state: int
    from_state_bit_patt: str
    to_state: int
    to_state_bit_patt: str
    to_state2: int
    to_state2_bit_patt: str
    parent_node_name: str
    child_node_name: str
    child2_node_name: str

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
        self.from_state_bit_patt = from_state_bit_patt
        self.to_state = to_state
        self.to_state2_bit_patt = to_state_bit_patt
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

    # class members
    anagenetic_stoch_maps_dict: \
        ty.Dict[str, ty.List[StochMap]]  # node name as key
    cladogenetic_stoch_maps_dict: \
        ty.Dict[str, StochMap]  # node name as key
    iteration_idx: int
    str_representation: str
    state2bit_lookup: pjbio.State2BitLookup
    n_regions: int

    # for summary and printing
    n_total_maps: int
    maps_counts_dict: ty.Dict[ty.Tuple[str], int]
    n_anagenetic_maps: int
    n_range_exp_maps: int
    n_range_cont_maps: int
    pctg_anagenetic_maps: float
    forbidden_higher_order_maps_counts_dict: \
        ty.Dict[ty.Tuple[str], int]
    n_forbidden_higher_order_anagenetic_maps: int
    pctg_forbidden_higher_order_anagenetic_maps: float
    n_different_clado_maps: int
    n_identical_clado_maps: int

    # iteration
    index: int
    
    def __init__(self,
                 iteration_idx: int,
                 ann_tr: AnnotatedTree,
                 state2bit_lookup: pjbio.State2BitLookup) -> None:
        
        self.iteration_idx = iteration_idx
        self.state2bit_lookup = state2bit_lookup
        self.n_regions = self.state2bit_lookup.n_char

        self.anagenetic_stoch_maps_dict = dict()
        self.cladogenetic_stoch_maps_dict = dict()
        self.identical_cladogenetic_stoch_maps_dict = dict()
        
        # summary / counts
        self.maps_counts_dict = dict()
        self.forbidden_higher_order_maps_counts_dict = dict()
        self.n_total_maps = 0
        self.n_anagenetic_maps = 0
        self.n_range_exp_maps = 0
        self.n_range_cont_maps = 0
        self.pctg_anagenetic_maps = 0.0
        self.n_forbidden_higher_order_anagenetic_maps = 0
        self.pctg_forbidden_higher_order_anagenetic_maps = 0.0
        self.n_different_clado_maps = 0
        self.pctg_diff_clado_maps = 0.0
        self.n_identical_clado_maps = 0
        self.pctg_ident_clado_maps = 0.0

        # self._init_str_representation()

        # sets critical properties in stoch maps
        # self._preorder_traverse_tree_annotating_maps(ann_tr)

    # updates
    #     self.n_total_maps
    #     self.anagenetic_stoch_maps_on_tree_dict
    #     self.map_counts_dict
    #     self.n_anagenetic_maps
    #     self.forbidden_higher_order_maps_counts_dict
    #     self.n_forbidden_higher_order_anagenetic_maps
    def add_map(self, stoch_map_str_list) -> None:
        def _which_region_changed(from_state: int, to_state: int) \
                -> ty.Tuple[bool, int]:
            pass

        # print("n_forbidden before add_map=", self.n_forbidden_higher_order_anagenetic_maps)
        nd_idx, br_start_age, br_end_age, from_state, to_state, \
        event_age, event_type, parent_idx, ch1_idx, ch2_idx \
            = stoch_map_str_list
        from_state_bit_patt = self.state2bit_lookup.get_bit(int(from_state))
        to_state_bit_patt = self.state2bit_lookup.get_bit(int(to_state))

        if event_type == "anagenetic":
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

            # getting counts for printing later #
            if is_higher_order:
                # populate forbidden maps counts dict
                if (from_state_bit_patt, to_state_bit_patt) \
                        not in self.forbidden_higher_order_maps_counts_dict:
                    self.forbidden_higher_order_maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] = 1
                
                else:
                    self.forbidden_higher_order_maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] += 1

                self.n_forbidden_higher_order_anagenetic_maps += 1
            
            else:
                # populate maps counts dict
                if (from_state_bit_patt, to_state_bit_patt) \
                        not in self.maps_counts_dict:
                    self.maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] = 1
                
                else:
                    self.maps_counts_dict \
                        [(from_state_bit_patt, to_state_bit_patt)] += 1
                    
                self.n_anagenetic_maps += 1
            # done with getting counts #

                # initializing maps, finally
                if is_expansion:
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
                    
                    self.n_range_exp_maps += 1

                else:
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
                    
                    self.n_range_cont_maps += 1
                        
                # TODO: later, replace idx with node name
                if nd_idx not in self.anagenetic_stoch_maps_dict:
                    self.anagenetic_stoch_maps_dict[nd_idx] = [stoch_map]

                else:
                    self.anagenetic_stoch_maps_dict[nd_idx].append(stoch_map)

        elif event_type == "cladogenetic":
            if nd_idx not in self.cladogenetic_stoch_maps_dict:
                self.n_total_maps += 1
                self.n_different_clado_maps += 1

                stoch_map = RangeSplit(self.n_regions,
                                       event_age,
                                       int(from_state),
                                       from_state_bit_patt,
                                       int(to_state),
                                       to_state_bit_patt,
                                       parent_idx,
                                       ch1_idx,
                                       ch2_idx)
                
                self.cladogenetic_stoch_maps_dict[nd_idx] = stoch_map

            else:
                # just counting here
                if self.cladogenetic_stoch_maps_dict[nd_idx].from_state \
                    == int(from_state) and \
                    self.cladogenetic_stoch_maps_dict[nd_idx].to_state \
                        == int(to_state):
                    
                    if self.iteration_idx \
                            not in self.identical_cladogenetic_stoch_maps_dict:
                        self.identical_cladogenetic_stoch_maps_dict \
                            [self.iteration_idx] = [nd_idx]
                    
                    else:
                        self.identical_cladogenetic_stoch_maps_dict \
                            [self.iteration_idx].append(nd_idx)
                    
                    self.n_identical_clado_maps += 1
                    self.n_different_clado_maps -= 1

                # updating the cladogenetic map with info from the second child
                # else:
                self.cladogenetic_stoch_maps_dict[nd_idx] \
                    .to_state2 = int(to_state)
                self.cladogenetic_stoch_maps_dict[nd_idx] \
                    .to_state2_bit_patt = to_state_bit_patt
                    
                    # debugging
                    # range_split = self.cladogenetic_stoch_maps_on_tree_dict[nd_idx]
                    # if range_split.from_state == 2 and range_split.to_state == 0 and range_split.to_state2 == 0:
                    #     exit("found weird cladogenetic event")
                    # if range_split.from_state == 2 and range_split.to_state == 1 and range_split.to_state2 == 1:
                    #     exit("found weird cladogenetic event")
                    # if range_split.from_state == 2 and range_split.to_state == 0 and range_split.to_state2 == 1:
                    #     print("found good cladogenetic event")
                    # if range_split.from_state == 2 and range_split.to_state == 1 and range_split.to_state2 == 0:
                    #     print("found good cladogenetic event")


    def get_clado_distribution(self) -> ty.Dict[int, int]:
        widespread_range_sizes_dict = dict()
        for nd_name, clado_stoch_map \
                in self.cladogenetic_stoch_maps_dict.items():
            
            # print(clado_stoch_map.from_state, clado_stoch_map.from_state_bit_patt)
            
            widespread_range_size = \
                sum(int(i) for i in list(clado_stoch_map.from_state_bit_patt))
            
            if not widespread_range_size in widespread_range_sizes_dict:
                widespread_range_sizes_dict[widespread_range_size] = 1

            else:
                widespread_range_sizes_dict[widespread_range_size] += 1

        return widespread_range_sizes_dict

    # internal
    def _preorder_traverse_tree_annotating_maps(self, ann_tr):
        pass

    def _init_str_representation(self) -> None:
        self._stoch_maps_summary()

        self.str_representation = "Stochastic maps in MCMC iteration " \
            + str(self.iteration_idx) + ":"
        
        self.str_representation += \
            "\n    Number of anagenetic changes = " \
            + str(self.n_anagenetic_maps) \
            + " (" + str(self.pctg_anagenetic_maps) + "%)" \
            + "\n        Range contractions = " + str(self.n_range_cont_maps) \
            + "\n        Range expansions = " + str(self.n_range_exp_maps) \
            + "\n    Number of higher-order anagenetic changes = " \
            + str(self.n_forbidden_higher_order_anagenetic_maps) \
            + " (" + str(self.pctg_forbidden_higher_order_anagenetic_maps) + "%)" \
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
        self.pctg_anagenetic_maps = round(100.0 * (self.n_anagenetic_maps \
            / self.n_total_maps), 2)
        
        self.pctg_forbidden_higher_order_anagenetic_maps = \
            round(100.0 * (self.n_forbidden_higher_order_anagenetic_maps \
            / self.n_total_maps), 2)
        
        self.pctg_diff_clado_maps = round(100.0 * (self.n_different_clado_maps \
            / self.n_total_maps), 2)
        
        self.pctg_ident_clado_maps = round(100.0 * (self.n_identical_clado_maps \
            / self.n_total_maps), 2)

    def __str__(self) -> str:
        self._init_str_representation()

        return self.str_representation


class StochMapCollection():

    stoch_maps_tree_dict: ty.Dict[int, StochMapsOnTree] = dict()

    def __init__(self,
                 stoch_maps_file_path: str,
                 ann_tr: AnnotatedTree,
                 state2bit_lookup: pjbio.State2BitLookup) -> None:

        # initializes stoch_maps_tree_list
        self._read_stoch_maps_file(stoch_maps_file_path,
                                   ann_tr,
                                   state2bit_lookup)
        
        self._init_clado_map_issue_str()

    # internal
    def _read_stoch_maps_file(self,
                              stoch_maps_file_path: str,
                              ann_tr: AnnotatedTree,
                              state2bit_lookup: pjbio.State2BitLookup) -> None:

        if not os.path.isfile(stoch_maps_file_path):
            exit("ERROR: Could not find " + stoch_maps_file_path + ". Exiting.")
        
        with open(stoch_maps_file_path, "r") as infile:
            infile.readline() # skip header

            # iterating over MCMC iterations logging stochastic maps
            for line in infile:
                stoch_map_on_tree_str_list = line.rstrip().split("\t")
                iteration_idx = int(stoch_map_on_tree_str_list[0])

                if iteration_idx not in self.stoch_maps_tree_dict: 
                    self.stoch_maps_tree_dict[iteration_idx] \
                        = StochMapsOnTree(iteration_idx,
                                          ann_tr,
                                          state2bit_lookup)      

                # we add this map to it
                self.stoch_maps_tree_dict[iteration_idx] \
                    .add_map(stoch_map_on_tree_str_list[1:])
                
                # print("just added map", iteration_idx, "n_forbidden =", self.stoch_maps_tree_dict[iteration_idx].n_forbidden_higher_order_anagenetic_maps)

    def _init_clado_map_issue_str(self) -> None:
        clado_issue_str = str()
        single_entry_cladogenetic_maps = 0
        single_entry_cladogenetic_maps_dict = dict()
        identical_cladogenetic_maps_dict = dict()
        nonwidespread_cladogenetic_maps_dict = dict()
        
        for it_idx, stoch_maps_on_tree \
                in self.stoch_maps_tree_dict.items():
            for nd_name, range_split_map \
                in stoch_maps_on_tree \
                    .cladogenetic_stoch_maps_dict.items():
                
                # did not find a second line for cladogenetic map
                # so the information from child 2 is missing!
                if range_split_map.to_state2 == None:
                    single_entry_cladogenetic_maps += 1

                    if it_idx not in single_entry_cladogenetic_maps_dict:
                        single_entry_cladogenetic_maps_dict[it_idx] = [nd_name]

                    else:
                        single_entry_cladogenetic_maps_dict[it_idx].append(nd_name)

                # range that splits is not widespread!
                if sum(int(i) for i in \
                       range_split_map.from_state_bit_patt) <= 1:
                    if it_idx not in nonwidespread_cladogenetic_maps_dict:
                        nonwidespread_cladogenetic_maps_dict[it_idx] = [nd_name]

                    else:
                        nonwidespread_cladogenetic_maps_dict[it_idx].append(nd_name)

            # now we interrogate the StochMapsOnTree instance
            # for (repeated) identical RangeSplit maps
            if it_idx in stoch_maps_on_tree.identical_cladogenetic_stoch_maps_dict:
                nd_idx_list = stoch_maps_on_tree \
                    .identical_cladogenetic_stoch_maps_dict[it_idx]
                identical_cladogenetic_maps_dict[it_idx] = nd_idx_list
 
        if len(single_entry_cladogenetic_maps_dict) > 0 \
            or len(identical_cladogenetic_maps_dict) > 0:
                clado_issue_str = ("\nThere were issues with the "
                                   "cladogenetic stochastic maps\n\n"
                                   "Iterations with single entries "
                                   "(and the affected nodes)\n")
                
                for it_idx, nd_list \
                        in single_entry_cladogenetic_maps_dict.items():
                    clado_issue_str += str(it_idx) \
                        + "\t" + ", ".join(i for i in nd_list) + "\n"
                
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
        
        print(clado_issue_str)

            





if __name__ == "__main__":
    
    tr = pjr.read_nwk_tree_str(
        "examples/trees_maps_files/turtle.tre",
        "read_tree",
        node_names_attribute="index")
    
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

    state2bit_lookup = pjbio.State2BitLookup(8, 2, geosse=True)
    # first = False
    # for i, (k, v) in enumerate(state2bit_lookup.int2bit_dict.items()):
    #     print(k, v)
    
    
    stoch_mapcoll = \
        StochMapCollection("examples/trees_maps_files/turtle_maps.tsv",
        # StochMapCollection("examples/trees_maps_files/geosse_maps.tsv",
                           tr,
                           state2bit_lookup)

    # for smot in stoch_mapcoll.stoch_maps_tree_dict.values():
    #     print(smot)
    #     widespread_range_size_dict = smot.get_clado_distribution()

    #     print("    Size of range being split, count:")
    #     for k, v in widespread_range_size_dict.items():
    #         print("    ", k, ", ", v, sep="")
