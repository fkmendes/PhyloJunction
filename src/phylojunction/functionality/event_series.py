import os
import math
import enum
import copy
import random
import collections
import typing as ty
import dendropy as dp
from natsort import natsorted
from numpy.random import choice

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.tree as pjt
import phylojunction.functionality.evol_event as pjev
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.biogeo as pjbio
import phylojunction.functionality.feature_io as pjfio
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class Hypothesis(enum.Enum):
    VICARIANCE = 0
    FOUNDER_EVENT = 1
    SPECIATION_BY_EXT = 2


class FromRegionSampler():

    _n_char: int
    _n_time_slices: int
    _param_log_dir: str

    # keys: iterations idx
    # values: 2d-list with parameter values
    _param_value_dict: ty.Dict[int, ty.List[ty.List[float]]]

    def __init__(self,
                 n_char: int,
                 param_log_dir: str,
                 file_name_prefix: str,
                 file_name_suffix: str,
                 param_name: str) -> None:

        self._n_char = n_char
        self._param_name = param_name
        self._time_slice_dict_list = list()

        # reading log files
        if not param_log_dir.endswith("/"):
            param_log_dir += "/"
        self._param_log_dir = param_log_dir

        log_fp_list = natsorted(
            [param_log_dir + f for f in os.listdir(param_log_dir) \
             if f.startswith(file_name_prefix) and \
             f.endswith(file_name_suffix + ".log")])
        self._n_time_slices = len(log_fp_list)

        # side-effect: populates
        #     (i)  self._n_its
        #     (ii) self._param_value_dict
        self.populate_param_value_dict(
            log_fp_list,
            file_name_prefix,
            file_name_suffix)

        # debugging
        # print(self._param_value_dict)

        print(("Finished reading .log files containing parameter "
               "values used to decode stochastic maps."))

        pass

    def populate_param_value_dict(self,
                                  param_log_fp_list: ty.List[str],
                                  prefix: str,
                                  suffix: str) -> None:
        for param_log_fp in param_log_fp_list:
            time_slice_idx = param_log_fp.\
                replace(prefix, "").\
                replace(suffix, "").\
                replace(".log", "").\
                replace(self._param_log_dir, "")

            # debugging
            # print('time_slice_idx', time_slice_idx)

            param_value_dict = dict()
            with open(param_log_fp, "r") as infile:
                tokens = infile.readline().rstrip().split("\t")

                # parameter name as key, column index as value
                param_name_idx_dict = dict()
                for region1_idx in range(1, self._n_char + 1):
                    for region2_idx in range(1, self._n_char + 1):
                        # name of parameter we care about
                        # (it assumes it is a triple-nested parameter
                        # 1D: epoch
                        # 2D: from region
                        # 3D: to region)
                        param_name = \
                            self._param_name \
                            + "[" + time_slice_idx + "]" \
                            + "[" + str(region1_idx) + "]" \
                            + "[" + str(region2_idx) + "]"

                        # this dictionary gives us the column index
                        # of the parameters we care about
                        param_name_idx_dict[param_name] = \
                            tokens.index(param_name)

                        # debugging
                        # print('param_name', param_name, 'pos', tokens.index(param_name))

                # reading rest of .log file
                for line in infile:
                    # empty matrix
                    param_val_mat = [[math.inf] * self._n_char for i in range(self._n_char)]
                    vals = line.rstrip().split("\t")
                    it = vals[0]

                    # we populate param_val_mat
                    for region1_idx in range(1, self._n_char + 1):
                        for region2_idx in range(1, self._n_char + 1):
                            param_name = \
                                self._param_name \
                                + "[" + time_slice_idx + "]" \
                                + "[" + str(region1_idx) + "]" \
                                + "[" + str(region2_idx) + "]"

                            val_idx = param_name_idx_dict[param_name]
                            param_val_mat[region1_idx-1][region2_idx-1] = float(vals[val_idx])

                    param_value_dict[int(it)] = param_val_mat

            # adding this time slice's dict to list
            self._time_slice_dict_list.append(param_value_dict)

    @property
    def time_slice_dict_list(self) -> \
            ty.List[ty.Dict[int, ty.List[ty.List[float]]]]:
        return self._time_slice_dict_list

    def sample_from_region_idx(self,
                               it_idx: int,
                               time_slice_idx: int,
                               potential_from_region_idx: ty.List[int],
                               to_region_idx: int) -> int:
        """Sample a region for an iteration.

        Parameter values in a .log file can be normalized and serve
        as weights for sampling a 'from' region for dispersal
        stochastic maps.

        Args:
            it_idx (int): Index of iteration we are looking at.
                This is parsed from the stochastic map table file.
            time_slice_idx (int): Index of time slice.
            potential_from_region_idx: List of indices for regions
                from which dispersal (range expansion) may have
                happened.
            to_region_idx: Index for region to which dispersal (range
                expansion) happened.

        Returns:
            (int): Sampled region index.
        """

        epoch_dict_it_param_mat = \
            self._time_slice_dict_list[time_slice_idx][it_idx]
        weight_list = list()
        for from_region_idx in potential_from_region_idx:
            weight_list.append(
                epoch_dict_it_param_mat[from_region_idx][to_region_idx]
            )

        tot = sum(weight_list)
        weight_list = [i / tot for i in weight_list]

        sampled_region_idx = \
            choice(potential_from_region_idx,
               1,
                p=weight_list).tolist()[0]

        return sampled_region_idx


class RegionDispersalAncestry():
    """Track all regions gained by dispersal

    This class keeps track of all regions gained by dispersal (range
    expansion) both (i) over a barrier, (ii) without a barrier. The
    purpose of this class is quantifying what proportion (in number
    of regions) of the splitting range had something to do with a
    dispersal over a barrier.
    """

    # keys are the region index that was gained at dispersal (range
    # expansion)
    # values are tuple containing:
    #   (i)  the 'to-region' index of the over-barrier dispersal
    #        that is the (direct or indirect) ancestor of the focal
    #        dispersal
    #   (ii) the time of the "seeding" over-barrier dispersal
    _regions_idx_gained_over_barrier: ty.Dict[int, ty.Tuple[int, float]]

    # each tuple corresponds to a "seeding" over-barrier dispersal,
    # where each of these is unique
    #
    # (these are the values of dictionary
    # _region_idx_gained_over_barrier, where they may appear
    # multiple times, for different gained regions)
    _over_barrier_seeding_dispersals: ty.Set[ty.Tuple[int, float]]

    # all region indices that were gained through standard no-barrier
    # dispersals
    _regions_idx_gained_w_out_barrier: ty.Set[int]

    _n_regions_gained_over_barrier: int
    _n_regions_gained_w_out_barrier: int

    def __init__(self) -> None:
        self._regions_idx_gained_over_barrier = set()
        self._over_barrier_seeding_dispersals = set()
        self._regions_idx_gained_w_out_barrier = set()

    @property
    def regions_idx_gained_over_barrier(self) -> \
            ty.Set[ty.Tuple[int, float]]:
        return self._regions_idx_gained_over_barrier

    @property
    def n_regions_gained_over_barrier(self) -> int:
        return len(self._regions_idx_gained_over_barrier)

    @property
    def n_regions_gained_w_out_barrier(self) -> int:
        return len(self._regions_idx_gained_w_out_barrier)


class EvolRelevantEventSeries:
    """Series of evolution-relevant events.

    This class has member methods that interrogate a series of event,
    truncating it depending on what the user specifies.

    """

    _event_list: ty.List[pjev.EvolRelevantEvent]
    _n_events: int
    _series_type: str
    _it_idx: int  # MCMC iteration from which we get an event

    # not initialized at instantiation, but populated
    # by external functions
    _trunc_event_list: ty.List[pjev.EvolRelevantEvent]
    _trunc_n_events: int
    _supported_hyp: Hypothesis

    def __init__(self,
                 event_list: ty.List[pjev.EvolRelevantEvent] = [],
                 series_type: str = "",
                 it_idx: int = -1) -> None:

        self._it_idx = it_idx
        self._series_type = series_type

        self.health_check()

        self._event_list = event_list
        self._n_events = len(event_list)

    def health_check(self) -> None:
        """Confirm validity of event list given series type."""

        if self._series_type == "speciation":
            last_event = self._event_list[-1]

            if not isinstance(last_event, pjsmap.RangeSplitOrBirth):
                raise ec.SequenceInvalidLastError(type(last_event))

    def truncate_upstream(self, truncator_fn: ty.Callable) -> None:
        """Deep-copy series and truncate copy upstream.

        Truncation of an event series is done according to some
        user-specified function. Truncation may be carried out, for
        example, at the first event that destabilizes a range (i.e.,
        the first "destabilizer" event).

        Side-effects:
            (i)  Populate self._trunc_event_list
            (ii) Populate self._trunc_n_events

        Parameters:
            truncator_fn (function): Method of static Truncate class
                that executes truncation.
        """

        self._trunc_event_list = truncator_fn(
            copy.deepcopy(self._event_list))

        self._trunc_n_events = len(self._trunc_event_list)


    @property
    def event_list(self) -> ty.List[pjev.EvolRelevantEvent]:
        return self._event_list

    def add_events(self, event_or_events: \
            ty.Union[pjev.EvolRelevantEvent,
            ty.List[pjev.EvolRelevantEvent]]):
        if isinstance(event_or_events, pjev.EvolRelevantEvent):
            self._event_list.append(event_or_events)

        elif isinstance(event_or_events, list):
            self._event_list += event_or_events

        else:
            exit("Cannot add event to event series. Exiting...")

    @property
    def n_events(self) -> int:
        return self._n_events

    @property
    def series_type(self) -> str:
        return self._series_type

    @property
    def supported_hyp(self) -> Hypothesis:
        return self._supported_hyp

    @supported_hyp.setter
    def supported_hyp(self, a_hyp: Hypothesis) -> None:
        if (a_hyp == Hypothesis.VICARIANCE or \
                a_hyp == Hypothesis.FOUNDER_EVENT) and \
                self._series_type != "speciation":
            raise ec.SequenceHypothesisTypeError(a_hyp, self._series_type)

        self._supported_hyp = a_hyp


class EvolRelevantEventSeriesTabulator():
    """Classify and tabulate event series.

    Parameters:
        ann_tr_list (AnnotatedTree): List of AnnotatedTree objects,
            one per iteration in 'smap_collection', on which character-
            state changes are being mapped.
        smap_collection (StochMapsOnTreeCollection): Collection of
            stochastic maps from multiple iterations.
        event_series_dict (dict): Dictionary with node labels as keys
            and a list of EvolRelevantEventSeries objects as values.
        hyp_support_dict (dict): Dictionary with Hypothesis objects
            as keys, and counts (int) as values.
    """

    _ann_tr_list: ty.List[pjt.AnnotatedTree]
    _smap_collection: pjsmap.StochMapsOnTreeCollection
    _geofeat_query: pjfio.GeoFeatureQuery

    # key of outer dict: node label
    # key of inner dict is the iteration index
    _event_series_dict: ty.Dict[str, ty.Dict[int, EvolRelevantEventSeries]]
    # same as above, but event series is truncated at its oldest end
    # at the first destabilizer dispersal
    _truncated_event_series_dict: ty.Dict[str, ty.Dict[int, EvolRelevantEventSeries]]

    # key of outer dict: node label
    # key of inner dict is the iteration index
    _region_dispersal_ancestry_dict: ty.Dict[str, ty.Dict[int, RegionDispersalAncestry]]

    _hyp_support_dict: ty.Dict[Hypothesis, int]
    _from_region_sampler: FromRegionSampler

    def __init__(self,
                 ann_tr_list: ty.List[pjt.AnnotatedTree],
                 smap_collection: pjsmap.StochMapsOnTreeCollection,
                 geofeat_query: pjfio.GeoFeatureQuery,
                 from_region_sampler: ty.Optional[FromRegionSampler] = None,
                 verbose: ty.Optional[bool] = False) -> None:

        self._event_series_dict = pjh.autovivify(2)
        self._truncated_event_series_dict = pjh.autovivify(2)
        self._region_dispersal_ancestry_dict = pjh.autovivify(2)
        self._from_region_sampler = from_region_sampler

        self._ann_tr_list = ann_tr_list

        self._smap_collection = smap_collection

        self._geofeat_query = geofeat_query

        # side-effect:
        # initializes self._event_series_dict
        # initializes self._trunc_event_series_dict
        #
        # range expansion stochastic maps are disambiguated in here
        if verbose:
            print("Beginning event series parsing:")

        self.initialize_event_series_dict()

        if verbose:
            print("  ... finished building all event series.")

        self.initialize_truncated_event_series_dict()

        if verbose:
            print("  ... finished truncating event series.")

        # this method further annotates the last dispersal
        # event of a truncated event series as stabilizing
        # the splitting range or not
        # self.initialize_region_dispersal_ancestry_dict()

        if verbose:
            print(("  ... finished tracking down the dispersal ancestry of every "
                   "gained region."))

        # for each event series, update each event component
        # with respect to its char_status_dict member --
        # this member keeps track of a region's connected or
        # disconnected status
        # self.update_event_char_member()
        # print("  Finished going through all event series and annotating it.")

        # for each truncated event series, update its supported_hyp
        # member depending on the nature of the truncated series
        # self.update_event_series_hyp_member()
        # print("  Finished classifying event series as supporting each hypothesis.")

        # side-effect:
        # initializes self._hyp_support_dict
        # self.tabulate_hyp_support()
        # print("  Finished tabulating hypothesis support.")

    def disambiguate_range_expansion(self,
                                     smap: pjsmap.StochMap,
                                     it_idx: int) -> None:
        """

        The side-effect of this method is to populate the
        member 'from_region_idx' of StochMap of RangeExpansion
        type.

        This disambiguation is necessary because....

        Args:
            smap (StochMap): StochMap object instance, carrying all
                info on a stochastic map.
            it_idx (int): Iteration index for which to disambiguate
                stochastic map. Note that the index is whatever is
                written in the stochastic map table and in the .log
                file of parameter values used to disambiguate the
                stochastic map. The indices for those iterations
                must match!
        """

        if smap.map_type == "expansion":
            bp = smap.from_state_bit_patt
            to_region_idx = smap.region_gained_idx
            potential_from_region_idx = \
                [idx for idx, b in enumerate(bp) if b == '1']

            # from 0 to (number of time slices - 1)
            time_slice_idx = self._geofeat_query.find_epoch_idx(smap.age)

            # if the class member exists, nothing needs to be done,
            # otherwise, we need to disambiguate it
            if smap.from_region_idx == None:
                # sample proportional to some scheme (e.g., proportional
                # to FIG rate scalers, m_d)
                if self._from_region_sampler != None:
                    sampled_idx = \
                        self._from_region_sampler.sample_from_region_idx(
                            it_idx,
                            time_slice_idx,
                            potential_from_region_idx,
                            to_region_idx)

                    # debugging
                    print('bit pattern', bp,
                          'potential_from_region_idx',
                          potential_from_region_idx,
                          'to_region_idx',
                          to_region_idx,
                          'sampled_idx',
                          sampled_idx)

                    # disambiguation!
                    smap.from_region_idx = sampled_idx

                else:
                    # disambiguation!
                    smap.from_region_idx = \
                        random.choice(potential_from_region_idx)

    def initialize_event_series_dict(self) -> None:

        def recursively_populate_event_series_dict(nd: dp.Node,
                                                   it_idx: int):
            """Populate self._event_series_dict recursively.

            Disambiguation of range expansion (dispersal) events
            happens here.
            """

            nd_name = nd.label
            event_series = EvolRelevantEventSeries()

            # retrieve list of events from parent node
            if nd.parent_node is not None:
                # and nd.parent_node.label in self._event_series_dict:
                parent_nd = nd.parent_node
                parent_nd_name = parent_nd.label

                is_root = parent_nd.num_child_nodes() == 1

                # NOTE: CAREFUL! this call here creates an empty dictionary even
                # if parent_nd_name is not a key in the outer _event_series_dict!
                # so we make sure to check that this node is younger than the root,
                # otherwise the following block adds the origin to the dictionary
                # (we don't want that, because there is no range splitting at the
                # origin)
                if not is_root:
                    parent_event_series = \
                        self._event_series_dict[parent_nd_name][it_idx]

                    # only if parent has a proper event series
                    # (it may be an empty dictionary if the parent node has no
                    # event series)
                    if isinstance(parent_event_series, EvolRelevantEventSeries):
                        # we will add events to those from the parent node
                        # (has to be deep copy!)
                        event_series = copy.deepcopy(parent_event_series)

            # do current node
            smap_on_tree = smap_coll.stoch_maps_tree_dict[it_idx]

            # only internal nodes (speciation!)
            if nd_name in smap_on_tree.clado_stoch_maps_dict:
                clado_smap = smap_on_tree.clado_stoch_maps_dict[nd_name]

                # only speciation events with range splitting
                if clado_smap.to_state2_bit_patt is not None:
                    # anagenetic
                    # NOTE: assumes stochastic maps are sorted in chronological order!!!
                    # (old first, young later)
                    if nd_name in smap_on_tree.anag_stoch_maps_dict:
                        anagenetic_smaps_list = smap_on_tree.anag_stoch_maps_dict[nd_name]

                        ############################################
                        # Further annotation of range expansions #
                        # (dispersals)                             #
                        ############################################
                        for ana_smap in anagenetic_smaps_list:
                            # first we disambiguate the source region
                            self.disambiguate_range_expansion(ana_smap,
                                                              it_idx)
                            # then we annotate the over barrier status
                            time_slice_idx = \
                                self._geofeat_query.find_epoch_idx(ana_smap.age)
                            conn_graph = \
                                self._geofeat_query.conn_graph_list[time_slice_idx]

                            if ana_smap.map_type == "expansion":
                                from_region_idx = ana_smap.from_region_idx
                                to_region_idx = ana_smap.region_gained_idx

                                are_connected = \
                                    conn_graph.are_connected(from_region_idx,
                                                             to_region_idx)
                                ana_smap.over_barrier = not are_connected

                        event_series.add_events(anagenetic_smaps_list)

                    # cladogenetic (has to be the last one in the event series)
                    event_series.add_events(clado_smap)

                # update or create value in dictionary
                self._event_series_dict[nd_name][it_idx] = event_series

            # recur
            for ch_nd in nd.child_nodes():
                recursively_populate_event_series_dict(ch_nd, it_idx)

        # iterating over each MCMC iteration when stochastic maps were logged
        for it_idx, smap in smap_coll.stoch_maps_tree_dict.items():
            ann_tr = smap.ann_tr
            seed_nd = smap.ann_tr.origin_node if smap.ann_tr.with_origin \
                else smap.ann_tr.root_node

            # populate self._self._event_series_dict
            #
            # key1: node name
            # value1 iteration dict:
            # key2: iteration index
            # value2 event series object
            recursively_populate_event_series_dict(seed_nd,
                                                   it_idx)

            # debugging
            for nd_label, it_event_series_dict in self._event_series_dict.items():
                # if nd_label in ("nd5"):
                for it_idx, event_series in it_event_series_dict.items():
                    if isinstance(event_series, EvolRelevantEventSeries):
                        if it_idx == 1:
                            print(nd_label)
                            for ev in event_series.event_list:
                                print(ev)
                print("\n")


    def initialize_truncated_event_series_dict(self) -> None:
        """Populate _truncated_event_series_dict

        This method visits all event series (for all internal nodes
        and all iterations), makes a deep copy of it, and truncates
        it at the first split-relevant over-barrier dispersal
        event.
        """

        for nd_label, it_event_series_dict in self._event_series_dict.items():
            rda = RegionDispersalAncestry()

            for it_idx, event_series in it_event_series_dict.items():
                # event_series will be an empty dictionary if range
                # did not split at node -- we do not care about those!
                # (even if subtending branch has maps!)
                if isinstance(event_series, EvolRelevantEventSeries):
                    smap_list = event_series.event_list

                    # we make deep copy because the annotation we carry out below
                    # (.split_relevant) will be different for each node in the tree,
                    # for the same dispersal event
                    new_smap_list = [copy.deepcopy(smap) for smap in smap_list]

                    has_seen_first_split_relevant_over_barrier = False
                    # first element must be both split-relevant and over barrier
                    relevant_smap_idx_list = list()

                    # root/origin will have an empty smap_list
                    if len(smap_list) > 0:
                        # first we find out the mutually exclusive sets of regions
                        range_split_smap = new_smap_list[-1] # last event is range split
                        ch1_bp = range_split_smap.to_state_bit_patt
                        ch2_bp = range_split_smap.to_state2_bit_patt
                        ch1_set = set([idx for idx, b in enumerate(ch1_bp) if b == "1"])
                        ch2_set = set([idx for idx, b in enumerate(ch2_bp) if b == "1"])

                        if not ch1_set.isdisjoint(ch2_set):
                            exit(("Error during truncation of event series. "
                                  "At range split, the two mutually exclusive ranges shared"
                                  " regions. Exiting..."))

                        # now we go through each dispersal and determine whether it is
                        # relevant to the focal range split
                        non_split_smap_list = smap_list[:-1]
                        for smap_idx, smap in enumerate(non_split_smap_list):
                            new_smap = copy.deepcopy(smap)

                            if smap.map_type == "expansion":
                                from_region_idx = smap.from_region_idx
                                to_region_idx = smap.region_gained_idx

                                # dispersal involves a from- and a to- regions
                                # that are being split between at speciation
                                if (from_region_idx in ch1_set and to_region_idx in ch2_set) or \
                                   (from_region_idx in ch2_set and to_region_idx in ch1_set):
                                        smap.split_relevant = True

                                        if smap.over_barrier:
                                            has_seen_first_split_relevant_over_barrier = True

                                        if has_seen_first_split_relevant_over_barrier:
                                            relevant_smap_idx_list.append(smap_idx)

                                else:
                                    smap.split_relevant = False

                                    if has_seen_first_split_relevant_over_barrier:
                                        relevant_smap_idx_list.append(smap_idx)

                            # extinctions
                            else:
                                if has_seen_first_split_relevant_over_barrier:
                                    relevant_smap_idx_list.append(smap_idx)

                    relevant_event_list = list()
                    relevant_smap_idx_list.append(-1)

                    # debugging
                    # print('nd_label', nd_label, 'relevant_smap_idx_list', relevant_smap_idx_list)

                    # getting only stochastic maps (shallow copies!) that matter and that
                    # will form the truncated event series
                    for relevant_smap_idx in relevant_smap_idx_list:
                        if len(smap_list) > 0:
                            relevant_event_list.append(smap_list[relevant_smap_idx])

                    # enter truncated event series into class member
                    if len(relevant_event_list) > 0:
                        truncated_event_series = \
                            EvolRelevantEventSeries(event_list=relevant_event_list)
                        self._truncated_event_series_dict[nd_label][it_idx] = \
                            truncated_event_series

                    # if there are no relevant stochastic maps, we use
                    # the original event series, which is probably empty anyway
                    else:
                        self._truncated_event_series_dict[nd_label][it_idx] = \
                            event_series

        # debugging
        # for nd_label, it_event_series_dict in self._truncated_event_series_dict.items():
        # # for nd_label, it_event_series_dict in self._event_series_dict.items():
        #     # if nd_label in ("nd9"):
        #     for it_idx, event_series in it_event_series_dict.items():
        #         if it_idx == 2:
        #             print(nd_label, 'it', it_idx)
        #             for ev in event_series.event_list:
        #                 print(ev)
        #         print("\n")

                # for ev in event_series.event_list:
                #     if ev.map_type == "expansion":
                #         # update rda members
                #         pass
                #
                #     elif ev.map_type == "contraction":
                #         # update rda members
                #         pass
                #
                # self._region_dispersal_ancestry_dict[nd_label][it_idx] = rda




    def update_event_char_member(self) -> None:
        """
        """

        pass

    def update_event_series_hyp_member(self) -> None:
        """
        """

        pass

    def tabulate_hyp_support(self) -> None:
        """
        """

        pass

    @property
    def hyp_support_dict(self) -> ty.Dict[Hypothesis, int]:
        return self._hyp_support_dict


##################################
# Event series parsing functions #
##################################

def truncate_at_split_relevant_event(event_series: EvolRelevantEventSeries):
    def _is_split_relevant_event():
        pass


if __name__ == "__main__":

    # n_chars = 2
    n_chars = 4

    state2bit_lookup = pjbio.State2BitLookup(n_chars, 2, geosse=True)

    n_states = state2bit_lookup.n_states

    # print(state2bit_lookup.bit2int_dict)
    # '1000': 0
    # '0100': 1
    # '0010': 2
    # '0001': 3
    # '1100': 4
    # '1010': 5
    # '0110': 6
    # '1001': 7
    # '0101': 8
    # '0011': 9
    # '1110': 10
    # '1101': 11
    # '1011': 12
    # '0111': 13
    # '1111': 14

    # ann_tr_list = [pjr.read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree2.tre",
    ann_tr_list = [pjr.read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree3.tre",
                                          "read_tree",
                                          node_names_attribute="index",
                                          n_states=n_states,
                                          in_file=True)]

    # node_states_file_path = "examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv"
    # pjsmap.StochMapsOnTreeCollection("examples/trees_maps_files/geosse_dummy_tree2_maps.tsv",
    smap_coll = \
        pjsmap.StochMapsOnTreeCollection("examples/trees_maps_files/geosse_dummy_tree3_maps.tsv",
                                         ann_tr_list,
                                         state2bit_lookup,
                                         node_states_file_path="examples/trees_maps_files/geosse_dummy_tree3_tip_states.tsv",
                                         stoch_map_attr_name="state")

    # frs = FromRegionSampler(
    #     n_chars,
    #     "examples/feature_files/two_regions_feature_set_event_series",
    #     "epoch_",
    #     "_rel_rates",
    #     "m_d"
    # )

    # feature_summary_fp = "examples/feature_files/two_regions_feature_set_event_series/feature_summary.csv"
    # age_summary_fp = "examples/feature_files/two_regions_feature_set_event_series/age_summary.csv"
    feature_summary_fp = "examples/feature_files/four_regions_feature_set_event_series/feature_summary.csv"
    age_summary_fp = "examples/feature_files/four_regions_feature_set_event_series/age_summary.csv"

    fc = pjfio.GeoFeatureCollection(
        feature_summary_fp,
        age_summary_fp=age_summary_fp)

    fq = pjfio.GeoFeatureQuery(fc)

    requirement_fn = \
        pjfio.GeoFeatureQuery.cb_feature_equals_value_is_connected(fc, 0, feat_id=1)

    # all members, including graph, are populated here!
    fq.populate_geo_cond_member_dicts("land_bridge", requirement_fn)

    event_series_tabulator = \
        EvolRelevantEventSeriesTabulator(
            ann_tr_list,
            smap_coll,
            fq#,
            #from_region_sampler=frs
        )