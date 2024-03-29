import os
import math
import enum
import copy
import random
import typing as ty
import dendropy as dp
from natsort import natsorted
from numpy.random import choice
from bisect import insort

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
    VICARIANCE = "Vicariance"
    FOUNDER_EVENT = "Founder event"
    SPECIATION_BY_EXT = "Speciation by extinction"
    AMBIGUOUS = "Ambiguous"

    def __str__(self) -> str:
        return str(self.value)


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
            # print("time_slice_idx", time_slice_idx)

            param_value_dict = dict()
            with open(param_log_fp, "r") as infile:
                tokens = infile.readline().rstrip().split("\t")

                # parameter name as key, column index as value
                param_name_idx_dict = dict()
                for region1_idx in range(1, self._n_char + 1):
                    for region2_idx in range(1, self._n_char + 1):
                        # name of parameter we care about
                        # (it assumes it is a triple-nested parameter
                        # 1D: time slice (i.e., epoch)
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

            # NOTE:
            # everything in the event series parsing we do, we want
            # old -> young, but note that a larger age means older!
            # so we need to place higher time_slice_idx in the beginning
            # of _time_slice_dict_list, otherwise this object will have
            # young -> old parameter values that are later used in
            # weighing during disambiguation
            self._time_slice_dict_list.insert(0, param_value_dict)

    @property
    def time_slice_dict_list(self) -> \
            ty.List[ty.Dict[int, ty.List[ty.List[float]]]]:
        return self._time_slice_dict_list

    @property
    def n_char(self) -> int:
        return self._n_char

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
            time_slice_idx (int): Index of time slice, with 0 being
                the present epoch.
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

    Parameters:
        event_list (list):
        n_events (int):
        it_idx (int):
        trunc_event_list (list):
    """

    _event_list: ty.List[pjev.EvolRelevantEvent]
    _n_events: int
    _it_idx: int  # MCMC iteration from which we get an event

    # not initialized at instantiation, but populated
    # by external functions
    _trunc_event_list: ty.List[pjev.EvolRelevantEvent]
    _trunc_n_events: int
    _supported_hyp: Hypothesis

    def __init__(self,
                 event_list: ty.Optional[ty.List[pjev.EvolRelevantEvent]] = None,
                 it_idx: int = -1) -> None:

        self._it_idx = it_idx

        # self.health_check()

        if event_list is None:
            self._event_list = list()

        else:
            self._event_list = event_list

        self._n_events = len(self._event_list)
        self._supported_hyp = None
        self._str_representation = "Event series"

    # def health_check(self) -> None:
    #     """Confirm validity of event list given series type."""
    #
    #     if self._series_type == "speciation":
    #         last_event = self._event_list[-1]
    #
    #         if not isinstance(last_event, pjsmap.RangeSplitOrBirth):
    #             raise ec.SequenceInvalidLastError(type(last_event))


    @property
    def event_list(self) -> ty.List[pjev.EvolRelevantEvent]:
        return self._event_list

    @property
    def trunc_event_list(self) -> ty.List[pjev.EvolRelevantEvent]:
        return self._trunc_event_list

    @trunc_event_list.setter
    def trunc_event_list(self,
                         trunc_event_list: ty.List[pjev.EvolRelevantEvent]) \
            -> None:
        self._trunc_event_list = trunc_event_list

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
        return len(self._event_list)

    # @property
    # def series_type(self) -> str:
    #     return self._series_type

    @property
    def supported_hyp(self) -> Hypothesis:
        return self._supported_hyp

    @supported_hyp.setter
    def supported_hyp(self, a_hyp: Hypothesis) -> None:
        hyp2replace = ""
        if self._supported_hyp is not None:
            hyp2replace = self._supported_hyp.value

        self._supported_hyp = a_hyp
        self._update_str_representation_hyp(hyp2replace)

    def _update_str_representation_hyp(self, hyp2replace: str = "") -> None:
        if hyp2replace not in (None, ""):
            self._str_representation = \
                self._str_representation.replace(hyp2replace,
                                                 self._supported_hyp.value)
            # debugging
            # print("replacing", hyp2replace, "with", self._supported_hyp.value)

        else:
            self._str_representation += "\n  Supported hypothesis: " \
                                        + self._supported_hyp.value

    def __str__(self) -> str:
        return self._str_representation


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
        directed_edges (bool): Flag specifying if edges in connectivity
            graph are to be treated as directed. Defaults to 'False'.
    """

    _ann_tr_list: ty.List[pjt.AnnotatedTree]
    _smap_collection: pjsmap.StochMapsOnTreeCollection
    _geofeat_query: pjfio.GeoFeatureQuery

    # key of outer dict: node label
    # key of inner dict is the iteration index
    _event_series_dict: ty.Dict[str, ty.Dict[int, EvolRelevantEventSeries]]
    # # same as above, but event series is truncated at its oldest end
    # # at the first destabilizer dispersal
    # _truncated_event_series_dict: ty.Dict[str, ty.Dict[int, EvolRelevantEventSeries]]

    # key of outer dict: node label
    # key of inner dict is the iteration index
    _region_dispersal_ancestry_dict: ty.Dict[str, ty.Dict[int, RegionDispersalAncestry]]

    _hyp_support_dict: ty.Dict[Hypothesis, int]
    _from_region_sampler: FromRegionSampler
    _n_char: int # number of regions

    _directed_edges: bool

    def __init__(self,
                 ann_tr_list: ty.List[pjt.AnnotatedTree],
                 smap_collection: pjsmap.StochMapsOnTreeCollection,
                 geofeat_query: ty.Optional[pjfio.GeoFeatureQuery] = None,
                 from_region_sampler: ty.Optional[FromRegionSampler] = None,
                 directed_edges: ty.Optional[bool] = False,
                 verbose: ty.Optional[bool] = False) -> None:

        self._event_series_dict = pjh.autovivify(2)
        # self._truncated_event_series_dict = pjh.autovivify(2)
        self._region_dispersal_ancestry_dict = pjh.autovivify(2)

        self._from_region_sampler = from_region_sampler
        self._n_char = from_region_sampler.n_char

        self._ann_tr_list = ann_tr_list

        self._smap_collection = smap_collection

        self._geofeat_query = geofeat_query

        self._directed_edges = directed_edges

        # side-effect:
        # initializes self._event_series_dict
        # initializes self._trunc_event_series_dict
        #
        # range expansion stochastic maps are disambiguated in here
        if verbose:
            print("Beginning event series parsing:")

        self.initialize_event_series_dict()

        if verbose:
            print("  ... finished reading stochastic maps.")

        # model may not have time-het paleogeo. features
        if self._geofeat_query is not None:
            self.add_paleogeo_event_series_dict(
                self._geofeat_query.geo_cond_name)

            if verbose:
                print("  ... finished reading paleogeographic events.")

        self.update_trunc_event_series()

        if verbose:
            print("  ... finished truncating event series.")

        self.annotate_event_series_hypothesis()

        if verbose:
            print(("  ... finished classifying event series according"
                   " to the hypothesis they support"))

    def disambiguate_range_expansion(self,
                                     smap: pjsmap.StochMap,
                                     it_idx: int) -> None:
        """Randomly disambiguate source region of range expansion.

        The side-effect of this method is to populate the
        member 'from_region_idx' of StochMap of RangeExpansion
        type.

        This disambiguation is necessary in cases where the source
        range includes multiple atomic regions, e.g., ABC -> ABCD.
        We do not know if the dispersal to D came from A, B or C.

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

        # this is already checked in the code calling this method,
        # but doing it again just to be sure...
        if smap.map_type == "expansion":
            bp = smap.from_state_bit_patt
            to_region_idx = smap.region_gained_idx
            potential_from_region_idx = \
                [idx for idx, b in enumerate(bp) if b == '1']

            # from 0 to (number of time slices - 1)
            # with 0 being the oldest
            time_slice_idx = self._geofeat_query.find_epoch_idx(smap.age)

            # if the class member exists, nothing needs to be done,
            # otherwise, we need to disambiguate it
            if smap.from_region_idx == None:
                # sample proportional to some scheme (e.g., proportional
                # to FIG rate scalers, m_d)
                #
                # note that time_slice_idx of 0 means present!
                # so inside sample_from_region_idx, we should take that
                # into account!
                if self._from_region_sampler != None:
                    sampled_idx = \
                        self._from_region_sampler.sample_from_region_idx(
                            it_idx,
                            time_slice_idx,
                            potential_from_region_idx,
                            to_region_idx)

                    # debugging
                    # print('bit pattern', bp,
                    #       'potential_from_region_idx',
                    #       potential_from_region_idx,
                    #       'to_region_idx',
                    #       to_region_idx,
                    #       'sampled_idx',
                    #       sampled_idx)

                    # disambiguation!
                    smap.from_region_idx = sampled_idx

                else:
                    # disambiguation!
                    smap.from_region_idx = \
                        random.choice(potential_from_region_idx)

    def is_range_expansion_split_relevant(self,
                                          from_region_idx: int,
                                          to_region_idx: int,
                                          splitting_range1_idxs: ty.Set[int],
                                          splitting_range2_idxs: ty.Set[int]) -> bool:
        """Determine if range expansion split-relevant.

        If the regions involved in the dispersal are
        both on the same side of a range split, they are not split-
        relevant.

        Args:
            from_region_idx (int):
            to_region_idx (int): Index of region receiving migrants.
            splitting_range1_idxs (set): Set of region indices on one side of
                range split.
            splitting_range2_idxs (set): Set of region indices on other side of
                range split.

        Returns:
            (bool): Boolean for whether range expasion is
                split-relevant.
        """

        # the concept of split-relevance does not make sense if
        # the range is not splitting!
        if splitting_range1_idxs == splitting_range2_idxs:
            return False

        elif (from_region_idx in splitting_range1_idxs \
              and to_region_idx in splitting_range2_idxs) \
                or \
                (from_region_idx in splitting_range2_idxs \
                 and to_region_idx in splitting_range1_idxs):
            return True

        return False

    def is_fragile_wrt_split(self,
                             splitting_range1_idxs: ty.Set[int],
                             splitting_range2_idxs: ty.Set[int],
                             conn_graph: pjfio.GeoGraph,
                             range_idxs: \
                                     ty.Optional[ty.List[int]] = None) \
            -> bool:
        """Determine if splitting range is fragile.

        Checks that every region in range1_idxs is in a different
        communicating class from every region in range2_idxs.

        Args:
            splitting_range1_idxs (list): List of indices (int) of one of
                ranges resulting from the split.
            splitting_range2_idxs (list): List of indices (int) of the other
                range resulting from the split.
            conn_graph (GeoGraph): Connectivity graph, with each
                node being a region, and each edge representing
                the possibility of migration between two regions
                (i.e., gene flow).
            range_idxs (list, optional). List of indices
                of all regions constituting a (i) expanding range, (ii)
                contracting range, or (iii) contracted range. Defaults
                to None.
        """

        def find_edge_pairwise(range1_idxs: ty.Set[int],
                               range2_idxs: ty.Set[int],
                               conn_graph: pjfio.GeoGraph):

            for region_idx1 in range1_idxs:
                for region_idx2 in range2_idxs:
                    edge_one_way = (region_idx1, region_idx2)
                    edge_another_way = (region_idx2, region_idx1)

                    # if a single edge is found between comm
                    # classes, the range is not fragile
                    if edge_one_way in conn_graph.edge_set or \
                            edge_another_way in conn_graph.edge_set:
                        return True

            return False

        # this part of the method finds out if range at speciation
        # is fragile at the cladogenetic event
        if range_idxs is None:
            # if there are no edges at all
            # the range must be unstable
            if len(conn_graph.edge_set) == 0:
                return True

            # if both children have the same range
            # the range cannot be unstable, it just does
            # not make sense
            if splitting_range1_idxs == splitting_range2_idxs:
                return False

            # if there are edges, we check pairwise
            at_least_one_edge = \
                find_edge_pairwise(splitting_range1_idxs,
                                   splitting_range2_idxs,
                                   conn_graph)
            is_fragile = not at_least_one_edge

            return is_fragile

        # this part of the method is used to determine if
        # a range prior to a dispersal is already fragile
        # with respect to a splitting event happening in
        # the future along this branch
        elif range_idxs is not None:
            range1_idx = set([])
            range2_idx = set([])

            for region_idx in range_idxs:
                if region_idx in splitting_range1_idxs:
                    range1_idx.add(region_idx)

                elif region_idx in splitting_range2_idxs:
                    range2_idx.add(region_idx)

            # debugging
            # print("range1_idx", range1_idx, "range2_idx", range2_idx)

            # at least one region on both sides of the split
            # must be occupied for a range to be considered
            # 'previously fragile'
            if len(range1_idx) == 0 or \
                len(range2_idx) == 0:

                is_fragile = False

                return is_fragile

            at_least_one_edge = \
                find_edge_pairwise(range1_idx,
                                   range2_idx,
                                   conn_graph)

            is_fragile = not at_least_one_edge

            return is_fragile


    def initialize_event_series_dict(self) -> None:

        def parse_range_expansion(ch1_set,
                                  ch2_set,
                                  ana_smap,
                                  ana_conn_graph):
            # disambiguate the source region
            # (this updates from_region_idx inside stoch map)
            self.disambiguate_range_expansion(ana_smap,
                                              it_idx)

            from_region_idx = ana_smap.from_region_idx
            to_region_idx = ana_smap.region_gained_idx

            # check if dispersal is split-relevant
            is_split_relevant = \
                self.is_range_expansion_split_relevant(
                    from_region_idx,
                    to_region_idx,
                    ch1_set,
                    ch2_set
                )
            ana_smap.split_relevant = is_split_relevant

            # check if regions involved in dispersals
            # belong to the same communication class
            # (i.e., there is a path of connectivity
            # edges between them)
            are_connected = \
                ana_conn_graph.are_connected(from_region_idx,
                                             to_region_idx,
                                             ana_smap.from_state_bit_patt)
            ana_smap.within_comm_class = \
                are_connected

            # check if dispersal was over barrier
            # (note that a dispersal may be over a
            # barrier, but still happen within a
            # communicating class)
            over_barrier = \
                True if (from_region_idx, to_region_idx) \
                        not in ana_conn_graph.edge_set else \
                    False

            if not self._directed_edges:
                over_barrier = \
                    over_barrier or \
                    (to_region_idx, from_region_idx) not in \
                    ana_conn_graph.edge_set

            ana_smap.over_barrier = over_barrier

            # see if pre-dispersal range is fragile,
            # and annotate ana_smap accordingly
            expanding_range_bp = ana_smap.from_state_bit_patt
            expanding_range_set = \
                set([idx for idx, b \
                     in enumerate(expanding_range_bp) if b == "1"])
            expanding_range_fragile = \
                self.is_fragile_wrt_split(ch1_set,
                                          ch2_set,
                                          ana_conn_graph,
                                          expanding_range_set)

            ana_smap.previously_fragile_wrt_split = \
                expanding_range_fragile

            # see if post-dispersal range is stable,
            # and annotate ana_smap accordingly
            expanded_range_bp = ana_smap.to_state_bit_patt
            expanded_range_set = \
                set([idx for idx, b \
                     in enumerate(expanded_range_bp) if b == "1"])
            expanded_range_fragile = \
                self.is_fragile_wrt_split(ch1_set,
                                          ch2_set,
                                          ana_conn_graph,
                                          expanded_range_set)

            ana_smap.stabilized_range_wrt_split = \
                not expanded_range_fragile


        def parse_range_contraction(ch1_set,
                                    ch2_set,
                                    ana_smap,
                                    ana_conn_graph,
                                    ana_conn_graph_prev_time_slice):
            # check if range before contraction
            # was unstable (and annotate ana_smap
            # accordingly)
            pre_contr_range_bp = ana_smap.from_state_bit_patt
            pre_contr_range_set = \
                set([idx for idx, b \
                     in enumerate(pre_contr_range_bp) if b == "1"])
            pre_contr_range_fragile = \
                self.is_fragile_wrt_split(ch1_set,
                                          ch2_set,
                                          ana_conn_graph_prev_time_slice,
                                          pre_contr_range_set)
            ana_smap.previously_fragile_wrt_split = \
                pre_contr_range_fragile

            # check if range after contraction
            # was unstable (and annotate ana_smap
            # accordingly)
            post_contr_range_bp = ana_smap.to_state_bit_patt
            post_contr_range_set = \
                set([idx for idx, b \
                     in enumerate(post_contr_range_bp) if b == "1"])
            post_contr_range_fragile = \
                self.is_fragile_wrt_split(ch1_set,
                                          ch2_set,
                                          ana_conn_graph,
                                          post_contr_range_set)
            ana_smap.fragile_after_contr_wrt_split = \
                post_contr_range_fragile

            # NOTE:
            # if range is stable before the contraction and unstable
            # after, then this is 'speciation by extinction'

        def do_anagenetic_smap(ch1_set, ch2_set, ana_smap):
            # gathering initial information to annotate
            # anagenetic stochastic map

            # if no features, assume one epoch
            ana_smap_time_slice_idx = 0
            if self._geofeat_query is not None:
                ana_smap_time_slice_idx = \
                    self._geofeat_query. \
                        find_epoch_idx(ana_smap.age)

            # connectivity graph in the time slice (i.e., epoch)
            # where the stochastic map is situated
            ana_conn_graph = \
                self._geofeat_query. \
                    conn_graph_list[ana_smap_time_slice_idx]

            if ana_smap.map_type == "contraction":
                # to determine if the contraction may have led
                # to speciation, we need to look at the connectivity
                # graph in the previous time slice as well
                ana_conn_graph_prev_time_slice = \
                    self._geofeat_query. \
                        conn_graph_list[ana_smap_time_slice_idx - 1]

                parse_range_contraction(
                    ch1_set,
                    ch2_set,
                    ana_smap,
                    ana_conn_graph,
                    ana_conn_graph_prev_time_slice
                )

            if ana_smap.map_type == "expansion":
                parse_range_expansion(
                    ch1_set,
                    ch2_set,
                    ana_smap,
                    ana_conn_graph)

        def recursively_populate_event_series_dict(nd: dp.Node,
                                                   it_idx: int):
            """Populate self._event_series_dict recursively.

            Disambiguation of range expansion (dispersal) events
            happens here.
            """

            # do current node
            smap_on_tree = self._smap_collection.stoch_maps_tree_dict[it_idx]
            nd_name = nd.label
            event_series = EvolRelevantEventSeries()
            parent_nd_name = ""
            parent_anagenetic_smaps_list: ty.List[pjsmap.StochMap] = list()

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

            # only internal nodes (speciation!)
            if nd_name in smap_on_tree.clado_stoch_maps_dict:
                clado_smap = smap_on_tree.clado_stoch_maps_dict[nd_name]

                # if no features, assume one epoch
                clado_smap_time_slice_idx = 0
                if self._geofeat_query is not None:
                    clado_smap_time_slice_idx = \
                        self._geofeat_query. \
                            find_epoch_idx(clado_smap.age)

                clado_conn_graph = \
                    self._geofeat_query. \
                        conn_graph_list[clado_smap_time_slice_idx]

                # child 1
                ch1_bp = clado_smap.to_state_bit_patt
                # child 2 (will be None if no range split at speciation)
                ch2_bp = clado_smap.to_state2_bit_patt
                # get sets of region indices for the two mutually
                # exclusive ranges
                ch1_set = \
                    set([idx for idx, b in enumerate(ch1_bp) if b == "1"])

                ch2_set = ch1_set
                if ch2_bp is not None:
                    ch2_set = \
                        set([idx for idx, b in enumerate(ch2_bp) if b == "1"])

                # if not ch1_set.isdisjoint(ch2_set):
                #     exit(("Error during truncation of event series. "
                #           "At range split, the two mutually exclusive ranges shared"
                #           " regions. Exiting..."))

                # only speciation events with range splitting
                # if clado_smap.to_state2_bit_patt is not None:
                # cladogenetic
                #
                # annotate cladogenetic stochastic map depending on
                # whether splitting range is or not fragile at splitting
                # moment
                splitting_range_is_fragile = \
                    self.is_fragile_wrt_split(ch1_set,
                                              ch2_set,
                                              clado_conn_graph)
                clado_smap.splitting_range_fragile = \
                    splitting_range_is_fragile

                # anagenetic
                # NOTE: assumes stochastic maps are sorted in chronological order!!!
                # (old first, young later)
                if nd_name in smap_on_tree.anag_stoch_maps_dict:
                    #############################################
                    # Further annotation of anagenetic changes  #
                    # for this node, after we have re-annotated #
                    # anagenetic changes inherited from parent  #
                    #############################################
                    anagenetic_smaps_list = \
                        smap_on_tree.anag_stoch_maps_dict[nd_name]

                    for ana_smap in anagenetic_smaps_list:
                        do_anagenetic_smap(ch1_set, ch2_set, ana_smap)

                    # annotating parent events as well
                    for parent_ana_smap in event_series.event_list:
                        if isinstance(parent_ana_smap,
                                      (pjsmap.RangeExpansion,
                                       pjsmap.RangeContraction)):
                            do_anagenetic_smap(ch1_set,
                                               ch2_set,
                                               parent_ana_smap)

                    event_series.add_events(anagenetic_smaps_list)

                    # cladogenetic (has to be the last one in the event series)
                    event_series.add_events(clado_smap)

                    # debugging
                    # if it_idx == 3 and nd_name in ("nd5", "nd7"):
                    #     print("\nDoing", nd_name)
                    #     for ev in event_series.event_list:
                    #         print(ev)
                    #     print("\n")

                # update or create value in dictionary
                self._event_series_dict[nd_name][it_idx] = event_series

            # recur
            for ch_nd in nd.child_nodes():
                recursively_populate_event_series_dict(ch_nd, it_idx)

        # iterating over each MCMC iteration when stochastic maps were logged
        for it_idx, smap in self._smap_collection.stoch_maps_tree_dict.items():
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


    def add_paleogeo_event_series_dict(self,
                                       geo_cond_name: str,
                                       directed: bool = False,
                                       user_young2old: bool = True) -> None:
        """Update self.event_series_dict with paleogeographic events.

        Only barriers that matter for a specific range split are
        considered.

        Args:
            geo_cond_name (str):
            directed (bool): Flag specifying if connectivity cares
                about edge direction or not. Defaults to 'False'.
            user_young2old (bool): Defaults to 'False'.
        """

        conn_gained_epoch_start_ages_mat = self._geofeat_query.\
            geo_cond_change_times_dict[geo_cond_name]
        conn_lost_epoch_start_ages_mat = self._geofeat_query.\
            geo_cond_change_back_times_dict[geo_cond_name]

        for nd_label, it_event_series_dict in self._event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series will be an empty dictionary if range
                # did not split at node -- we do not care about those!
                # (even if subtending branch has maps!)
                if isinstance(event_series, EvolRelevantEventSeries):
                    event_list = event_series.event_list

                    # root/origin will have an empty smap_list
                    if len(event_list) > 0:
                        # first we find out the mutually exclusive sets of regions
                        range_split_smap = event_list[-1]  # last event is range split
                        range_split_age = range_split_smap.age
                        ch1_bp = range_split_smap.to_state_bit_patt
                        ch2_bp = range_split_smap.to_state2_bit_patt
                        ch1_set = set([idx for idx, b in enumerate(ch1_bp) if b == "1"])

                        ch2_set = ch1_set
                        if ch2_bp is not None:
                            ch2_set = set([idx for idx, b in enumerate(ch2_bp) if b == "1"])

                        # if not ch1_set.isdisjoint(ch2_set):
                        #     exit(("Error during truncation of event series. "
                        #           "At range split, the two mutually exclusive ranges shared"
                        #           " regions. Exiting..."))

                        # update event_list, going over connectivity graph for
                        # each epoch, and we want old to young (conn_graph_list must also
                        # be old -> young!)
                        idx_sequence = list()
                        if user_young2old:
                            idx_sequence = range(self._geofeat_query.n_time_slices)

                        else:
                            idx_sequence = reversed(range(self._geofeat_query.n_time_slices))

                        for epoch_idx in idx_sequence:
                            epoch_start_age = self._geofeat_query.feat_coll. \
                                epoch_age_start_list_old2young[epoch_idx]

                            # if epoch starts after the range split at speciation,
                            # we do not care about it
                            if epoch_start_age < range_split_age:
                                break

                            this_epoch_conn_graph = \
                                self._geofeat_query.conn_graph_list[epoch_idx]
                            edge_set = this_epoch_conn_graph.edge_set

                            # check if for at least one pair of nodes in
                            # the splitting (mutually exclusive) set of
                            # regions is connected in the connectivity
                            # graph
                            for from_region_idx in ch1_set:
                                for to_region_idx in ch2_set:
                                    e = (from_region_idx, to_region_idx)
                                    rev_e = (to_region_idx, from_region_idx)

                                    # if range is stable
                                    if e in edge_set or \
                                            (not directed and rev_e in edge_set):

                                        # and if barrier disappeared in this epoch
                                        # (note that we are looking at change_epoch instead
                                        # of change_back because of the behavior of
                                        # _populate_conn_graph_list() in GeoFeatureQuery,
                                        # which assumes 1 means connectivity!)
                                        if epoch_start_age in \
                                                conn_gained_epoch_start_ages_mat\
                                                        [from_region_idx][to_region_idx]:
                                            bd = pjfio.BarrierDisappearance(
                                                self._n_char,
                                                epoch_start_age,
                                                from_region_idx,
                                                to_region_idx
                                            )

                                            # add barrier disappearance as event
                                            # according to its age (see __lt__ of
                                            # EvolRelevantEvent)
                                            insort(event_list, bd)

                                    # if range has at least one barrier
                                    if e not in edge_set or \
                                            (not directed and rev_e not in edge_set):

                                        # and if barrier appeared in this epoch
                                        if epoch_start_age in \
                                                conn_lost_epoch_start_ages_mat\
                                                        [from_region_idx][to_region_idx]:
                                            ba = pjfio.BarrierAppearance(
                                                self._n_char,
                                                epoch_start_age,
                                                from_region_idx,
                                                to_region_idx
                                            )

                                            # add barrier appearance as event
                                            # according to its age (see __lt__ of
                                            # EvolRelevantEvent)
                                            insort(event_list, ba)

    def update_trunc_event_series(self) -> None:
        """Populate _truncated_event_series_dict

        This method visits all event series (for all internal nodes
        and all iterations) and makes a deep copy of it. Then it
        truncates each event series at the first split-relevant
        over-barrier dispersal event after the last split-relevant
        stable range.
        """

        # truncated event series data structure
        #
        # it    splitting_nd_idx    event_abbrev_csv_list   hyp
        #  1    5   b+_0_1,d_0_1_o_r,b-_0_1 vicariance
        #  2    5   ... founder-event

        for nd_label, it_event_series_dict in self._event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():

                # all events, paleogeographic and biogeographic
                #
                # root or origin may not have event series,
                # in which case event_series will be empty dictionary
                if isinstance(event_series, dict):
                    break

                event_list = event_series.event_list

                # we make deep copy because the annotation we carry out below
                # (.split_relevant) will be different for each node in the tree,
                # for the same dispersal event
                trunc_event_list = [copy.deepcopy(event) for event in event_list]

                self.event_series_dict[nd_label][it_idx].\
                    trunc_event_list = trunc_event_list

                # now we update truncated list
                if len(event_list) > 1:
                    for ev_idx, ev in enumerate(event_list):
                        # barrier disappearance is always split-relevant,
                        # but we must check if the actual range was made
                        # stable by the barrier disappearance (connectivity
                        # changes are irrelevant if a species is occupying
                        # a single region, for example)
                        #
                        # we need to know the range right at the moment of
                        # re-connection due to a barrier disappearing
                        know_range_before = False
                        range_before_bp = str()

                        if ev_idx == 0:
                            # there are no stoch maps older than this paleogeo
                            # event, so we need to find the to_range_bit_patt
                            # of a younger stoch map
                            for next_ev_idx in range((ev_idx + 1), len(event_list)):
                                next_ev = event_list[next_ev_idx]

                                # could be anagenetic, cladogenetic or no change!
                                if isinstance(next_ev,
                                              pjsmap.StochMap):
                                    range_before_bp = next_ev.from_state_bit_patt
                                    know_range_before = True

                        elif ev_idx > 0:
                            # first we try to find a stochastic map happening
                            # before the paleogeographic event
                            for rev_ev_idx in range((ev_idx-1), 0, -1):
                                prev_ev = event_list[rev_ev_idx]

                                if isinstance(prev_ev, (pjsmap.RangeExpansion,
                                                        pjsmap.RangeContraction,
                                                        pjsmap.NoChange)):
                                    range_before_bp = prev_ev.to_state_bit_patt
                                    know_range_before = True
                                    break

                                elif isinstance(prev_ev, pjsmap.RangeSplitOrBirth):
                                    # could not find range before if we have
                                    # a speciation event right before paleogeo
                                    # event -- can't know if to use left or right
                                    # 'to_range' bit pattern
                                    break

                            # if we could not find a stochastic map before the
                            # paleogeo event, it must be because we're deep in the
                            # tree, or the previous event was a range split and the
                            # 'to' range is ambiguous (no idea if left or right)
                            #
                            # we must now try to see what the 'from' range is of
                            # the next stochastic map
                            if not know_range_before:
                                for next_ev_idx in range((ev_idx+1), len(event_list)):
                                    next_ev = event_list[next_ev_idx]

                                    # could be anagenetic, cladogenetic or no change!
                                    if isinstance(next_ev,
                                                  pjsmap.StochMap):
                                        range_before_bp = next_ev.from_state_bit_patt
                                        know_range_before = True

                        # Truncation #
                        truncate_here = False

                        ###########################################
                        # Truncation case 1: a barrier disappears #
                        ###########################################

                        # if a barrier disappears, we need to know that
                        # the reconnected regions
                        if isinstance(ev, pjfio.BarrierDisappearance):
                            assert(range_before_bp)
                            occ_regions_idxs = \
                                set([idx for idx, b in \
                                     enumerate(range_before_bp) if b == "1"])

                            reg1_idx = ev.from_node_idx
                            reg2_idx = ev.to_node_idx

                            # if the 2 regions reconnected by the
                            # barrier disappearance (which is split-relevant,
                            # so no need to check!) are occupied, then
                            # we have a stable range (even if the range is not
                            # fully connected!)
                            if reg1_idx in occ_regions_idxs and \
                                reg2_idx in occ_regions_idxs:
                                truncate_here = True

                        ############################################
                        # Truncation case 2: a split-irrelevant    #
                        # dispersal happens, but it reconnects the #
                        # splitting components                     #
                        ############################################

                        elif isinstance(ev, pjsmap.RangeExpansion):
                            if ev.stabilized_range_wrt_split:
                                truncate_here = True

                        # actually truncate!
                        if truncate_here:
                            event_series.trunc_event_list = \
                                trunc_event_list[ev_idx:]

    def annotate_event_series_hypothesis(self) -> None:
        """Look at truncated event series and annotate hypothesis."""

        for nd_label, it_event_series_dict in self._event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # root or origin may not have event series,
                # in which case event_series will be empty dictionary
                if isinstance(event_series, dict):
                    break

                trunc_event_list = event_series.trunc_event_list
                n_events = len(trunc_event_list)
                destabilizing_dispersal = False
                destabilizing_extinction = False
                first_destabilizing_event = ""

                # no events means speciation within region
                if n_events == 0:
                    event_series.supported_hyp = Hypothesis.AMBIGUOUS

                    # debugging
                    # print("annotating ambiguous for ", nd_label, "it", it_idx)

                # at least one event means speciation between regions
                else:
                    clado_smap = trunc_event_list[-1]

                    # if splitting range is stable, nothing special happened
                    if not clado_smap.splitting_range_fragile:
                        event_series.supported_hyp = Hypothesis.AMBIGUOUS

                        # debugging
                        # print("range is stable, annotating ambiguous for ", nd_label, "it", it_idx)

                        # next iteration, same node!
                        continue

                    elif n_events > 1:
                        #########################################
                        # Scanning events                       #
                        # (further classification comes later!) #
                        #########################################
                        for ev in trunc_event_list:
                            # dispersal (range expansion)
                            if isinstance(ev, pjsmap.RangeExpansion):
                                if ev.split_relevant and \
                                        not ev.stabilized_range_wrt_split:
                                    # we have seen a destabilizing dispersal!
                                    destabilizing_dispersal = True

                                    # if destabilizing dispersal is the
                                    # first destabilizing event and range
                                    # was not previously unstable
                                    if not ev.previously_fragile_wrt_split and \
                                            first_destabilizing_event == "":
                                        first_destabilizing_event = "d"

                            # extinction (range contraction)
                            elif isinstance(ev, pjsmap.RangeContraction):
                                if ev.fragile_after_contr_wrt_split:
                                    destabilizing_extinction = True

                                    if not ev.previously_fragile_wrt_split and \
                                            first_destabilizing_event == "":
                                        first_destabilizing_event = "e"

                    ################
                    # Classifying! #
                    ################

                    if first_destabilizing_event == "e" and not \
                            destabilizing_dispersal:
                        event_series.supported_hyp = Hypothesis.SPECIATION_BY_EXT

                        # debugging
                        # print("annotating sp by ext for ", nd_label, "it", it_idx)

                    elif first_destabilizing_event == "d" and not \
                            destabilizing_extinction:
                        event_series.supported_hyp = Hypothesis.FOUNDER_EVENT

                        # debugging
                        # print("annotating founder event for ", nd_label, "it", it_idx)

                    # range was already unstable when either destabilizing dispersal
                    # or destabilizing extinction happened
                    elif first_destabilizing_event == "":
                        if destabilizing_dispersal and not \
                                destabilizing_extinction:
                            event_series.supported_hyp = Hypothesis.AMBIGUOUS

                            # debugging
                            # print("annotating vic for ", nd_label, "it", it_idx)

                        elif not destabilizing_dispersal and \
                                destabilizing_extinction:
                            event_series.supported_hyp = Hypothesis.AMBIGUOUS

                        elif not destabilizing_dispersal and not \
                                destabilizing_extinction:
                            event_series.supported_hyp = Hypothesis.VICARIANCE

                    else:
                        # debugging
                        print("first_destabilizing_event", first_destabilizing_event)
                        print("destabilizing_dispersal", destabilizing_dispersal)
                        print("destabilizing_extinction", destabilizing_extinction)
                        exit("found truncated event series that I do not know how to parse")


    @property
    def event_series_dict(self) -> ty.Dict[str, ty.Dict[int, EvolRelevantEventSeries]]:
        return self._event_series_dict

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

    node_states_file_path = "examples/trees_maps_files/geosse_dummy_tree3_tip_states.tsv"
    stoch_maps_file_path = "examples/trees_maps_files/geosse_dummy_tree3_maps.tsv"
    # node_states_file_path = "examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv"
    # stoch_maps_file_path = "examples/trees_maps_files/geosse_dummy_tree2_maps.tsv"
    smap_coll = \
        pjsmap.StochMapsOnTreeCollection(stoch_maps_file_path,
                                         ann_tr_list,
                                         state2bit_lookup,
                                         node_states_file_path=node_states_file_path,
                                         stoch_map_attr_name="state")

    # param_log_dir = "examples/feature_files/two_regions_feature_set_event_series"
    param_log_dir = "examples/feature_files/four_regions_feature_set_event_series"
    frs = FromRegionSampler(
        n_chars,
        param_log_dir,
        "epoch_age_",
        "_rel_rates",
        "m_d"
    )

    # feature_summary_fp = "examples/feature_files/two_regions_feature_set_event_series/feature_summary.csv"
    # age_summary_fp = "examples/feature_files/two_regions_feature_set_event_series/age_summary.csv"
    feature_summary_fp = "examples/feature_files/four_regions_feature_set_event_series/feature_summary.csv"
    age_summary_fp = "examples/feature_files/four_regions_feature_set_event_series/age_summary.csv"

    fc = pjfio.GeoFeatureCollection(
        feature_summary_fp,
        age_summary_fp=age_summary_fp)

    fq = pjfio.GeoFeatureQuery(fc)

    requirement_fn = \
        pjfio.GeoFeatureQuery.cb_feature_equals_value_is_connected(fc, 1, feat_id=1)

    # all members, including graph, are populated here!
    fq.populate_geo_cond_member_dicts("land_bridge", requirement_fn)

    est = \
        EvolRelevantEventSeriesTabulator(
            ann_tr_list,
            smap_coll,
            fq,
            from_region_sampler=frs
        )

    # looking at things
    # in epoch_age_2_rel_rates, make iteration 1's A -> C 666 to force scenario (ii)
    it_to_look_at = [3]
    for nd_label, it_event_series_dict in est.event_series_dict.items():
        for it_idx, event_series in it_event_series_dict.items():
            if isinstance(event_series, EvolRelevantEventSeries):
                if it_idx in it_to_look_at:
                    print(nd_label)
                    # for ev in event_series.trunc_event_list:
                    for ev in event_series.event_list:
                        # print(ev)
                        print(ev.short_str())
                    print(event_series)
        print("\n")