import os
import math
import enum
import copy
import typing as ty
from natsort import natsorted
import dendropy as dp

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.tree as pjt
import phylojunction.functionality.evol_event as pjev
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.biogeo as pjbio
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

    @property
    def sample_from_region_idx(self,
                               time_slice_idx: int,
                               potential_from_region_idx: ty.List[int],
                               to_region_idx: int) -> int:
        """Sample a region according to parameter values.

        Parameter values in a .log file can be normalized and serve
        as weights for sampling a 'from' region for dispersal
        stochastic maps.

        Args:
            time_slice_idx: Index of time slice.
            potential_from_region_idx: List of indices for regions
                from which dispersal (range expansion) may have
                happened.
            to_region_idx: Index for region to which dispersal (range
                expansion) happened.

        Returns:
            (int): A sampled region index.
        """

        param_values = list()
        for param_mat in self._param_value_dict[time_slice_idx]:
            for from_idx in potential_from_region_idx:
                param_values.append(param_mat[from_idx][to_region_idx])

        # 1. obtain weight list
        # 2. sample indices according to the weight list
        # 3. return sampled index


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

    # key of inner dict is the iteration index
    _event_series_dict: ty.Dict[str, ty.Dict[int, ty.List[EvolRelevantEventSeries]]]

    _hyp_support_dict: ty.Dict[Hypothesis, int]
    _from_region_sampler: FromRegionSampler

    def __init__(self,
                 ann_tr_list: ty.List[pjt.AnnotatedTree],
                 smap_collection: pjsmap.StochMapsOnTreeCollection,
                 from_region_sampler: ty.Optional[FromRegionSampler] = None) -> None:

        self._event_series_dict = pjh.autovivify(2)
        self._from_region_sampler = from_region_sampler

        self._ann_tr_list = ann_tr_list
        print("Read trees.")

        self._smap_collection = smap_collection
        print("Read stochastic maps.")

        print("Beginning event series parsing:")
        # side-effect:
        # initializes self._event_series_dict
        # initializes self._trunc_event_series_dict
        self.initialize_event_series_dict()
        print("  Finished building all event series and truncating them.")

        # for each event series, classify events as dispersal
        # or local extinction
        self.scan_update_events()

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

    def initialize_event_series_dict(self) -> None:

        def recursively_populate_event_series_dict(nd: dp.Node,
                                                   it_idx: int):
            """Populate self._event_series_dict recursively."""

            nd_name = nd.label
            event_series = EvolRelevantEventSeries()

            # retrieve list of events from parent node
            if nd.parent_node is not None:
                parent_nd = nd.parent_node
                parent_nd_name = parent_nd.label
                parent_event_series = \
                    self._event_series_dict[parent_nd_name][it_idx]

                # we will add events to those from the parent node
                # (has to be deep copy!)
                event_series = copy.deepcopy(parent_event_series)

            # do current node
            smap_on_tree = smap_coll.stoch_maps_tree_dict[it_idx]

            # anagenetic
            # NOTE: assumes stochastic maps are sorted in chronological order!!!
            # (old first, young later)
            if nd_name in smap_on_tree.anag_stoch_maps_dict:
                anagenetic_smaps_list = smap_on_tree.anag_stoch_maps_dict[nd_name]
                event_series.add_events(anagenetic_smaps_list)

            # TODO:
            # decode anagenetic smaps here!!! (i.e., find 'from' regions)

            # cladogenetic (has to be the last one in the event series)
            if nd_name in smap_on_tree.clado_stoch_maps_dict:
                clado_smap = smap_on_tree.clado_stoch_maps_dict[nd_name]
                event_series.add_events(clado_smap)

            # update or create value in dictionary
            self._event_series_dict[nd_name][it_idx] = event_series

            # recur
            for ch_nd in nd.child_nodes():
                recursively_populate_event_series_dict(ch_nd, it_idx)

        # iterating over each MCMC iteration when stochastic maps were logged
        for it_idx, smap in smap_coll.stoch_maps_tree_dict.items():
            ann_tr = smap.ann_tr
            root_nd = ann_tr.root_node

            # populate self._self._event_series_dict
            # key: node name
            # value: event series object
            recursively_populate_event_series_dict(root_nd,
                                                   it_idx)

            # debugging
            # for nd_label, it_event_series_dict in self._event_series_dict.items():
            #     print(nd_label)
            #     for it_idx, event_series in it_event_series_dict.items():
            #         for ev in event_series.event_list:
            #             print(ev)
            #     print("\n")

    def scan_update_events(self) -> None:
        n_regions = self._smap_collection.n_char
        print('n_regions', n_regions)

        for nd_label, it_event_series_dict \
                in self._event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                for ev in event_series.event_list:
                    # if range expansion (dispersal), we need to know which
                    # region the dispersal happened from
                    if ev.map_type == "expansion":
                        bp = ev.from_state_bit_patt
                        potential_from_region_idx = \
                            [idx for idx, b in enumerate(bp) if b == '1']

                        # if the class member exists, nothing needs to be done,
                        # otherwise, we need to determine it
                        if ev.from_region_idx == None:
                            # sample proportional to some scheme (e.g., proportional
                            # to FIG rate scalers, m_d)
                            if self._from_region_sampler != None:
                                sampled_idx = \
                                    self._from_region_sampler(potential_from_region_idx)

                            else:
                                # uniformly pick element of list
                                # random_pick(potential_from_region_idx)
                                pass

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

    ann_tr_list = [pjr.read_nwk_tree_str("examples/trees_maps_files/geosse_dummy_tree2.tre",
                                          "read_tree",
                                          node_names_attribute="index",
                                          n_states=3,
                                          in_file=True)]

    n_chars = 2
    state2bit_lookup = pjbio.State2BitLookup(n_chars, 2, geosse=True)

    smap_coll = \
        pjsmap.StochMapsOnTreeCollection("examples/trees_maps_files/geosse_dummy_tree2_maps.tsv",
                                         ann_tr_list,
                                         state2bit_lookup,
                                         node_states_file_path="examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv",
                                         stoch_map_attr_name="state")

    event_series_tabulator = \
        EvolRelevantEventSeriesTabulator(
            ann_tr_list,
            smap_coll
        )