import typing as ty
import enum
import copy
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

    def __init__(self,
                 ann_tr_list: ty.List[pjt.AnnotatedTree],
                 smap_collection: pjsmap.StochMapsOnTreeCollection) -> None:

        self._event_series_dict = pjh.autovivify(2)

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
            for nd_label, it_event_series_dict in self._event_series_dict.items():
                print(nd_label)
                for it_idx, event_series in it_event_series_dict.items():
                    for ev in event_series.event_list:
                        print(ev)
                print("\n")


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