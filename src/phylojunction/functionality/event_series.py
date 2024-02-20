import typing as ty
import enum
import copy

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.tree as pjt
import phylojunction.functionality.evol_event as pjev
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.biogeo as pjbio
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class Hypothesis(enum.Enum):
    VICARIANCE = 0
    FOUNDER_EVENT = 1


class EvolRelevantEventSeries:
    """Series of evolution-relevant events.

    This class has member methods that interrogate a series of event,
    truncating it depending on what the user specifies.

    """

    _event_list: ty.List[pjev.EvolRelevantEvent]
    _n_events: int
    _series_type: str
    _it_idx: int

    # not initialized at instantiation, but populated
    # by external functions
    _trunc_event_list: ty.List[pjev.EvolRelevantEvent]
    _trunc_n_events: int
    _supported_hyp: Hypothesis

    def __init__(self,
                 event_list: ty.List[pjev.EvolRelevantEvent],
                 series_type: str,
                 it_idx: int) -> None:

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
    _event_series_dict: ty.Dict[str, ty.List[EvolRelevantEventSeries]]
    _hyp_support_dict: ty.Dict[Hypothesis, int]

    def __init__(self,
                 ann_tr_list: ty.List[pjt.AnnotatedTree],
                 smap_collection: pjsmap.StochMapsOnTreeCollection) -> None:

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
        self.update_event_char_member()
        print("  Finished going through all event series and annotating it.")

        # for each truncated event series, update its supported_hyp
        # member depending on the nature of the truncated series
        self.update_event_series_hyp_member()
        print("  Finished classifying event series as supporting each hypothesis.")

        # side-effect:
        # initializes self._hyp_support_dict
        self.tabulate_hyp_support()
        print("  Finished tabulating hypothesis support.")

    def initialize_event_series_dict(self) -> None:
        # (1) traverse the tree
        # (2) collect stoch maps for b/w-region speciation event
        # (3) initialize event series
        # (4)

        # iterating over each MCMC iteration when stochastic maps were logged
        for it_idx in smap_coll.sorted_it_idxs:
            ann_tr = self._ann_tr_list[it_idx]

            # traversing the tree
            for nd in self.ann_tr.tree.preorder_traversal():

                # collect full series

                # truncate series
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