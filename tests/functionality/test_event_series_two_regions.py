import unittest

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.event_series as pjev
import phylojunction.functionality.feature_io as pjfio

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestEventSeriesTwoRegions(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        n_chars = 2

        state2bit_lookup = \
            pjbio.State2BitLookup(n_chars,
                                  2,
                                  geosse=True)

        # number of ranges
        n_states = state2bit_lookup.n_states

        tr_fp = "examples/trees_maps_files/geosse_dummy_tree2.tre"
        ann_tr_list = [pjr.read_nwk_tree_str(tr_fp,
                                             "read_tree",
                                             node_names_attribute="index",
                                             n_states=n_states,
                                             in_file=True)]

        node_states_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree2_tip_states.tsv"
        stoch_maps_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree2_maps.tsv"

        smap_coll = \
            pjsmap.StochMapsOnTreeCollection(
                stoch_maps_file_path,
                ann_tr_list,
                state2bit_lookup,
                node_states_file_path=node_states_file_path,
                stoch_map_attr_name="state")

        param_log_dir = \
            "examples/feature_files/two_regions_feature_set_event_series"

        frs = pjev.FromRegionSampler(
            n_chars,
            param_log_dir,
            "epoch_age_",
            "_rel_rates",
            "m_d"
        )

        feature_summary_fp = \
            ("examples/feature_files/two_regions_feature_set_event_series/"
             "feature_summary.csv")
        age_summary_fp = \
            ("examples/feature_files/two_regions_feature_set_event_series/"
             "/age_summary.csv")

        fc = pjfio.GeoFeatureCollection(
            feature_summary_fp,
            age_summary_fp=age_summary_fp)

        fq = pjfio.GeoFeatureQuery(fc)

        requirement_fn = \
            pjfio.GeoFeatureQuery.\
                cb_feature_equals_value_is_connected(
                fc, 1, feat_id=1)

        fq.populate_geo_cond_member_dicts(
            "land_bridge",
            requirement_fn)

        cls.est = \
            pjev.EvolRelevantEventSeriesTabulator(
                ann_tr_list,
                smap_coll,
                fq,
                from_region_sampler=frs
            )

    def test_event_list_it1(self) -> None:
        """Test two-region, 5-taxon tree event list."""

        nd6_exp = ("d(4.75)_sr|nb|wc|st(10)>st(11)"
                   "b+(4.25)"
                   "s(4.0)_un(11>10_01)")
        nd8_exp = ("b+(4.25)"
                   "d(3.25)_sr|ob|oc|st(10)>un(11)"
                   "s(3.0)_un(11>10_01)")
        nd7_exp = ("b+(4.25)"
                   "d(3.25)_sr|ob|oc|st(10)>un(11)"
                   "s(3.0)_un(11>10_01)"
                   "d(2.75)_sr|ob|oc|st(01)>un(11)"
                   "b-(1.5)"
                   "s(1.0)_st(11>10_01)")

        nd6_obs = str()
        nd8_obs = str()
        nd7_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            if nd_label == "nd6":
                                nd6_obs += ev.short_str()

                            elif nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd8":
                                nd8_obs += ev.short_str()

        self.assertEqual(nd6_obs, nd6_exp)
        self.assertEqual(nd7_obs, nd7_exp)
        self.assertEqual(nd8_obs, nd8_exp)

    def test_trunc_event_list_it1(self) -> None:
        """Test two-region, 5-taxon tree truncated event list."""

        nd6_exp = ("d(4.75)_sr|nb|wc|st(10)>st(11)"
                   "b+(4.25)"
                   "s(4.0)_un(11>10_01)")
        nd8_exp = ("b+(4.25)"
                   "d(3.25)_sr|ob|oc|st(10)>un(11)"
                   "s(3.0)_un(11>10_01)")
        nd7_exp = "b-(1.5)s(1.0)_st(11>10_01)"

        nd6_obs = str()
        nd8_obs = str()
        nd7_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            if nd_label == "nd6":
                                nd6_obs += ev.short_str()

                            elif nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd8":
                                nd8_obs += ev.short_str()

    def test_hyp(self) -> None:
        """Test two-region, 5-taxon tree hypothesis."""

        nd6_exp = str("Vicariance")
        nd8_exp = str("Founder event")
        nd7_exp = str("Ambiguous")

        nd6_obs = str()
        nd8_obs = str()
        nd7_obs = str()

        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if nd_label == "nd6":
                        nd6_obs = str(event_series.supported_hyp)

                    elif nd_label == "nd7":
                        nd7_obs = str(event_series.supported_hyp)

                    elif nd_label == "nd8":
                        nd8_obs = str(event_series.supported_hyp)

        self.assertEqual(nd6_obs, nd6_exp)
        self.assertEqual(nd8_obs, nd8_exp)
        self.assertEqual(nd7_obs, nd7_exp)

if __name__ == "__main__":
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_event_series_two_regions.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_event_series_two_regions
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_event_series_two_regions.TestEventSeriesTwoRegions.test_it1

    unittest.main()