import unittest

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.event_series as pjev
import phylojunction.functionality.feature_io as pjfio

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestEventSeriesABCDTricky(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        n_chars = 4

        state2bit_lookup = \
            pjbio.State2BitLookup(n_chars,
                                  2,
                                  geosse=True)

        # number of ranges
        n_states = state2bit_lookup.n_states

        tr_fp = "examples/trees_maps_files/geosse_dummy_tree4.tre"
        ann_tr_list = [pjr.read_nwk_tree_str(tr_fp,
                                             "read_tree",
                                             node_names_attribute="index",
                                             n_states=n_states,
                                             in_file=True)]

        node_states_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree4_tip_states.tsv"
        stoch_maps_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree4_maps.tsv"

        smap_coll = \
            pjsmap.StochMapsOnTreeCollection(
                stoch_maps_file_path,
                ann_tr_list,
                state2bit_lookup,
                node_states_file_path=node_states_file_path,
                stoch_map_attr_name="state")

        param_log_dir = \
            "examples/feature_files/feature_set_event_series_ABCD_tricky"

        frs = pjev.FromRegionSampler(
            n_chars,
            param_log_dir,
            "epoch_age_",
            "_rel_rates",
            "m_d"
        )

        feature_summary_fp = \
            ("examples/feature_files/feature_set_event_series_ABCD_tricky/"
             "feature_summary.csv")
        age_summary_fp = \
            ("examples/feature_files/feature_set_event_series_ABCD_tricky/"
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

    def test_event_list_4reg_tricky_it1(self) -> None:
        """Test four-region, 4-taxon tree tricky event list, iteration 1."""

        nd3_exp = ("b+(5.0)_reg(0|3)_N/A(1100)"
                   "d(3.5)_sr|ob|oc|st(1100)>un(1110)"
                   "b+(3.0)_reg(0|1)"
                   "b-(3.0)_reg(2|3)_N/A(1110)"
                   "d(2.5)_si|nb|wc|un(1110)>un(1111)"
                   "b-(2.0)_reg(0|1)"
                   "e(1.5)_un(1111)>un(1101)"
                   "d(0.5)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_event_list_4reg_tricky_it2(self) -> None:
        """Test four-region, 4-taxon tree tricky event list, iteration 2."""

        nd3_exp = ("b+(5.0)_reg(0|3)_destab(1101)"
                   "d(3.5)_sr|ob|oc|un(1101)>un(1111)"
                   "b+(3.0)_reg(0|1)"
                   "b-(3.0)_reg(2|3)"
                   "e(2.5)_un(1111)>un(1101)"
                   "b-(2.0)_reg(0|1)"
                   "d(1.5)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [2]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_event_list_4reg_tricky_it3(self) -> None:
        """Test four-region, 4-taxon tree tricky event list, iteration 3."""

        nd3_exp = ("b+(5.0)_reg(0|3)_destab(1001)"
                   "d(3.5)_sr|ob|oc|un(1001)>un(1011)"
                   "b+(3.0)_reg(0|1)_N/A(1011)"
                   "b-(3.0)_reg(2|3)"
                   "b-(2.0)_reg(0|1)_N/A(1011)"
                   "d(1.5)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [3]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_trunc_event_list_4reg_tricky_it1(self) -> None:
        """Test four-region, 4-taxon tree truncated tricky event list, iteration 1."""

        nd3_exp = ("b+(5.0)_reg(0|3)_N/A(1100)"
                   "d(3.5)_sr|ob|oc|st(1100)>un(1110)"
                   "b+(3.0)_reg(0|1)"
                   "b-(3.0)_reg(2|3)_N/A(1110)"
                   "d(2.5)_si|nb|wc|un(1110)>un(1111)"
                   "b-(2.0)_reg(0|1)"
                   "e(1.5)_un(1111)>un(1101)"
                   "d(0.5)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_trunc_event_list_4reg_tricky_it2(self) -> None:
        """Test four-region, 4-taxon tree truncated tricky event list, iteration 2."""

        nd3_exp = ("b+(5.0)_reg(0|3)_destab(1101)"
                   "d(3.5)_sr|ob|oc|un(1101)>un(1111)"
                   "b+(3.0)_reg(0|1)"
                   "b-(3.0)_reg(2|3)"
                   "e(2.5)_un(1111)>un(1101)"
                   "b-(2.0)_reg(0|1)"
                   "d(1.5)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [2]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_trunc_event_list_4reg_tricky_it3(self) -> None:
        """Test four-region, 4-taxon tree truncated tricky event list, iteration 3."""

        nd3_exp = ("b+(5.0)_reg(0|3)_destab(1001)"
                   "d(3.5)_sr|ob|oc|un(1001)>un(1011)"
                   "b+(3.0)_reg(0|1)_N/A(1011)"
                   "b-(3.0)_reg(2|3)"
                   "b-(2.0)_reg(0|1)_N/A(1011)"
                   "d(1.5)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd3_obs = str()

        it_to_look_at = [3]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            nd3_obs += ev.short_str()

        self.assertEqual(nd3_exp, nd3_obs)

    def test_hyp_ann_all_its_4reg_tricky(self) -> None:
        """Test four-region, 4-taxon tree hypothesis annotation, tricky, all iterations."""

        # iteration 1, 2 and 3
        nd3_exp = {1: "Founder event", 2: "Vicariance", 3: "Ambiguous"}

        nd3_obs = dict((i, "") for i in [1, 2, 3])

        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    nd3_obs[it_idx] = str(event_series.supported_hyp)

        self.assertEqual(nd3_exp, nd3_obs)

    def test_node_count_supporting_hyp_dict_4reg_tricky(self) -> None:
        """Test four-region, 4-taxon tree, tricky, hypothesis node count."""

        wr_exp = 0
        amb_exp = 1
        fe_exp = 1
        sbe_exp = 0
        vic_exp = 1

        hsd = self.est.node_count_supporting_hyp_dict

        wr_obs = 0
        amb_obs = 0
        fe_obs = 0
        sbe_obs = 0
        vic_obs = 0
        for it_idx, hyp_dict in hsd.items():
            for hyp, node_count in hyp_dict.items():
                if hyp == "Within-region speciation":
                    wr_obs += node_count

                elif hyp == "Ambiguous":
                    amb_obs += node_count

                elif hyp == "Founder event":
                    fe_obs += node_count

                elif hyp == "Speciation by extinction":
                    sbe_obs += node_count

                elif hyp == "Vicariance":
                    vic_obs += node_count

        self.assertEqual(wr_exp, wr_obs)
        self.assertEqual(amb_exp, amb_obs)
        self.assertEqual(fe_exp, fe_obs)
        self.assertEqual(sbe_exp, sbe_obs)
        self.assertEqual(vic_exp, vic_obs)

    def test_hyp_support_by_node_dict_4reg(self) -> None:
        """Test four-region, 4-taxon tree, node hypothesis."""

        nd3_exp = \
            {'Vicariance': 1,
             'Founder event': 1,
             'Speciation by extinction': 0,
             'Ambiguous': 1,
             'Within-region speciation': 0}

        nd3_obs = self.est.hyp_support_by_node_dict["nd3"]

        self.assertEqual(nd3_exp, nd3_obs)


if __name__ == "__main__":
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_event_series_four_regions_tricky.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_event_series_four_regions_tricky
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_event_series_four_regions_tricky.TestEventSeriesABCDTricky.test_hyp_support_by_node_dict_4reg

    unittest.main()