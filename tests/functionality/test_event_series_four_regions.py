import unittest

# pj imports
import phylojunction.readwrite.pj_read as pjr
import phylojunction.functionality.biogeo as pjbio
import phylojunction.functionality.stoch_map as pjsmap
import phylojunction.functionality.event_series as pjev
import phylojunction.functionality.feature_io as pjfio

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestEventSeriesFourRegions(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        n_chars = 4

        state2bit_lookup = \
            pjbio.State2BitLookup(n_chars,
                                  2,
                                  geosse=True)

        # number of ranges
        n_states = state2bit_lookup.n_states

        tr_fp = "examples/trees_maps_files/geosse_dummy_tree3.tre"
        ann_tr_list = [pjr.read_nwk_tree_str(tr_fp,
                                             "read_tree",
                                             node_names_attribute="index",
                                             n_states=n_states,
                                             in_file=True)]

        node_states_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree3_tip_states.tsv"
        stoch_maps_file_path = \
            "examples/trees_maps_files/geosse_dummy_tree3_maps.tsv"

        smap_coll = \
            pjsmap.StochMapsOnTreeCollection(
                stoch_maps_file_path,
                ann_tr_list,
                state2bit_lookup,
                node_states_file_path=node_states_file_path,
                stoch_map_attr_name="state")

        param_log_dir = \
            "examples/feature_files/four_regions_feature_set_event_series"

        frs = pjev.FromRegionSampler(
            n_chars,
            param_log_dir,
            "epoch_age_",
            "_rel_rates",
            "m_d"
        )

        feature_summary_fp = \
            ("examples/feature_files/four_regions_feature_set_event_series/"
             "feature_summary.csv")
        age_summary_fp = \
            ("examples/feature_files/four_regions_feature_set_event_series/"
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

    def test_event_list_4reg_it1(self) -> None:
        """Test four-region, 4-taxon tree event list, iteration 1."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)")
        nd5_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)"
                   "b+(3.0)_reg(0|3)_destab(1001)"
                   "d(2.5)_si|nb|wc|un(1001)>un(1101)"
                   "b+(2.0)_reg(1|2)_N/A(1101)"
                   "d(0.75)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)"
                   "b+(3.0)_reg(0|3)_destab(1001)"
                   "b+(2.0)_reg(1|2)_N/A(1001)"
                   "d(1.5)_sr|ob|oc|un(1001)>un(1011)"
                   "d(0.75)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_event_list_4reg_it2(self) -> None:
        """Test four-region, 4-taxon tree event list, iteration 2."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_N/A(1101)"
                   "s(3.25)_st(1101>1101_1101)")
        nd5_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_destab(1101)"
                   "s(3.25)_st(1101>1101_1101)"
                   "b+(3.0)_reg(0|3)_destab(1101)"
                   "e(2.5)_|st(1101)>un(1100)"
                   "b+(2.0)_reg(1|2)_N/A(1100)"
                   "d(1.5)_sr|ob|oc|st(1100)>un(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_destab(1101)"
                   "s(3.25)_st(1101>1101_1101)"
                   "b+(3.0)_reg(0|3)_destab(1101)"
                   "e(2.5)_|st(1101)>st(1001)"
                   "b+(2.0)_reg(1|2)_N/A(1001)"
                   "d(1.5)_si|nb|wc|un(1001)>un(1011)"
                   "d(0.75)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [2]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_event_list_4reg_it3(self) -> None:
        """Test four-region, 4-taxon tree event list, iteration 3."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1000)"
                   "d(4.5)_si|ob|wc|st(1000)>st(1010)"
                   "b+(4.0)_reg(1|3)_N/A(1010)"
                   "s(3.25)_st(1010>1010_1010)")
        nd5_exp = ("b+(5.0)_reg(0|2)_N/A(1000)"
                   "d(4.5)_sr|ob|wc|st(1000)>un(1010)"
                   "b+(4.0)_reg(1|3)_N/A(1010)"
                   "s(3.25)_st(1010>1010_1010)"
                   "b+(3.0)_reg(0|3)_N/A(1010)"
                   "d(2.5)_si|nb|wc|un(1010)>st(1110)"
                   "b+(2.0)_reg(1|2)_destab(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1000)"
                   "d(4.5)_sr|ob|wc|st(1000)>un(1010)"
                   "b+(4.0)_reg(1|3)_N/A(1010)"
                   "s(3.25)_st(1010>1010_1010)"
                   "b+(3.0)_reg(0|3)_N/A(1010)"
                   "b+(2.0)_reg(1|2)_N/A(1010)"
                   "d(1.5)_si|nb|wc|un(1010)>un(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [3]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_trunc_event_list_4reg_it1(self) -> None:
        """Test four-region, 4-taxon tree truncated event list, iteration 1."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)")
        nd5_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)"
                   "b+(3.0)_reg(0|3)_destab(1001)"
                   "d(2.5)_si|nb|wc|un(1001)>un(1101)"
                   "b+(2.0)_reg(1|2)_N/A(1101)"
                   "d(0.75)_si|nb|wc|un(1101)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1001)"
                   "b+(4.0)_reg(1|3)_N/A(1001)"
                   "s(3.25)_st(1001>1001_1001)"
                   "b+(3.0)_reg(0|3)_destab(1001)"
                   "b+(2.0)_reg(1|2)_N/A(1001)"
                   "d(1.5)_sr|ob|oc|un(1001)>un(1011)"
                   "d(0.75)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [1]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_trunc_event_list_4reg_it2(self) -> None:
        """Test four-region, 4-taxon tree truncated event list, iteration 2."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_N/A(1101)"
                   "s(3.25)_st(1101>1101_1101)")
        nd5_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_destab(1101)"
                   "s(3.25)_st(1101>1101_1101)"
                   "b+(3.0)_reg(0|3)_destab(1101)"
                   "e(2.5)_|st(1101)>un(1100)"
                   "b+(2.0)_reg(1|2)_N/A(1100)"
                   "d(1.5)_sr|ob|oc|st(1100)>un(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1101)"
                   "b+(4.0)_reg(1|3)_destab(1101)"
                   "s(3.25)_st(1101>1101_1101)"
                   "b+(3.0)_reg(0|3)_destab(1101)"
                   "e(2.5)_|st(1101)>st(1001)"
                   "b+(2.0)_reg(1|2)_N/A(1001)"
                   "d(1.5)_si|nb|wc|un(1001)>un(1011)"
                   "d(0.75)_si|nb|wc|un(1011)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [2]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_trunc_event_list_4reg_it3(self) -> None:
        """Test four-region, 4-taxon tree truncated event list, iteration 3."""

        nd7_exp = ("b+(5.0)_reg(0|2)_N/A(1000)"
                   "d(4.5)_si|ob|wc|st(1000)>st(1010)"
                   "b+(4.0)_reg(1|3)_N/A(1010)"
                   "s(3.25)_st(1010>1010_1010)")
        nd5_exp = ("b+(2.0)_reg(1|2)_destab(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")
        nd6_exp = ("b+(5.0)_reg(0|2)_N/A(1000)"
                   "d(4.5)_sr|ob|wc|st(1000)>un(1010)"
                   "b+(4.0)_reg(1|3)_N/A(1010)"
                   "s(3.25)_st(1010>1010_1010)"
                   "b+(3.0)_reg(0|3)_N/A(1010)"
                   "b+(2.0)_reg(1|2)_N/A(1010)"
                   "d(1.5)_si|nb|wc|un(1010)>un(1110)"
                   "d(0.75)_si|nb|wc|un(1110)>un(1111)"
                   "s(0.25)_un(1111>1100_0011)")

        nd7_obs = str()
        nd5_obs = str()
        nd6_obs = str()

        it_to_look_at = [3]
        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if it_idx in it_to_look_at:
                        for ev in event_series.trunc_event_list:
                            if nd_label == "nd7":
                                nd7_obs += ev.short_str()

                            elif nd_label == "nd5":
                                nd5_obs += ev.short_str()

                            elif nd_label == "nd6":
                                nd6_obs += ev.short_str()

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_hyp_ann_all_its_4reg(self) -> None:
        """Test four-region, 4-taxon tree hypothesis annotation, all iterations."""

        # iteration 1, 2 and 3
        nd7_exp = {1: "Within-region speciation", 2: "Within-region speciation", 3: "Within-region speciation"}
        nd5_exp = {1: "Vicariance", 2: "Founder event", 3: "Vicariance"}
        nd6_exp = {1: "Ambiguous", 2: "Speciation by extinction", 3: "Founder event"}

        nd7_obs = dict((i, "") for i in [1, 2, 3])
        nd5_obs = dict((i, "") for i in [1, 2, 3])
        nd6_obs = dict((i, "") for i in [1, 2, 3])

        for nd_label, it_event_series_dict in self.est.event_series_dict.items():
            for it_idx, event_series in it_event_series_dict.items():
                # event_series of root will be an empty dictionary, need
                # to check
                if isinstance(event_series, pjev.EvolRelevantEventSeries):
                    if nd_label == "nd7":
                        nd7_obs[it_idx] = str(event_series.supported_hyp)

                    elif nd_label == "nd5":
                        nd5_obs[it_idx] = str(event_series.supported_hyp)

                    elif nd_label == "nd6":
                        nd6_obs[it_idx] = str(event_series.supported_hyp)

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)

    def test_node_count_supporting_hyp_dict_4reg(self) -> None:
        """Test four-region, 4-taxon tree, hypothesis node count."""

        wr_exp = 3
        amb_exp = 1
        fe_exp = 2
        sbe_exp = 1
        vic_exp = 2

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

        nd7_exp = \
            {'Vicariance': 0,
             'Founder event': 0,
             'Speciation by extinction': 0,
             'Ambiguous': 0,
             'Within-region speciation': 3}
        nd5_exp = \
            {'Vicariance': 2,
             'Founder event': 1,
             'Speciation by extinction': 0,
             'Ambiguous': 0,
             'Within-region speciation': 0}
        nd6_exp = \
            {'Vicariance': 0,
             'Founder event': 1,
             'Speciation by extinction': 1,
             'Ambiguous': 1,
             'Within-region speciation': 0}

        nd7_obs = self.est.hyp_support_by_node_dict["nd7"]
        nd5_obs = self.est.hyp_support_by_node_dict["nd5"]
        nd6_obs = self.est.hyp_support_by_node_dict["nd6"]

        self.assertEqual(nd7_exp, nd7_obs)
        self.assertEqual(nd5_exp, nd5_obs)
        self.assertEqual(nd6_exp, nd6_obs)


if __name__ == "__main__":
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_event_series_four_regions.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_event_series_four_regions
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_event_series_four_regions.TestEventSeriesFourRegions.test_hyp_support_by_node_dict_4reg

    unittest.main()