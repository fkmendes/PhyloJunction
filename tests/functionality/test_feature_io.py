import unittest
import math

# pj imports
import phylojunction.functionality.feature_io as pjgeo

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestFeatureIO(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        two_region_feat_summary_fp = \
            ("examples/feature_files/two_regions_feature_set_event_series/"
             "feature_summary.csv")
        two_region_age_summary_fp = \
            ("examples/feature_files/two_regions_feature_set_event_series/"
             "age_summary.csv")
        four_region_feat_summary_fp = \
            ("examples/feature_files/four_regions_feature_set/"
             "feature_summary.csv")
        four_region_age_summary_fp = \
            ("examples/feature_files/four_regions_feature_set/"
             "age_summary.csv")

        cls.two_region_geo_coll = pjgeo.GeoFeatureCollection(
            two_region_feat_summary_fp,
            two_region_age_summary_fp)

        cls.four_region_geo_coll = pjgeo.GeoFeatureCollection(
            four_region_feat_summary_fp,
            four_region_age_summary_fp)

        cls.two_geo_query = pjgeo.GeoFeatureQuery(cls.two_region_geo_coll)
        cls.four_geo_query = pjgeo.GeoFeatureQuery(cls.four_region_geo_coll)

        cls.two_cb_is_1_requirement_fn = \
            pjgeo.GeoFeatureQuery.cb_feature_equals_value_is_connected(cls.two_region_geo_coll,
                                                                       1,
                                                                       feat_name="cb_1")

        cls.four_cb_is_1_requirement_fn = \
            pjgeo.GeoFeatureQuery.cb_feature_equals_value_is_connected(cls.four_region_geo_coll,
                                                                       0,
                                                                       feat_name="cb_1")

    def test_read_features_epoch_ages_two_regions(self) -> None:
        """
        Test that GeoFeatureCollection (specifically its initialization
        method _read_age_summary()) is populating the right ages into
        the class' appropriate members. This test uses
        feature_files/two_regions_feature_set_event_series.
        """

        self.assertEqual(self.two_region_geo_coll.epoch_age_end_list_young2old,
                         [0.0, 1.5, 4.25])
        self.assertEqual(self.two_region_geo_coll.epoch_age_end_list_old2young,
                         [4.25, 1.5, 0.0])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_young2old,
                         [1.5, 4.25, math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_old2young,
                         [math.inf, 4.25, 1.5])
        self.assertEqual(self.two_region_geo_coll.epoch_mid_age_list_young2old,
                         [0.75, 2.875, math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_mid_age_list_old2young,
                         [math.inf, 2.875, 0.75])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_young2old,
                         [1.5, 4.25, math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_old2young,
                         [math.inf, 4.25, 1.5])

    def test_read_features_epoch_ages_four_regions(self) -> None:
        """
        Test that GeoFeatureCollection (specifically its initialization
        method _read_age_summary()) is populating the right ages into
        the class' appropriate members. This test uses
        feature_files/four_regions_feature_set.
        """

        self.assertEqual(
            self.four_region_geo_coll.epoch_age_end_list_young2old,
            [0.0, 10.0, 20.0])
        self.assertEqual(
            self.four_region_geo_coll.epoch_age_end_list_old2young,
            [20.0, 10.0, 0.0])
        self.assertEqual(
            self.four_region_geo_coll.epoch_age_start_list_young2old,
            [10.0, 20.0, math.inf])
        self.assertEqual(
            self.four_region_geo_coll.epoch_age_start_list_old2young,
            [math.inf, 20.0, 10.0])
        self.assertEqual(
            self.four_region_geo_coll.epoch_mid_age_list_young2old,
            [5.0, 15.0, math.inf])
        self.assertEqual(
            self.four_region_geo_coll.epoch_mid_age_list_old2young,
            [math.inf, 15.0, 5.0])

    def test_epoch_find(self) -> None:
        """Test that auxiliary method finds the epoch a time belongs.

        This test also checks that the epoch index cache is built
        correctly (some times will be repeatedly queried if the
        tree never changes and its internal node ages remain
        constant).
        """

        # present epoch (2)
        age1 = 1.0
        epoch_idx1 = self.four_geo_query.find_epoch_idx(age1)

        # second youngest epoch (1)
        age2 = 11.0
        epoch_idx2 = self.four_geo_query.find_epoch_idx(age2)

        # oldest epoch (0)
        age3 = 21.0
        epoch_idx3 = self.four_geo_query.find_epoch_idx(age3)

        # present epoch (2)
        age4 = 10.0
        epoch_idx4 = self.four_geo_query.find_epoch_idx(age4)

        self.assertEqual([epoch_idx1, epoch_idx2, epoch_idx3, epoch_idx4],
                         [2, 1, 0, 2])

    def test_read_features_names_two_regions(self) -> None:
        """
        Test that GeoFeaturesCollection (specifically its initialization
        method _read_feat_summary_init_feats is populating the right
        region names and indices into the class' appropriate members.
        This test uses
        feature_files/two_regions_feature_set_event_series.
        """

        self.assertEqual(self.two_region_geo_coll.region_name_idx_dict,
                         {'reg1':0, 'reg2':1})

        self.assertEqual(self.two_region_geo_coll.region_idx_name_dict,
                         {0:'reg1', 1:'reg2'})

    def test_read_features_names_four_regions(self) -> None:
        """
        Test that GeoFeaturesCollection (specifically its initialization
        method _read_feat_summary_init_feats is populating the right
        region names and indices into the class' appropriate members.
        This test uses
        feature_files/four_regions_feature_set.
        """

        self.assertEqual(
            self.four_region_geo_coll.region_name_idx_dict,
            {'reg1':0, 'reg2':1, 'reg3':2, 'reg4':3})

        self.assertEqual(
            self.four_region_geo_coll.region_idx_name_dict,
            {0:'reg1', 1:'reg2', 2:'reg3', 3:'reg4'})

    def test_read_features_values_and_types_two_regions(self) -> None:
        """
        Test that GeoFeaturesCollection (specifically its initialization
        methods _init_feat_name_epochs_dict and
        _init_feat_type_rel_featid_epochs_dictis are populating the right
        info inside the class' appropriate members.
        """

        # from feature_summary.csv, epoch 1 is youngest, epoch 3 is oldest
        exp1 = ('Feature (cb_1) | between-categorical 1 | Epoch 1 '
                '| 2 regions\n[[0 1]\n [1 0]]')
        exp2 = ('Feature (cb_1) | between-categorical 1 | Epoch 2 '
                '| 2 regions\n[[0 0]\n [0 0]]')
        exp3 = ('Feature (cb_1) | between-categorical 1 | Epoch 3 '
                '| 2 regions\n[[0 1]\n [1 0]]')

        self.assertEqual(str(self.two_region_geo_coll.feat_name_epochs_dict['cb_1'][1]),
                         exp1)
        self.assertEqual(str(self.two_region_geo_coll.feat_name_epochs_dict['cb_1'][2]),
                         exp2)
        self.assertEqual(str(self.two_region_geo_coll.feat_name_epochs_dict['cb_1'][3]),
                         exp3)

        cb1_1_str = str(self.two_region_geo_coll.feat_type_rel_featid_epochs_dict \
                            [pjgeo.GeoFeatureType.CATEGORICAL] \
                            [pjgeo.GeoFeatureRelationship.BETWEEN][1][1])
        cb2_1_str = str(self.two_region_geo_coll.feat_type_rel_featid_epochs_dict \
                            [pjgeo.GeoFeatureType.CATEGORICAL] \
                            [pjgeo.GeoFeatureRelationship.BETWEEN][1][2])
        cb3_1_str = str(self.two_region_geo_coll.feat_type_rel_featid_epochs_dict \
                            [pjgeo.GeoFeatureType.CATEGORICAL] \
                            [pjgeo.GeoFeatureRelationship.BETWEEN][1][3])

        self.assertEqual(cb1_1_str, exp1)
        self.assertEqual(cb2_1_str, exp2)
        self.assertEqual(cb3_1_str, exp3)

    def test_geo_query_two_regions(self) -> None:
        """
        Test that GeoFeaturesQuery is capable of initializing its
        members correctly, namely:
            (i)   geo_cond_bit_dict
            (ii)  geo_condition_change_times
            (iii) geo_condition_change_back_times
            (iv)  conn_graph_list

        Using feature_files/two_regions_feature_set_event_series/
        """

        self.two_geo_query.populate_geo_cond_member_dicts(
            "land_bridge",
            self.two_cb_is_1_requirement_fn,
            True)

        geo_cond_bit_patterns_list = list()
        for k in self.two_geo_query.geo_cond_bit_dict["land_bridge"]:
            # * unpacks list
            geo_cond_bit_patterns_list.append(" ".join(k))

        # A -> A, A -> B, B -> A, B -> B
        # each bit is for an epoch, from young to old (see feature_summary.csv)
        geo_cond_bit_patterns = " ".join(i for i in geo_cond_bit_patterns_list)
        self.assertEqual(geo_cond_bit_patterns, "000 101 101 000")

        # self.assertEqual(
        # self.two_geo_query.geo_oldest_cond_bit_dict,
        #     {'land_bridge': [['0', '1'], ['0', '1']]}
        # )

        # A -> A, A -> B, B -> A, B -> B
        # when connectivity exists!
        self.assertEqual(
            self.two_geo_query.get_geo_condition_change_times("land_bridge"),
            [[[], [1.5]], [[1.5], []]]
        )

        # when connectivity does not exist!
        self.assertEqual(
            self.two_geo_query.get_geo_condition_change_back_times("land_bridge"),
            [[[], [4.25]], [[4.25], []]]
        )

        # epoch 1, idx = 0 (youngest, the index comes from feat_summary.csv)
        self.assertSetEqual(
            self.two_geo_query.conn_graph_list[0].edge_set,
            {(0, 1), (1, 0)}
        )

        # epoch 2, idx = 1 (no edges!)
        self.assertSetEqual(
            self.two_geo_query.conn_graph_list[1].edge_set,
            set()
        )

        # epoch 3, idx = 2 (oldest, the index comes from feat_summary.csv)
        self.assertSetEqual(
            self.two_geo_query.conn_graph_list[2].edge_set,
            {(0, 1), (1, 0)}
        )

        self.assertEqual(self.two_geo_query.get_comm_classes(0.1),
                         [{0, 1}])
        self.assertEqual(self.two_geo_query.get_comm_classes(1.0),
                         [{0, 1}])
        self.assertEqual(self.two_geo_query.get_comm_classes(1.6),
                         [{0}, {1}])
        self.assertEqual(self.two_geo_query.get_comm_classes(4.35),
                         [{0, 1}])

    def test_geo_query_four_regions(self) -> None:
        """
        Test that GeoFeaturesQuery is capable of initializing its
        members correctly, namely:
            (i)   geo_cond_bit_dict
            (ii)  geo_condition_change_times
            (iii) geo_condition_change_back_times
            (iv)  conn_graph_list

        Using feature_files/four_regions_feature_set/
        """

        self.four_geo_query.populate_geo_cond_member_dicts(
            "land_bridge",
            self.four_cb_is_1_requirement_fn,
            True)

        geo_cond_bit_patterns_list = list()
        for k in self.four_geo_query.geo_cond_bit_dict["land_bridge"]:
            # * unpacks list
            geo_cond_bit_patterns_list.append(" ".join(k))

        # each bit is an epoch, and it goes from young -> old
        # as specified by the time indices in the feature summary file
        self.assertEqual(geo_cond_bit_patterns_list[0],
                         '111 000 010 111')
        self.assertEqual(geo_cond_bit_patterns_list[1],
                         '000 111 010 111')
        self.assertEqual(geo_cond_bit_patterns_list[2],
                         '111 010 111 110')
        self.assertEqual(geo_cond_bit_patterns_list[3],
                         '111 111 011 111')

        # youngest epoch: 1-> 4 (111), 2 -> 4 (111), 3 -> 4 (110), so all connected
        # middle epoch:   1-> 4 (111), 2 -> 4 (111), 3 -> 4 (110), so all connected
        # oldeset epoch:  1-> 4 (111), 2 -> 4 (111), 4 -> 3 (011), so all connected
        # youngest epoch: 1-> 4 (111), 2 -> 4 (111), 3 -> 4 (110), so all connected
        self.assertEqual(self.four_geo_query.get_comm_classes(1.0),
                         [{0, 1, 2, 3}])
        self.assertEqual(self.four_geo_query.get_comm_classes(11.0),
                         [{0, 1, 2, 3}])
        self.assertEqual(self.four_geo_query.get_comm_classes(21.0),
                         [{0, 1, 2, 3}])
        self.assertEqual(self.four_geo_query.get_comm_classes(10.0),
                         [{0, 1, 2, 3}])


if __name__ == '__main__':
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_feature_io.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_feature_io
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_feature_io.TestFeatureIO.test_read_features

    unittest.main()