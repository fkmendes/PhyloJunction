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
            ("examples/feature_files/two_regions_feature_set/"
             "feature_summary.csv")
        two_region_age_summary_fp = \
            ("examples/feature_files/two_regions_feature_set/"
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
            pjgeo.GeoFeatureQuery.cb_feature_equals_value(cls.two_region_geo_coll,
                                                          0,
                                                          feat_name="cb_1")

        cls.four_cb_is_1_requirement_fn = \
            pjgeo.GeoFeatureQuery.cb_feature_equals_value(cls.four_region_geo_coll,
                                                          0,
                                                          feat_name="cb_1")

    def test_read_features_epoch_ages_two_regions(self) -> None:
        """
        Test that GeoFeatureCollection (specifically its initialization
        method _read_age_summary()) is populating the right ages into
        the class' appropriate members. This test uses
        feature_files/two_regions_feature_set.
        """

        self.assertEqual(self.two_region_geo_coll.epoch_age_end_list_young2old,
                         [0.0, 1.5, 4.25])
        self.assertEqual(self.two_region_geo_coll.epoch_age_end_list_old2young,
                         [4.25, 1.5, 0.0])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_young2old,
                         [1.5, 4.25, -math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_old2young,
                         [-math.inf, 4.25, 1.5])
        self.assertEqual(self.two_region_geo_coll.epoch_mid_age_list_young2old,
                         [0.75, 2.875, -math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_mid_age_list_old2young,
                         [-math.inf, 2.875, 0.75])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_young2old,
                         [1.5, 4.25, -math.inf])
        self.assertEqual(self.two_region_geo_coll.epoch_age_start_list_old2young,
                         [-math.inf, 4.25, 1.5])

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
            [10.0, 20.0, -math.inf])
        self.assertEqual(
            self.four_region_geo_coll.epoch_age_start_list_old2young,
            [-math.inf, 20.0, 10.0])
        self.assertEqual(
            self.four_region_geo_coll.epoch_mid_age_list_young2old,
            [5.0, 15.0, -math.inf])
        self.assertEqual(
            self.four_region_geo_coll.epoch_mid_age_list_old2young,
            [-math.inf, 15.0, 5.0])

    def test_read_features_names_two_regions(self) -> None:
        """
        Test that GeoFeaturesCollection (specifically its initialization
        method _read_feat_summary_init_feats is populating the right
        region names and indices into the class' appropriate members.
        This test uses
        feature_files/two_regions_feature_set.
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

        exp1 = ('Feature (cb_1) | between-categorical 1 | Epoch 1 '
                '| 2 regions\n[[0 0]\n [0 0]]')
        exp2 = ('Feature (cb_1) | between-categorical 1 | Epoch 2 '
                '| 2 regions\n[[0 1]\n [1 0]]')
        exp3 = ('Feature (cb_1) | between-categorical 1 | Epoch 3 '
                '| 2 regions\n[[0 0]\n [0 0]]')

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
        # two regions
        self.two_geo_query.populate_geo_cond_member_dicts(
            "land_bridge",
            self.two_cb_is_1_requirement_fn)

        geo_cond_bit_patterns_list = list()
        for k in self.two_geo_query.geo_cond_bit_dict["land_bridge"]:
            # * unpacks list
            geo_cond_bit_patterns_list.append(" ".join(k))

        # A -> A, A -> B, B -> A, B -> B
        # each bit is for an epoch
        geo_cond_bit_patterns = " ".join(i for i in geo_cond_bit_patterns_list)
        self.assertEqual(geo_cond_bit_patterns, "111 101 101 111")

        self.assertEqual(
        self.two_geo_query.geo_oldest_cond_bit_dict,
            {'land_bridge': [['1', '1'], ['1', '1']]}
        )

        self.assertEqual(
            self.two_geo_query.get_geo_condition_change_times("land_bridge"),
            [[[], [1.5]], [[1.5], []]]
        )

        self.assertEqual(
            self.two_geo_query.get_geo_condition_change_back_times("land_bridge"),
            [[[], [4.25]], [[4.25], []]]
        )


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