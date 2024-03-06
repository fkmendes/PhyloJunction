import unittest

# pj imports
import phylojunction.functionality.event_series as pjes

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestRegionSampler(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.frs = pjes.FromRegionSampler(
            2,
            "examples/feature_files/two_regions_feature_set_event_series",
            "epoch_",
            "_rel_rates",
            "m_d"
        )

    def test_region_sampler_reading_methods(self) -> None:
        """Test parameters in .log file are read correctly."""

        # there should be 3 dictionaries, one for each epoch
        self.assertEqual(len(self.frs.time_slice_dict_list), 3)

        # epoch with index 1 (user-specified, should be youngest)
        self.assertEqual(
            self.frs.time_slice_dict_list[0],
            {1: [[1.0, 2.0], [0.1, 1.0]], 10: [[1.0, 4.0], [0.5, 1.0]]}
        )

        # epoch with index 2 (user-specified, should be youngest)
        self.assertEqual(
            self.frs.time_slice_dict_list[1],
            {1: [[1.0, 0.1], [0.1, 1.0]], 10: [[1.0, 0.2], [0.2, 1.0]]}
        )

        # epoch with index 3 (user-specified, should be youngest)
        self.assertEqual(
            self.frs.time_slice_dict_list[2],
            {1: [[1.0, 0.1], [2.1, 1.0]], 10: [[1.0, 0.2], [2.2, 1.0]]}
        )


    def test_region_sampler_sample_region(self) -> None:
        """Test regions are sampled correctly."""

        sampled_regions = list()
        for i in range(20000):
            # note that with two regions only, there is nothing to disambiguate!
            # this test is artificial in that we are forcing the potential 'from region' to be 0 and 1!
            # it should nevertheless still serve its purpose
            sampled_regions.append(
                self.frs.sample_from_region_idx(1, 0, [0,1], 1)
            )

        n_region_idx0_it0_epoch0 = sum(1 for i in sampled_regions if i == 0)
        n_region_idx1_it0_epoch0 = sum(1 for i in sampled_regions if i == 1)
        self.assertAlmostEqual(
            float(n_region_idx0_it0_epoch0 / n_region_idx1_it0_epoch0),
            2.0 / 1.0,
            delta=0.05
        )

        sampled_regions = list()
        for i in range(20000):
            # note that with two regions only, there is nothing to disambiguate!
            # this test is artificial in that we are forcing the potential 'from region' to be 0 and 1!
            # it should nevertheless still serve its purpose
            sampled_regions.append(
                self.frs.sample_from_region_idx(1, 1, [0, 1], 1)
            )

        n_region_idx0_it0_epoch1 = sum(1 for i in sampled_regions if i == 0)
        n_region_idx1_it0_epoch1 = sum(1 for i in sampled_regions if i == 1)
        self.assertAlmostEqual(
            float(n_region_idx0_it0_epoch1 / n_region_idx1_it0_epoch1),
            0.1 / 1.0,
            delta=0.05
        )

        sampled_regions = list()
        for i in range(20000):
            # note that with two regions only, there is nothing to disambiguate!
            # this test is artificial in that we are forcing the potential 'from region' to be 0 and 1!
            # it should nevertheless still serve its purpose
            sampled_regions.append(
                self.frs.sample_from_region_idx(10, 0, [0, 1], 1)
            )

        n_region_idx0_it10_epoch1 = sum(1 for i in sampled_regions if i == 0)
        n_region_idx1_it10_epoch1 = sum(1 for i in sampled_regions if i == 1)
        self.assertAlmostEqual(
            float(n_region_idx0_it10_epoch1 / n_region_idx1_it10_epoch1),
            4.0 / 1.0,
            delta=0.1
        )


if __name__ == '__main__':
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_region_sampler.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_region_sampler
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_region_sampler.TestRegionSampler.test_region_sampler_reading_methods

    unittest.main()