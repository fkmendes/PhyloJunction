import unittest

# pj imports
import phylojunction.functionality.feature_io as pjgeo

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestConnGraph(unittest.TestCase):

    def test_graph1(self) -> None:
        """Test multiple comm. classes are correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0}, {1}, {2}, {3}, {4}, {5}, {6}
        # (no edges!)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0}, {1}, {2}, {3}, {4}, {5}, {6}])
        self.assertFalse(g.are_connected(0, 1))
        self.assertFalse(g.are_connected(2, 3))
        self.assertFalse(g.are_connected(4, 5))
        self.assertFalse(g.are_connected(0, 2))
        self.assertFalse(g.are_connected(0, 4))
        self.assertFalse(g.are_connected(2, 4))
        self.assertFalse(g.are_connected(6, 0))
        self.assertFalse(g.are_connected(6, 3))
        self.assertFalse(g.are_connected(6, 4))

    def test_graph2(self) -> None:
        """Test multiple comm. classes are correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0,1}, {2,3}, {4,5}, {6}
        g.add_edge(0, 1, is_directed=True)
        g.add_edge(2, 3, is_directed=True)
        g.add_edge(4, 5, is_directed=True)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0, 1}, {2, 3}, {4, 5}, {6}])
        self.assertTrue(g.are_connected(0, 1))
        self.assertTrue(g.are_connected(2, 3))
        self.assertTrue(g.are_connected(4, 5))
        self.assertFalse(g.are_connected(0, 2))
        self.assertFalse(g.are_connected(0, 4))
        self.assertFalse(g.are_connected(2, 4))
        self.assertFalse(g.are_connected(6, 0))
        self.assertFalse(g.are_connected(6, 3))
        self.assertFalse(g.are_connected(6, 4))

    def test_graph3(self) -> None:
        """Test multiple single comm. classe is correctly built."""

        n_regions = 7

        g = pjgeo.GeoGraph(n_regions)

        # comm classes: {0,1,2,3,4,5,6}
        g.add_edge(0, 1, is_directed=True)
        g.add_edge(2, 3, is_directed=True)
        g.add_edge(4, 5, is_directed=True)
        g.add_edge(6, 0, is_directed=True)
        g.add_edge(6, 2, is_directed=True)
        g.add_edge(6, 4, is_directed=True)

        g.populate_comm_class_members()

        # debugging
        # print(g.comm_class_set_list)

        self.assertEqual(g.comm_class_set_list,
                         [{0, 1, 2, 3, 4, 5, 6}])
        self.assertTrue(g.are_connected(0, 1))
        self.assertTrue(g.are_connected(2, 3))
        self.assertTrue(g.are_connected(4, 5))
        self.assertTrue(g.are_connected(0, 2))
        self.assertTrue(g.are_connected(0, 4))
        self.assertTrue(g.are_connected(2, 4))
        self.assertTrue(g.are_connected(6, 0))
        self.assertTrue(g.are_connected(6, 3))
        self.assertTrue(g.are_connected(6, 4))


if __name__ == '__main__':
    # From PhyloJunction/
    #
    # $ python3 tests/functionality/test_conn_graph.py
    #
    # or
    #
    # $ python3 -m tests.functionality.test_conn_graph
    #
    # or
    #
    # $ python3 -m unittest tests.functionality.test_conn_graph.TestConnGraph.test_graph1

    unittest.main()