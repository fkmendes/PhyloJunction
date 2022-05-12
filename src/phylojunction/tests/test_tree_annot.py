import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/tests/)
import unittest
from dendropy import Tree, Node, Taxon

# pj imports
from data.tree import AnnotatedTree

class TestAnnotateTree(unittest.TestCase): 
    
    @classmethod
    def setUpClass(cls):
        rootedgeless_tr_str = "((nd6:2.5,nd7:1.0)nd3:1.0,(nd4:1.0,nd5:2.5)nd2:1.0)root:0.0;"
        rootedge_tr_str = "(((nd6:3.0,nd7:1.0)nd3:1.0,(nd4:1.0,nd5:3.0)nd2:1.0)root:2.0)origin:0.0;"
        rootedgeless_no_spn_living_tr_str = "(root:1.0)origin:0.0;"
        rootedgeless_no_spn_dead_tr_str = "(root:0.5)origin:0.0;"

        tr_root = Tree.get(data=rootedgeless_tr_str, schema="newick")
        tr_origin = Tree.get(data=rootedge_tr_str, schema="newick")
        tr_no_spn = Tree.get(data=rootedgeless_no_spn_living_tr_str, schema="newick")
        tr_no_spn_dead = Tree.get(data=rootedgeless_no_spn_dead_tr_str, schema="newick")

        origin_node1 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node1.alive = False
        root_node1 = Node(taxon=Taxon(label="root"), label="root", edge_length=1.0)
        root_node1.alive = True
        origin_node1.add_child(root_node1)
        tr_no_spn_built = Tree(seed_node=origin_node1)

        origin_node2 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node2.alive = False
        root_node2 = Node(taxon=Taxon(label="root"), label="root", edge_length=0.5)
        root_node2.alive = True
        origin_node2.add_child(root_node2)
        tr_no_spn_dead_built = Tree(seed_node=origin_node2)

        total_state_count = 1

        cls.tree_root = AnnotatedTree(tr_root, total_state_count, start_at_origin=False, epsilon=1e-12)
        cls.tree_origin = AnnotatedTree(tr_origin, total_state_count, start_at_origin=True, epsilon=1e-12)
        cls.tree_no_spn = AnnotatedTree(tr_no_spn, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)
        cls.tree_no_spn_dead = AnnotatedTree(tr_no_spn_dead, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)
        cls.tree_no_spn_built = AnnotatedTree(tr_no_spn_built, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)
        cls.tree_no_spn_dead_built = AnnotatedTree(tr_no_spn_dead_built, total_state_count, start_at_origin=True, max_age=1.0, epsilon=1e-12)

        # for testing the parsing of states into AnnotatedTree member dictionary
        origin_node3 = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node3.alive = False
        origin_node3.state = 0
        root_node3 = Node(taxon=Taxon(label="root"), label="root", edge_length=0.5)
        root_node3.alive = False
        root_node3.state = 0
        origin_node3.add_child(root_node3)
        child_left = Node(taxon=Taxon(label="nd1"), label="nd1", edge_length=1.0)
        child_left.alive = True
        child_left.state = 0
        root_node3.add_child(child_left)
        child_right = Node(taxon=Taxon(label="nd2"), label="nd2", edge_length=1.0)
        child_right.alive = True
        child_right.state = 1
        root_node3.add_child(child_right)
        bif_tr_two_states = Tree(seed_node=origin_node3)

        total_state_count = 2

        cls.bifurcating_tree_two_states = AnnotatedTree(bif_tr_two_states, total_state_count, start_at_origin=True, epsilon=1e-12)

    def test_root_origin_age_rootedge_length(self):
        """
        Test annotation of root/origin age and root edge
        """
        
        # tree starting at root
        self.assertEqual(self.tree_root.root_edge_length, 0.0, "No root edge: should be 0.0.")
        self.assertIsNone(self.tree_root.origin_age, "Origin age should be 'none'.")
        self.assertEqual(self.tree_root.root_age, 3.5, "Root age should be 3.5.")

        # tree starting at origin
        self.assertEqual(self.tree_origin.root_edge_length, 2.0, "Root edge should have length 2.0.")
        self.assertEqual(self.tree_origin.origin_age, 6.0, "Origin age should be 6.0.")
        self.assertEqual(self.tree_origin.root_age, 4.0, "Root age should be 4.0.")
        self.assertEqual(self.tree_no_spn.root_edge_length, 1.0, "Root edge should have length 1.0.")
        self.assertEqual(self.tree_no_spn_dead.root_edge_length, 0.5, "Root edge should have length 0.5.")
        self.assertEqual(self.tree_no_spn_built.root_edge_length, 1.0, "Root edge should have length 1.0.")
        self.assertEqual(self.tree_no_spn_dead_built.root_edge_length, 0.5, "Root edge should have length 0.5.")

        # tree starting at origin
        self.assertFalse(self.tree_no_spn.tree_died, msg="Tree made all the way to maximum age (1.0), so it should not have been marked as dead.")
        self.assertTrue(self.tree_no_spn_dead.tree_died, msg="Tree did not make all the way to maximum age (1.0), so it should have been marked as dead.")
        self.assertFalse(self.tree_no_spn_built.tree_died, msg="Tree made all the way to maximum age (1.0), so it should not have been marked as dead.")
        self.assertTrue(self.tree_no_spn_dead_built.tree_died, msg="Tree did not make all the way to maximum age (1.0), so it should have been marked as dead.")

    def test_node_counting(self):
        """
        Test counting of extant and extinct nodes
        """
        
        # tree starting at root
        # extant nodes
        self.assertEqual(self.tree_root.n_extant_obs_nodes, 2, "Count of observable extant nodes should be 2.")

        # extant nodes label check
        self.assertSequenceEqual(self.tree_root.extant_obs_nodes_labels, ("nd6", "nd5"), "Labels should be (\"nd6\", \"nd5\").")

        # extinct nodes
        self.assertEqual(self.tree_root.n_extinct_obs_nodes, 2, "Count of observable extinct nodes should be 2.")
        
        # extinct nodes label check
        self.assertSequenceEqual(self.tree_root.extinct_obs_nodes_labels, ("nd7", "nd4"), "Labels should be (\"nd7\", \"nd4\").")

        # ----------------------- #

        # tree starting at origin
        # extant nodes
        self.assertEqual(self.tree_origin.n_extant_obs_nodes, 2, "Count of observable extant nodes should be 2.")
        self.assertEqual(self.tree_no_spn.n_extant_obs_nodes, 1, "Count of observable extant nodes should be 1.")
        self.assertEqual(self.tree_no_spn_built.n_extant_obs_nodes, 1, "Count of observable extant nodes should be 1.")
        self.assertEqual(self.tree_no_spn_dead.n_extant_obs_nodes, 0, "Count of observable extant nodes should be 0.")
        self.assertEqual(self.tree_no_spn_dead_built.n_extant_obs_nodes, 0, "Count of observable extant nodes should be 0.")
        
        # extant nodes label check
        self.assertSequenceEqual(self.tree_origin.extant_obs_nodes_labels, ("nd6", "nd5"), "Labels should be (\"nd6\", \"nd5\").")
        self.assertSequenceEqual(self.tree_no_spn.extant_obs_nodes_labels, ("root",), "Labels should be (\"root\").") # must have comma
        self.assertSequenceEqual(self.tree_no_spn_built.extant_obs_nodes_labels, ("root",), "Labels should be (\"root\").") # must have comma

        # extinct nodes
        self.assertEqual(self.tree_origin.n_extinct_obs_nodes, 2, "Count of observable extinct nodes should be 2.")
        self.assertEqual(self.tree_no_spn_dead.n_extinct_obs_nodes, 1, "Count of observable extinct nodes should be 1.")
        self.assertEqual(self.tree_no_spn_dead_built.n_extinct_obs_nodes, 1, "Count of observable extinct nodes should be 1.")
        
        # extinct node label check
        self.assertSequenceEqual(self.tree_origin.extinct_obs_nodes_labels, ("nd7", "nd4"), "Labels should be (\"nd7\", \"nd4\").")
        self.assertSequenceEqual(self.tree_no_spn_dead.extinct_obs_nodes_labels, ("root",), "Labels should be (\"root\").") # must have comma
        self.assertSequenceEqual(self.tree_no_spn_dead_built.extinct_obs_nodes_labels, ("root",), "Labels should be (\"root\").") # must have comma

    # def test_node_attr_dict_population(self):
    #     """
    #     Test counting of extant and extinct nodes
    #     """

    #     self.bifurcating_tree_two_states.populate_nd_attr_dict(["state"])
    #     self.bifurcating_tree_two_states.get_gcf(node_attr="state")

if __name__ == '__main__':
    # can be called from tests/
    # $ python3 test_tree_annot.py
    # 
    # can also be called from phylojunction/
    # $ python3 tests/test_tree_annot.py
    # or
    # $ python3 -m tests.test_tree_annot
    # or, for a specific test
    # $ python3 -m unittest tests.test_tree_annot.TestAnnotateTree.test_root_origin_age_rootedge_length
    #
    # can also be called from VS Code, if open folder is phylojuction/

    unittest.main()