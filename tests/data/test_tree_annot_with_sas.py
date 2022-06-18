import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa

class TestAnnotateTreeWithSAs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        ##############################################################################################
        # Tree 1: origin + SA then survive until with max age;                                       #
        # Note that the way "execute_sample_ancestor" is implemented when it happens before a proper #
        # birth event (giving rise to the root), the origin instead undergoes a "fake birth" event,  #
        # where a 'dummy node' is made its child, and two children are added to the dummy node:      # 
        #                                                                                            #
        # (1) SA,                                                                                    #
        # (2) A node between the origin and the root (the "brosc" node)                              #
        #                                                                                            #
        # Then (2) stays alive and potentially undergoes another birth event when the root would     #
        # be created. In this Tree 1, however, (2) stays alive and does not undergo a birth event    #
        ##############################################################################################       
        origin_node = Node(taxon=Taxon(label="origin"), label="origin", edge_length=0.0)
        origin_node.alive = False
        origin_node.is_sa = False
        origin_node.is_sa_dummy_parent = False
        origin_node.is_sa_lineage = False

        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.alive = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        
        origin_node.add_child(dummy_node)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.alive = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        brosc_node = Node(taxon=Taxon(label="brosc"), label="brosc", edge_length=1.0)
        brosc_node.alive = True
        brosc_node.is_sa = False
        brosc_node.is_sa_dummy_parent = False
        brosc_node.is_sa_lineage = True
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(brosc_node)

        # no spn: no root
        tr_sa_no_spn_built = Tree(seed_node=origin_node)

        print("tr_sa_no_spn_built.seed_age = " + str(tr_sa_no_spn_built.max_distance_from_root()))

        total_state_count = 1

        sa_global_time = 2.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor("sa1", "bw_origin_root", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "bw_origin_root": [sa] }
        
        cls.tree_sa_no_spn_built = pjtr.AnnotatedTree(
            tr_sa_no_spn_built,
            total_state_count,
            start_at_origin=True,
            max_age=2.0,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

    def test_node_counting(self):
        """
        Test counting of extant, extinct and SA nodes
        """
        
        # Tree 1: origin + SA then survive until with max age
        self.assertEqual(self.tree_sa_no_spn_built.n_extant_terminal_nodes, 1, "Count of terminal extant nodes should be 1.")
        self.assertEqual(self.tree_sa_no_spn_built.n_extinct_terminal_nodes, 0, "Count of terminal extinct nodes should be 0.")
        self.assertEqual(self.tree_sa_no_spn_built.n_sa, 1, "Count of sampled ancestor nodes should be 1.")

        # Tree 2: origin + SA + SA then survive until with max age

        # Tree 3: origin + SA then extinction before max age

        # Tree 4: origin + SA + SA then extinction before max age

if __name__ == "__main__":
    # Assuming you opened the PhyloJunction/ (repo root) folder
    # on VSCode and that you want to run this as a standalone script,
    # i.e., "Run Without Debugging", you will need to configure your
    # launch.json to have:
    # 
    # "env": {"PYTHONPATH": "${workspaceRoot}/src/phylojunction/"}
    # 
    # and your settings.json to have:
    #   
    # "python.analysis.extraPaths": [ "${workspaceFolder}/src/phylojunction/" ]
    # 
    # If you want to run this as a standalone from PhyloJunction/
    # on the terminal, remember to add "src/phylojunction" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3 tests/data/test_tree_annot_with_sas.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_annot_with_sas
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_annot_with_sas.TestAnnotateTreeWithSAs.test_node_counting

    unittest.main()