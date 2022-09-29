import unittest
from dendropy import Tree, Node, Taxon

# pj imports
import phylojunction.data.tree as pjtr
import phylojunction.data.sampled_ancestor as pjsa

class TestExtractReconstructedTree(unittest.TestCase):

    def test_extract_reconstructed_tree_one_extant_one_extinct_one_sa(self):
        """
        Test extraction of reconstructed tree from complete tree

        Tree has one extant taxon, one extinct taxon, one sampled ancestor

        To see it on icytree: ((sa1:0.0,sp1:1.0)dummy1:1.0,sp2:1.25)root:0.0;
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.alive = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens on who would have been the left child ("sp1")
        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.alive = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root node, goes extinct
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.25)
        extant_sp2.alive = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(dummy_node)
        root_node.add_child(extant_sp2)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.alive = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.5)
        extant_sp1.alive = True
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor("sa1", "sp1", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            total_state_count,
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        tr_rec = ann_tr.extract_reconstructed_tree()
        
        self.assertEqual(ann_tr.tree.as_string(schema="newick", suppress_annotations=True), "((sa1:0.0,sp1:0.5)dummy1_dummy1:1.0,sp2:1.25)root_root:0.0;\n")
        self.assertEqual(tr_rec.as_string(schema="newick", suppress_annotations=True), "((sa1:0.0,sp1:0.5)dummy1_dummy1:1.0)root_root:0.0;\n")

    
    def test_extract_reconstructed_tree_all_extinct_one_sa(self):
        """
        Test extraction of reconstructed tree from complete tree

        Tree has one extant taxon, two extinct taxa, one sampled ancestor

        To see it on icytree: ((sa1:0.0,sp1:1.0)dummy1:1.0,sp2:1.25)root:0.0;
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.alive = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node after ancestor sampling happens on who would have been the left child ("sp1")
        dummy_node = Node(taxon=Taxon(label="dummy1"), label="dummy1", edge_length=1.0)
        dummy_node.alive = False
        dummy_node.is_sa = False
        dummy_node.is_sa_dummy_parent = True
        dummy_node.is_sa_lineage = False

        # right child of root node, goes extinct
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.25)
        extant_sp2.alive = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(dummy_node)
        root_node.add_child(extant_sp2)
        
        # right child of dummy_node
        sa_node = Node(taxon=Taxon(label="sa1"), label="sa1", edge_length=0.0)
        sa_node.alive = False
        sa_node.is_sa = True
        sa_node.is_sa_dummy_parent = False
        sa_node.is_sa_lineage = False

        # left child of dummy node
        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=0.5)
        extant_sp1.alive = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False
        
        dummy_node.add_child(sa_node)
        dummy_node.add_child(extant_sp1)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(dummy_node.taxon)
        tr_complete.taxon_namespace.add_taxon(sa_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        sa_global_time = 1.0
        time_to_sa_lineage_node = 1.0
        sa = pjsa.SampledAncestor("sa1", "sp1", sa_global_time, time_to_lineage_node=time_to_sa_lineage_node)
        sa_lineage_dict = { "sp1": [sa] }
        
        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            total_state_count,
            start_at_origin=False,
            max_age=max_age,
            sa_lineage_dict=sa_lineage_dict,
            epsilon=1e-12)

        tr_rec = ann_tr.extract_reconstructed_tree()

        self.assertEqual(ann_tr.tree.as_string(schema="newick", suppress_annotations=True), "((sa1:0.0,sp1:0.5)dummy1_dummy1:1.0,sp2:1.25)root_root:0.0;\n")
        self.assertEqual(tr_rec.as_string(schema="newick", suppress_annotations=True), "((sa1:0.0)dummy1_dummy1:1.0)root_root:0.0;\n")


    def test_extract_reconstructed_tree_all_extinct_no_sa(self):
        """
        Test extraction of reconstructed tree from complete tree

        Tree has died
        """

        root_node = Node(taxon=Taxon(label="root"), label="root", edge_length=0.0)
        root_node.alive = False
        root_node.is_sa = False
        root_node.is_sa_dummy_parent = False
        root_node.is_sa_lineage = False

        # left child of root node
        extant_sp1 = Node(taxon=Taxon(label="sp1"), label="sp1", edge_length=1.0)
        extant_sp1.alive = False
        extant_sp1.is_sa = False
        extant_sp1.is_sa_dummy_parent = False
        extant_sp1.is_sa_lineage = False

        # right child of root node, goes extinct
        extant_sp2 = Node(taxon=Taxon(label="sp2"), label="sp2", edge_length=1.0)
        extant_sp2.alive = False
        extant_sp2.is_sa = False
        extant_sp2.is_sa_dummy_parent = False
        extant_sp2.is_sa_lineage = False

        root_node.add_child(extant_sp1)
        root_node.add_child(extant_sp2)

        tr_complete = Tree(seed_node=root_node)
        tr_complete.taxon_namespace.add_taxon(root_node.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp1.taxon)
        tr_complete.taxon_namespace.add_taxon(extant_sp2.taxon)

        total_state_count = 1

        max_age = 2.0

        ann_tr = pjtr.AnnotatedTree(
            tr_complete,
            total_state_count,
            start_at_origin=False,
            max_age=max_age,
            epsilon=1e-12)

        tr_rec = ann_tr.extract_reconstructed_tree()
        
        self.assertEqual(ann_tr.tree.as_string(schema="newick", suppress_annotations=True), "(sp1:1.0,sp2:1.0)root_root:0.0;\n")
        self.assertEqual(tr_rec.as_string(schema="newick", suppress_annotations=True), ";\n")


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
    # $ python3 tests/data/test_tree_extract_reconstructed.py
    # 
    # or
    #
    # $ python3 -m tests.data.test_tree_extract_reconstructed
    #
    # or 
    #
    # $ python3 -m unittest tests.data.test_tree_extract_reconstructed.TestExtractReconstructedTree.test_extract_reconstructed_tree_one_extant_one_extinct

    unittest.main()