import unittest
import re

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmdbox.cmd_parse as cmd
import phylojunction.interface.cmdbox.cmd_parse_utils as cmdu

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class TestVarAssignment(unittest.TestCase):
    def test_var_assignment(self):
        """
        Test if a series of different variable assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """

        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        cmd_line1 = "a <- 1"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("a")
        
        self.assertEqual(1, pgm_obj.n_nodes)
        self.assertEqual(type(a_node_pgm.value), list)
        self.assertEqual(len(a_node_pgm.value), 1)
        self.assertEqual(a_node_pgm.value, ["1"])
        
        # --- #

        cmd_line2 = "b <- [1, 2, 3]"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line2)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        b_node_pgm = pgm_obj.get_node_pgm_by_name("b")

        self.assertEqual(2, pgm_obj.n_nodes)
        self.assertEqual(type(b_node_pgm.value), list)
        self.assertEqual(len(b_node_pgm.value), 3)
        self.assertEqual(b_node_pgm.value, ["1", "2", "3"])
        
        # --- #
        
        cmd_line3 = "a <- [1]"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line3)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line3)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("a")

        self.assertEqual(2, pgm_obj.n_nodes)
        self.assertEqual(type(a_node_pgm.value), list)
        self.assertEqual(len(a_node_pgm.value), 1)
        self.assertEqual(a_node_pgm.value, ["1"])

        # --- #
        
        cmd_line4 = "c <- b"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line4)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line4)
        c_node_pgm = pgm_obj.get_node_pgm_by_name("c")

        self.assertEqual(3, pgm_obj.n_nodes)
        self.assertEqual(type(c_node_pgm.value), list)
        self.assertEqual(len(c_node_pgm.value), 3)
        self.assertEqual(c_node_pgm.value, ["1", "2", "3"])

        # --- # 

        cmd_line5 = "d <- [c, 4, 5, 6]"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line5)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line5)
        d_node_pgm = pgm_obj.get_node_pgm_by_name("d")

        self.assertEqual(4, pgm_obj.n_nodes)
        self.assertEqual(type(d_node_pgm.value), list)
        self.assertEqual(len(d_node_pgm.value), 6)
        self.assertEqual(d_node_pgm.value, ["1", "2", "3", "4", "5", "6"])

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
    # on the terminal, remember to add "src/" to
    # PYTHONPATH (system variable), or to set it if it does not
    # exist -- don't forget to export it!
    # 
    # Then you can do:
    # $ python3.9 tests/interface/test_cmd_var_assignment
    # 
    # or
    #
    # $ python3.9 -m tests.interface.test_cmd_var_assignment
    #
    # or 
    #
    # $ python3.9 -m unittest tests.interface.test_cmd_var_assignment.TestVarAssignment.test_var_assignment

    unittest.main()