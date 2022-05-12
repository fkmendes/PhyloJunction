import sys
sys.path.append('../')
import unittest
import re

# pj imports
import pgm.pgm as pgm
import interface.cmd_parse as cmd
import interface.cmd_parse_utils as cmdu

class TestVarAssignment(unittest.TestCase):
    def test_var_assignment(self):
        """
        Test if a series of different variable assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """
        pgm = pgm.ProbabilisticGraphicalModel()
        
        cmd_line1 = "a <- 1"

        rv_name, _, rv_spec = re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(pgm, rv_name, rv_spec, cmd_line1)
        a_node_pgm = pgm.get_node_pgm_by_name("a")
        
        self.assertEqual(1, pgm.n_nodes)
        self.assertEqual(type(a_node_pgm.value), list)
        self.assertEqual(len(a_node_pgm.value), 1)
        self.assertEqual(a_node_pgm.value, ["1"])
        
        # --- #

        cmd_line2 = "b <- [1, 2, 3]"

        rv_name, _, rv_spec = re.split(cmdu.assign_regex, cmd_line2)
        cmd.parse_variable_assignment(pgm, rv_name, rv_spec, cmd_line2)
        b_node_pgm = pgm.get_node_pgm_by_name("b")

        self.assertEqual(2, pgm.n_nodes)
        self.assertEqual(type(b_node_pgm.value), list)
        self.assertEqual(len(b_node_pgm.value), 3)
        self.assertEqual(b_node_pgm.value, ["1", "2", "3"])
        
        # --- #
        
        cmd_line3 = "a <- [1]"

        rv_name, _, rv_spec = re.split(cmdu.assign_regex, cmd_line3)
        cmd.parse_variable_assignment(pgm, rv_name, rv_spec, cmd_line3)
        a_node_pgm = pgm.get_node_pgm_by_name("a")

        self.assertEqual(2, pgm.n_nodes)
        self.assertEqual(type(a_node_pgm.value), list)
        self.assertEqual(len(a_node_pgm.value), 1)
        self.assertEqual(a_node_pgm.value, ["1"])

        # --- #
        
        cmd_line4 = "c <- b"

        rv_name, _, rv_spec = re.split(cmdu.assign_regex, cmd_line4)
        cmd.parse_variable_assignment(pgm, rv_name, rv_spec, cmd_line4)
        c_node_pgm = pgm.get_node_pgm_by_name("c")

        self.assertEqual(3, pgm.n_nodes)
        self.assertEqual(type(c_node_pgm.value), list)
        self.assertEqual(len(c_node_pgm.value), 3)
        self.assertEqual(c_node_pgm.value, ["1", "2", "3"])

        # --- # 

        cmd_line5 = "d <- [c, 4, 5, 6]"

        rv_name, _, rv_spec = re.split(cmdu.assign_regex, cmd_line5)
        cmd.parse_variable_assignment(pgm, rv_name, rv_spec, cmd_line5)
        d_node_pgm = pgm.get_node_pgm_by_name("d")

        self.assertEqual(4, pgm.n_nodes)
        self.assertEqual(type(d_node_pgm.value), list)
        self.assertEqual(len(d_node_pgm.value), 6)
        self.assertEqual(d_node_pgm.value, ["1", "2", "3", "4", "5", "6"])

if __name__ == "__main__":
    unittest.main()