import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/tests/)
import unittest
import re
import math
from statistics import mean

# pj imports
import pgm.pgm as pgm
import interface.cmd_parse as cmd
import interface.cmd_parse_utils as cmdu

class TestSamplingDnAssignment(unittest.TestCase):
    def test_sampling_unif_assignment(self):
        """
        Test if uniform sampling distribution assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        cmd_line = "u ~ unif(n=100000, min=-1.0, max=1.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("u")
  
        self.assertTrue(isinstance(a_node_pgm.value, list))
        self.assertEqual(len(a_node_pgm.value), 100000)
        self.assertAlmostEqual(0.0, mean(a_node_pgm.value), delta=1e-2)
        self.assertLessEqual(-1.0, min(a_node_pgm.value))
        self.assertGreater(1.0, max(a_node_pgm.value))
        self.assertEqual(1, pgm_obj.n_nodes)


    def test_sampling_unif_vectorized_assignment(self):
        """
        Test if uniform (with vectorized inputs) sampling distribution assignments
        are correctly evaluated and result in the right probabilistic graphical model
        """
        
        mins = [0.0, 1.0, 2.0, 3.0, 4.0]
        maxs = [1.0, 2.0, 3.0, 4.0, 5.0]
        tups = tuple((i, maxs[idx])for idx, i in enumerate(mins))

        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        cmd_line1 = "mins <- [0.0, 1.0, 2.0, 3.0, 4.0]"
        
        stoch_node_name, _, stoch_node_spec = re.split(cmdu.assign_regex, cmd_line1)
        cmd.parse_variable_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        
        cmd_line2 = "u ~ unif(n=5, nr=2, min=mins, max=[1.0, 2.0, 3.0, 4.0, 5.0])"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("u")
  
        self.assertTrue(isinstance(a_node_pgm.value, list))
        for idx, tup in enumerate(tups):
            self.assertTrue(tup[0] <= a_node_pgm.value[idx * 2] < tup[1])
            self.assertTrue(tup[0] <= a_node_pgm.value[idx * 2 + 1] < tup[1])
        self.assertEqual(2, pgm_obj.n_nodes)
        

    def test_sampling_ln_assignment(self):
        """
        Test if log-normal sampling distribution assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """
        
        ###################################
        # Log-normal, log space (default) #
        ###################################
        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        cmd_line1 = "ln1 ~ lognormal(n=100000, mean=-3.25, sd=0.25)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("ln1")
        
        self.assertTrue(isinstance(a_node_pgm.value, list))
        self.assertEqual(len(a_node_pgm.value), 100000)
        self.assertAlmostEqual(2.37, mean(a_node_pgm.value), delta=0.1)
        self.assertEqual(1, pgm_obj.n_nodes)
        
        ##########################
        # Log-normal, real space #
        ##########################

        exp_mean = str(math.exp(-3.25))
        exp_sd = str(math.exp(0.25))
        cmd_line2 = "ln2 ~ lognormal(n=100000, mean=" + exp_mean + ", sd=" + exp_sd + ", log_space=\"false\")"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("ln2")

        self.assertTrue(isinstance(a_node_pgm.value, list))
        self.assertEqual(len(a_node_pgm.value), 100000)
        self.assertAlmostEqual(2.37, mean(a_node_pgm.value), delta=0.1)
        self.assertEqual(2, pgm_obj.n_nodes)

if __name__ == "__main__":
    # can be called from tests/
    # $ python3 test_cmd_sampling_dn_assignment.py
    # 
    # can also be called from phylojunction/
    # $ python3 tests/test_cmd_sampling_dn_assignment.py
    # or
    # $ python3 -m tests.test_cmd_sampling_dn_assignment
    # or, for a specific test
    # $ python3 -m unittest tests.test_cmd_sampling_dn_assignment.TestSamplingDnAssignment.test_sampling_unif_assignment
    #
    # can also be called from VS Code, if open folder is phylojuction/

    unittest.main()