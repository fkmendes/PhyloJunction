import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/tests/)
import unittest
import re
import math
from statistics import mean, stdev

# pj imports
import pgm.pgm as pgm
import interface.cmd_parse as cmd
import interface.cmd_parse_utils as cmdu
import utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestInferenceRevBayesParametricDn(unittest.TestCase):
    def test_pj2rb_uniform(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        uniform distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        n_repl = 100000
        cmd_line = "u ~ unif(n=1, nr=" + str(n_repl) + ", min=-1.0, max=1.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("u")

        # rev_str_list = a_node_pgm.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    u[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"

        rb_sample_mean = -0.002172273 # mean(u) in RB
        rb_sample_sd = 0.5781801 # stdev(u) in RB
        rb_sample_min = -0.9999855 # min(u) in RB
        rb_sample_max = 0.9999579 # max(u) in RB

        pj_unif_values = a_node_pgm.value # 1000 floats in list

        self.assertAlmostEqual(rb_sample_mean, mean(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_sd, stdev(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_min, min(pj_unif_values), delta=1e-2)
        self.assertAlmostEqual(rb_sample_max, max(pj_unif_values), delta=1e-2)
        

    def test_pj2rb_exponential(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        exponential distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        n_repl = 100000
        cmd_line = "e1 ~ exponential(n=1, nr=" + str(n_repl) + ", rate=0.5)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("e1")

        # rev_str_list = a_node_pgm.sampling_dn.get_rev_inference_spec_info()
        # str_to_run_in_rb_for_unittest = "for (i in 1:" + str(n_repl) + ") {\n" + \
        #     "    e1[i] ~ " + rev_str_list[0] + "\n" + \
        #     "}"

        rb_sample_mean = 2.013223 # mean(e1) in RB
        rb_sample_sd = 2.017317 # stdev(e1) in RB

        pj_exponential_values = a_node_pgm.value # 1000 floats in list

        self.assertAlmostEqual(rb_sample_mean, mean(pj_exponential_values), delta=1e-1)
        self.assertAlmostEqual(rb_sample_sd, stdev(pj_exponential_values), delta=1e-1)


    def test_pj2rb_gamma(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        gamma distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        pass


    def test_pj2rb_normal(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        normal distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        pass


    def test_pj2rb_lognormal(self):
        """
        Test if PJ -> .Rev conversion functionality is correct for
        lognormal distribution, by comparing against expected strings
        and expected sampled values (from running RevBayes)
        """
        pass

if __name__ == "__main__":
    # can be called from tests/
    # $ python3 test_inference_rb_dn_parametric.py
    # 
    # can also be called from phylojunction/
    # $ python3 tests/test_inference_rb_dn_parametric.py
    # or
    # $ python3 -m tests.test_inference_rb_dn_parametric
    # or, for a specific test
    # $ python3 -m unittest tests.test_inference_rb_dn_parametric.TestInferenceRevBayesParametricDn.test_pj2rb_uniform
    #
    # can also be called from VS Code, if open folder is phylojuction/

    unittest.main()