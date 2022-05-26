import unittest
import re
import math
from statistics import mean

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.interface.cmd.cmd_parse as cmd
import phylojunction.interface.cmd.cmd_parse_utils as cmdu
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestParametricSamplingDnAssignment(unittest.TestCase):
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
        

    def test_unif_misspec(self):
        """
        Test uniform sampling distribution assignment throws exception
        if mandatory parameters are missing
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "u ~ unif(n=1, nr=1, max=1.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Uniform\" was not properly initialized. Parameter \"min\" is missing.")
        
        cmd_line2 = "u ~ unif(n=1, nr=1, min=0.0)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Uniform\" was not properly initialized. Parameter \"max\" is missing.")


    def test_sampling_exp_assignment(self):
        """
        Test if exponential sampling distribution assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """

        ################################################
        # Exponential, rate parameterization (default) #
        ################################################
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "e1 ~ exponential(n=100000, nr=1, rate=0.5)"
        
        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        a_node_pgm1 = pgm_obj.get_node_pgm_by_name("e1")

        self.assertTrue(isinstance(a_node_pgm1.value, list))
        self.assertEqual(len(a_node_pgm1.value), 100000)
        self.assertAlmostEqual(2.0, mean(a_node_pgm1.value), delta=0.05)
        self.assertEqual(1, pgm_obj.n_nodes)

        #######################################
        # Exponential, scale parameterization #
        #######################################
        cmd_line2 = "e2 ~ exponential(n=100000, nr=1, rate=0.5, rate_parameterization=\"false\")"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_pgm2 = pgm_obj.get_node_pgm_by_name("e2")

        self.assertTrue(isinstance(a_node_pgm2.value, list))
        self.assertEqual(len(a_node_pgm2.value), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_pgm2.value), delta=0.05)
        self.assertEqual(2, pgm_obj.n_nodes)


    def test_exp_misspec(self):
        """
        Test exponential sampling distribution assignment throws exception
        if mandatory parameters are missing
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line = "e ~ exponential(n=1, nr=1)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Exponential\" was not properly initialized. Parameter \"rate\" is missing.")


    def test_sampling_gamma_assignment(self):
        """
        Test if gamma sampling distribution assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """

        ###########################################
        # Gamma, scale parameterization (default) #
        ###########################################
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "g1 ~ gamma(n=100000, nr=1, shape=0.5, scale=0.5)"
        
        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        a_node_pgm1 = pgm_obj.get_node_pgm_by_name("g1")
        
        self.assertTrue(isinstance(a_node_pgm1.value, list))
        self.assertEqual(len(a_node_pgm1.value), 100000)
        self.assertAlmostEqual(0.25, mean(a_node_pgm1.value), delta=0.05)
        self.assertEqual(1, pgm_obj.n_nodes)

        ################################
        # Gamma, rate parameterization #
        ################################
        cmd_line2 = "g2 ~ gamma(n=100000, nr=1, shape=0.5, scale=0.5, rate_parameterization=\"true\")"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)
        a_node_pgm2 = pgm_obj.get_node_pgm_by_name("g2")

        self.assertTrue(isinstance(a_node_pgm2.value, list))
        self.assertEqual(len(a_node_pgm2.value), 100000)
        self.assertAlmostEqual(1.0, mean(a_node_pgm2.value), delta=0.05)
        self.assertEqual(2, pgm_obj.n_nodes)


    def test_gamma_misspec(self):
        """
        Test gamma sampling distribution assignment throws exception
        if mandatory parameters are missing
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "g ~ gamma(n=1, nr=1, shape=0.5)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Gamma\" was not properly initialized. Parameter \"scale\" is missing.")
        
        cmd_line2 = "g ~ gamma(n=1, nr=1, scale=0.5)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Gamma\" was not properly initialized. Parameter \"shape\" is missing.")


    def test_sampling_normal_assignment(self):
        """
        Test if normal sampling distribution assignments are correctly evaluated
        and result in the right probabilistic graphical model
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()
        
        cmd_line1 = "n ~ normal(n=100000, mean=0.5, sd=0.1)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)
        cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)
        a_node_pgm = pgm_obj.get_node_pgm_by_name("n")
        
        self.assertTrue(isinstance(a_node_pgm.value, list))
        self.assertEqual(len(a_node_pgm.value), 100000)
        self.assertAlmostEqual(0.5, mean(a_node_pgm.value), delta=0.1)
        self.assertEqual(1, pgm_obj.n_nodes)


    def test_normal_misspec(self):
        """
        Test normal sampling distribution assignment throws exception
        if mandatory parameters are missing
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "n ~ normal(n=1, nr=1, sd=0.1)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Normal\" was not properly initialized. Parameter \"mean\" is missing.")
        
        cmd_line2 = "n ~ normal(n=1, nr=1, mean=0.5)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Normal\" was not properly initialized. Parameter \"sd\" is missing.")


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


    def test_ln_misspec(self):
        """
        Test log-normal sampling distribution assignment throws exception
        if mandatory parameters are missing
        """
        
        pgm_obj = pgm.ProbabilisticGraphicalModel()

        cmd_line1 = "ln ~ lognormal(n=1, nr=1, sd=0.25)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line1)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line1)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Log-normal\" was not properly initialized. Parameter \"mean\" is missing.")
        
        cmd_line2 = "ln ~ lognormal(n=1, nr=1, mean=-3.25)"

        stoch_node_name, _, stoch_node_spec = re.split(cmdu.sampled_as_regex, cmd_line2)

        with self.assertRaises(ec.DnInitMisspec) as exc:
            cmd.parse_samp_dn_assignment(pgm_obj, stoch_node_name, stoch_node_spec, cmd_line2)        
        self.assertEqual(str(exc.exception), "\nERROR: Distribution \"Log-normal\" was not properly initialized. Parameter \"sd\" is missing.")

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
    # $ python3 tests/distribution/test_cmd_parametric_sampling_dn_assignment.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_cmd_parametric_sampling_dn_assignment
    #
    # or 
    #
    # $ python3 -m unittest tests.test_cmd_parametric_sampling_dn_assignment.TestParametricSamplingDnAssignment.test_sampling_unif_assignment
    
    unittest.main()