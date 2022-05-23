import sys
sys.path.extend(["../", "../phylojunction"])
import unittest
import math
from statistics import mean

# pj imports
import distribution.dn_parametric as dnpar

class TestParametricDns(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dummy_parent_node_tracker = dict()

    def test_unif(self):
        """Test uniform draws with non-vectorized input"""
        # unif_rv = dnpar.DnUnif([100000, 1, [-1.0], [1.0]]).generate()
        unif_rv = dnpar.DnUnif(100000, 1, [-1.0], [1.0], self.dummy_parent_node_tracker).generate()
        self.assertAlmostEqual(0.0, mean(unif_rv), delta=1e-2)
        self.assertLessEqual(-1.0, min(unif_rv))
        self.assertGreater(1.0, max(unif_rv))

    def test_unif_vectorized(self):
        """Test uniform draws with vectorized input"""
        repl_size = 2
        mins = [0.0, 1.0, 2.0, 3.0, 4.0]
        maxs = [1.0, 2.0, 3.0, 4.0, 5.0]
        tups = tuple((i, maxs[idx])for idx, i in enumerate(mins))
        # unif_rv = dnpar.DnUnif([5, repl_size, mins, maxs]).generate()
        unif_rv = dnpar.DnUnif(5, repl_size, mins, maxs, self.dummy_parent_node_tracker).generate()

        for idx, tup in enumerate(tups):
            self.assertTrue(tup[0] <= unif_rv[idx * repl_size] < tup[1])
            self.assertTrue(tup[0] <= unif_rv[idx * repl_size + 1] < tup[1])

    def test_exp(self):
        """Test exponential draws with non-vectorized input"""
        # exp_rv1 = dnpar.DnExponential([100000, 1, 0.5, True]).generate() # rate parameterization (True is default)
        # exp_rv2 = dnpar.DnExponential([100000, 1, 0.5, False]).generate()
        exp_rv1 = dnpar.DnExponential(100000, 1, 0.5, True, self.dummy_parent_node_tracker).generate() # rate parameterization (True is default)
        exp_rv2 = dnpar.DnExponential(100000, 1, 0.5, False, self.dummy_parent_node_tracker).generate()
        self.assertAlmostEqual(2.0, mean(exp_rv1), delta=0.05)
        self.assertAlmostEqual(0.5, mean(exp_rv2), delta=0.05)

    def test_gamma(self):
        """Test gamma draws with non-vectorized input"""
        # gamma_rv1 = dnpar.DnGamma([100000, 1, [0.5], [0.5], False]).generate() # rate parameterization (False is default)
        # gamma_rv2 = dnpar.DnGamma([100000, 1, [0.5], [0.5], True]).generate() # rate parameterization (False is default)
        gamma_rv1 = dnpar.DnGamma(100000, 1, [0.5], [0.5], False, self.dummy_parent_node_tracker).generate() # rate parameterization (False is default)
        gamma_rv2 = dnpar.DnGamma(100000, 1, [0.5], [0.5], True, self.dummy_parent_node_tracker).generate() # rate parameterization (False is default)
        self.assertAlmostEqual(0.25, mean(gamma_rv1), delta=0.05)
        self.assertAlmostEqual(1.0, mean(gamma_rv2), delta=0.05)

    def test_normal(self):
        """Test normal draws with non-vectorized input"""
        # normal_rv = dnpar.DnNormal([100000, 1, [0.5], [0.1]]).generate()
        normal_rv = dnpar.DnNormal(100000, 1, [0.5], [0.1], self.dummy_parent_node_tracker).generate()
        self.assertAlmostEqual(0.5, mean(normal_rv), delta=0.05)

    def test_lognormal(self):
        """Test log-normal draws with non-vectorized input"""
        # lognormal_rv1 = dnpar.DnLogNormal([100000, 1, [-3.25], [0.25], True]).generate() # log-space (True is default)
        # lognormal_rv2 = dnpar.DnLogNormal([100000, 1, [math.exp(-3.25)], [math.exp(0.25)], False]).generate() # real space (log-space is False), R's default when "mean" and "sd" are used
        lognormal_rv1 = dnpar.DnLogNormal(100000, 1, [-3.25], [0.25], True, self.dummy_parent_node_tracker).generate() # log-space (True is default)
        lognormal_rv2 = dnpar.DnLogNormal(100000, 1, [math.exp(-3.25)], [math.exp(0.25)], False, self.dummy_parent_node_tracker).generate() # real space (log-space is False), R's default when "mean" and "sd" are used
        self.assertAlmostEqual(2.37, mean(lognormal_rv1), delta=0.1)
        self.assertAlmostEqual(2.37, mean(lognormal_rv2), delta=0.1)

if __name__ == '__main__':
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
    # $ python3 tests/distribution/test_dn_parametric_object.py
    # 
    # or
    #
    # $ python3 -m tests.distribution.test_dn_parametric_object
    #
    # or 
    #
    # $ python3 -m unittest tests.distribution.test_dn_parametric_object.TestParametricDns.test_unif
    
    unittest.main()