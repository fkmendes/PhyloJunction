import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/tests/)
import unittest

# pj imports
import utility.helper_functions as pjh
import calculation.math_utils as pjmath

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestExpectations(unittest.TestCase):

    def test_yule_height(self):
        """
        Test if expectation of Yule tree height
        """

        exp_h = pjmath.exp_root_height_yule_ntaxa(0.3, 7)
        self.assertAlmostEqual(exp_h, 5.309524, msg="Expected height should be 5.309524.", delta=1e-6)

    def test_bd_extant_count(self):
        """
        Test count of extant species of birth-death tree
        """

        exp_n = pjmath.exp_extant_count_bd(0.9, 0.7, 7.5)
        self.assertAlmostEqual(exp_n, 4.481689, msg="Expected count of extant taxa should be 4.481689.", delta=1e-6)

if __name__ == '__main__':
    unittest.main()