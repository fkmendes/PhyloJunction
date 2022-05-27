import unittest

# pj imports
import phylojunction.calculation.discrete_sse as sseobj

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class TestDiscreteSSECheckingFunctions(unittest.TestCase):
    
    def test_state2pattern_converter(self):
        """
        Test if conversion between integers representing compound states and bit patterns is correct
        """
        
        n_characters = 3 # regions A, B, C
        n_states_per_char = 2 # presence/absence (bit)
        # compound-states integer coding -> bit pattern
        # Null: 0 -> 000
        # A:    1 -> 100
        # B:    2 -> 010
        # C:    3 -> 001
        # AB:   4 -> 110
        # AC:   5 -> 101
        # BC:   6 -> 011
        # ABC:  7 -> 111
        
        svc = sseobj.StateIntoPatternConverter(n_characters, n_states_per_char)
        
        self.assertEqual(list(svc.int2set_dict.values()), ["000", "100", "010", "001", "110", "101", "011", "111"])
        self.assertEqual(list(svc.set2int_dict.values()), ["0", "1", "2", "3", "4", "5", "6", "7"])

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
    # $ python3 tests/calculation/test_discrete_sse_checking_fns.py
    # 
    # or
    #
    # $ python3 -m tests.calculation.test_discrete_sse_checking_fns
    #
    # or 
    #
    # $ python3 -m unittest tests.calculation.test_discrete_sse_checking_fns.TestDiscreteSSECheckingFunctions.test_state2pattern_converter

    unittest.main()