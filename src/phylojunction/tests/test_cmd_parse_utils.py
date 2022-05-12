import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/tests/)
import unittest

# pj imports
import interface.cmd_parse_utils as cmdu

class TestCommandParse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tokenizer_exp_dict = { "name": "\"lambda0\"", "value": "1.0", "event":"\"w_speciation\"", "states":"[0,0,0]" }

    def test_tokenizer(self):
        """
        Test if tokenizer produces the correct dictionary outputs
        """
        fn_spec1 = "name=\"lambda0\", value=1.0, event=\"w_speciation\", states=[0, 0, 0]"
        fn_spec2 = "name=\"lambda0\", states=[0, 0, 0], value=1.0, event=\"w_speciation\""
        fn_spec3 = "states=[0, 0, 0], name=\"lambda0\", value=1.0, event=\"w_speciation\""
        
        cmd_line1 = "a <- sse_rate(" + fn_spec1 + ")"    
        cmd_line2 = "a <- sse_rate(" + fn_spec2 + ")" 
        cmd_line3 = "a <- sse_rate(" + fn_spec3 + ")"     
        
        token_dict1 = cmdu.tokenize_fn_spec(fn_spec1, cmd_line1)
        token_dict2 = cmdu.tokenize_fn_spec(fn_spec2, cmd_line2)
        token_dict3 = cmdu.tokenize_fn_spec(fn_spec3, cmd_line3)

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict1[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict2[k])

        for k, v in self.tokenizer_exp_dict.items():
            self.assertEqual(v, token_dict3[k])

if __name__ == "__main__":
    # can be called from tests/
    # $ python3 test_cmd_parse_utils.py
    # 
    # can also be called from phylojunction/
    # $ python3 tests/test_cmd_parse_utils.py
    # or
    # $ python3 -m tests.test_cmd_parse_utils
    # or, for a specific test
    # $ python3 -m unittest tests.test_cmd_parse_utils.TestCommandParse.test_tokenizer
    #
    # can also be called from VS Code, if open folder is phylojuction/

    unittest.main()