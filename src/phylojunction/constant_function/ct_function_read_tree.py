import os
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class CtFnTreeReader(pgm.ConstantFn):

    CT_FN_NAME = "read_tree"

    # both defaults to 1 in make_tree_reader()
    n_samples: int
    n_repl: int
    
    total_tr_count: int
    tr_fp_or_str: str
    tr_str_list: ty.List[str]

    # validation of pars happens in phylojunction_grammar
    def __init__(self,
                 n_samples: int,
                 n_repl: int,
                 tr_fp_or_str: str,
                 in_file: bool) -> None:

        self.n_samples = n_samples
        self.n_repl = n_repl
        self.tr_fp_or_str = tr_fp_or_str

        # side-effect:
        # (i) initializes self.tr_str_list
        # (ii) initializes self.total_tr_count
        # (ii) updates self.n_samples and self.n_repl
        tr_str_list = \
            self._initialize_tr_str_list(tr_fp_or_str, in_file)


    def _initialize_tr_str_list(self,
                                tr_fp_or_str: str,
                                in_file: bool = False):
        """Place tree strings into list and update class count members

        Args:
            tr_fp_or_str (str): String containing a single tree in Newick format
                or the path to a file containing one or more Newick-formatted
                tree strings
            in_file (bool): If tree string(s) are in a file being passed as argument.
                (True) or if Newick string is being passed directly (False).
            Defaults to 'False'.

        Returns:
            List of tree strings (with all replicates of all samples).
        """

        if in_file:
            with open(tr_fp_or_str, "r") as tr_in:
                self.tr_str_list = tr_in.read().splitlines()

        self.total_tr_count = len(self.tr_str_list)

        if self.n_samples != 1 and self.n_samples != self.total_tr_count:
            # TODO: add proper exception here
            raise RuntimeError

