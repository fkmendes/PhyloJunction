import os
import typing as ty

# dendropy imports for error checking
from dendropy.dataio.newickreader import NewickReader as dpnr
from dendropy.dataio.tokenizer import Tokenizer as dptk

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
import phylojunction.readwrite.pj_read as pjr
from phylojunction.data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class CtFnTreeReader(pgm.ConstantFn):

    CT_FN_NAME = "read_tree"

    # both defaults to 1 in make_tree_reader()
    n_samples: int
    n_repl: int
    node_name_attr: str
    
    in_file: bool
    total_tr_count: int
    tr_str_list: ty.List[str]

    # validation of pars happens in phylojunction_grammar
    def __init__(self,
                 n_samples: int,
                 n_repl: int,
                 tr_fp_or_str: str,
                 node_name_attr: str,
                 in_file: bool) -> None:

        self.n_samples = n_samples
        self.n_repl = n_repl
        self.node_name_attr = node_name_attr
        self.in_file = in_file
        self.tr_str_list = list()

        if in_file:
            with open(tr_fp_or_str, "r") as tr_in:
                self.tr_str_list = tr_in.read().splitlines()
        
        else:
             self.tr_str_list = [tr_fp_or_str]

        self.total_tr_count = len(self.tr_str_list)

        # if n_samples and n_repl are not provided, but there is more
        # than one tree, we update self.n_samples, and assume trees are
        # not replicates, but different simulations (samples)
        if self.n_samples == 1 and self.n_repl == 1 and \
                self.total_tr_count > 1:
            self.n_samples = self.total_tr_count

        self.init_check_vectorize_sample_size()

    def init_check_vectorize_sample_size(self) -> None:
        """Check provided tree strings and count members match"""

        arg: str = "string" 
        if self.in_file:
            arg = "file_path"

        if (self.n_samples != 1 or self.n_repl != 1) and \
            (self.n_samples * self.n_repl) != self.total_tr_count:
                raise ec.ObjInitIncorrectDimensionError(
                        "CtFnTreeReader",
                        arg,
                        self.total_tr_count,
                        exp_len=self.n_samples * self.n_repl
                    )

    def get_rev_inference_spec_info(self) -> ty.List[str]:
        pass

    def generate(self) -> ty.List[AnnotatedTree]:
        sampled_values: ty.List[AnnotatedTree] = list()
        
        try:
            for tr_str in self.tr_str_list:
                ann_tr = pjr.read_nwk_tree_str(
                    tr_str,
                    in_file=False,
                    node_names_attribute=self.node_name_attr)
                
                sampled_values.append(ann_tr)

        except (TypeError,
                AttributeError,
                dptk.UnexpectedEndOfStreamError,
                dpnr.NewickReaderMalformedStatementError,
                dpnr.NewickReaderDuplicateTaxonError) as e:
            if not self.in_file:
                raise ec.ParseInvalidNewickStringError("string")
            
            else:
                raise ec.ParseInvalidNewickStringError("file_path")
        
        return sampled_values
            

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
    # $ python3 src/phylojunction/constant_function/ct_function_read_tree.py
    # 
    # or
    #
    # $ python3 -m phylojunction.constant_function.ct_function_read_tree

    n_samples = 1
    n_repl = 1

    tr1 = CtFnTreeReader(1, 1, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    
    try:
        tr_error1 = CtFnTreeReader(1, 2, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    except Exception as e:
        print(e)

    try:
        tr_error2 = CtFnTreeReader(2, 1, "((A:1.0,B:1.0):1.0,C:2.0);", False)
    except Exception as e:
        print(e)
    
    tr2 = CtFnTreeReader(1, 1, "examples/trees_maps_files/tree_to_read.tre", True)
    tr3 = CtFnTreeReader(1, 1, "examples/trees_maps_files/trees_to_read.tre", True)
    
    try:
        tr_error3 = CtFnTreeReader(1, 3, "examples/trees_maps_files/trees_to_read.tre", True)
    except Exception as e:
        print(e)

    try:
        tr_error4 = CtFnTreeReader(3, 1, "examples/trees_maps_files/trees_to_read.tre", True)
    except Exception as e:
        print(e)

    print(tr1.generate())
    print(tr2.generate())
    print(tr3.generate())