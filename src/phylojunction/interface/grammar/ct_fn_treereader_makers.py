import os
import typing as ty

# dendropy imports for error checking
from dendropy.dataio.newickreader import NewickReader as dpnr
from dendropy.dataio.tokenizer import Tokenizer as dptk

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
import phylojunction.readwrite.pj_read as pjr

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def make_tree_reader(ct_fn_name: str,
                     ct_fn_param_dict: \
                        ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
        -> pgm.ConstantFn:
    #############################
    # IMPORTANT: Default values #
    #############################
    n_samples: int = 1
    n_repl: int = 1
    tree_fp: str = ""
    tree_str: str = ""

    # args were already grammar-checked in constant_fn_grammar
    for arg, val in ct_fn_param_dict.items():
        # if element in val is string, it remains unchanged,
        # if NodePGM, we get its string-fied value
        extracted_val: ty.List[str] = \
            pgm.extract_value_from_nodepgm(val)

        if arg == "n":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            # only one element always
            try:
                n_samples = int(extracted_val[0])

            except ValueError:
                raise ec.ParseRequireIntegerError("read_tree", arg)
        
        elif arg == "nr":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            # only one element always
            try:
                n_repl = int(extracted_val[0])

            except ValueError:
                raise ec.ParseRequireIntegerError("read_tree", arg)
        
        elif arg == "file_path":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            if "string" in ct_fn_param_dict:
                raise ec.ParseMutuallyExclusiveParametersError(arg, "string")

            if not os.path.isfile(val):
                raise ec.ParsePathDoesNotExistError(arg, extracted_val[0])
            
            tree_fp = extracted_val[0]

        elif arg == "string":
            if len(extracted_val) > 1:
                raise ec.ParseRequireSingleValueError("read_tree", arg)

            if "file_path" in ct_fn_param_dict:
                raise ec.ParseMutuallyExclusiveParametersError(arg, "file_path")

            tree_str = extracted_val[0]

        try:
            # TODO: call pjr.read_nwk_tree_str()
            pass

        except (TypeError,
                AttributeError,
                dptk.UnexpectedEndOfStreamError,
                dpnr.NewickReaderMalformedStatementError,
                dpnr.NewickReaderDuplicateTaxonError) as e:
            if tree_str:
                raise ec.ParseInvalidNewickStringError("string")
            
            elif tree_fp:
                raise ec.ParseInvalidNewickStringError("file_path")