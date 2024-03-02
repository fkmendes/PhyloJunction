import os
import typing as ty

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
from phylojunction.constant_function.ct_function_read_tree import CtFnTreeReader

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def make_tree_reader(ct_fn_name: str,
                     ct_fn_param_dict: \
                        ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) \
        -> pgm.ConstantFn:
    #############################
    # IMPORTANT: Default values #
    #############################
    n_samples: int = 1
    n_repl: int = 1
    tr_fp_or_str: str = str()
    node_name_attr: str = str()
    in_file: bool = False

    # args were already grammar-checked in constant_fn_grammar
    for arg, val in ct_fn_param_dict.items():
        # if element in val is string, it remains unchanged,
        # if NodeDAG, we get its string-fied value
        extracted_val: ty.List[str] = \
            pgm.extract_vals_as_str_from_node_dag(val)

        if True:
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
                
                tr_fp_or_str = extracted_val[0].replace('"', '')

                if "string" in ct_fn_param_dict:
                    raise ec.ParseMutuallyExclusiveParametersError(arg, "string")

                if not os.path.isfile(tr_fp_or_str):
                    print('arg', arg, 'tr_fp_or_str', tr_fp_or_str)
                    raise ec.ParsePathDoesNotExistError(
                        arg,
                        tr_fp_or_str)

                in_file = True

            elif arg == "string":
                if len(extracted_val) > 1:
                    raise ec.ParseRequireSingleValueError("read_tree", arg)

                if "file_path" in ct_fn_param_dict:
                    raise ec.ParseMutuallyExclusiveParametersError(arg, "file_path")

                tr_fp_or_str = extracted_val[0].replace('"', '')

            elif arg == "node_name_attr":
                node_name_attr = extracted_val[0].replace('"', '')

    return CtFnTreeReader(n_samples,
                        n_repl,
                        tr_fp_or_str,
                        node_name_attr,
                        in_file)