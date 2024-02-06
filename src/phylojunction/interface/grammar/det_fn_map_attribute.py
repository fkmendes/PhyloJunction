import typing as ty

# pj imports
import phylojunction.utility.exception_classes as ec
import phylojunction.pgm.pgm as pgm
import phylojunction.data.tree as pjtr
import phylojunction.calculation.discrete_sse as sseobj  # type: ignore

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def make_mapped_ann_tree(det_fn_name: str,
                         det_fn_param_dict:
                         ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) -> \
                            pjtr.AnnotatedTree:
    
    ann_tr_list: ty.List[pjtr.AnnotatedTree] = list()

    # val is a list of strings or nodes
    for arg, val in det_fn_param_dict.items():
        first_element = val[0]
        print(first_element)

        # if arg == "tree":
        #     if not isinstance(first_element, pgm.NodeDAG):
        #         raise(ec.ParseInvalidArgumentError(
        #             "tree",
        #             ""))