import os
import typing as ty

# pj imports
# import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.exception_classes as ec
import phylojunction.pgm.pgm as pgm
import phylojunction.readwrite.pj_read as pjr
from phylojunction.data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class PJCtFnGrammar():

    ct_fn_grammar_dict: ty.Dict[str, ty.Set[str]]

    #########################################
    # All available deterministic functions #
    #########################################
    ct_fn_grammar_dict = {
        "read_tree": set(["file_path", "string", "node_name_attr"])
    }

    @classmethod
    def grammar_check(cls, ct_fn_id: str, fn_param: str) -> bool:
        if fn_param in cls.ct_fn_grammar_dict[ct_fn_id]:
            return True
        return False

    @classmethod
    def init_return_ann_tr(
        cls,
        ct_fn_param_dict:
            ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> AnnotatedTree:

        tree_fp = ""
        tree_str = ""

        if not ct_fn_param_dict:
            raise ec.ParseMissingSpecificationError("read_tree")

        for arg, val in ct_fn_param_dict.items():
            if not cls.grammar_check("read_tree", arg):
                raise ec.ParseNotAParameterError(arg)

            if arg == "file_path":
                if "string" in ct_fn_param_dict:
                    raise ec.ParseMutuallyExclusiveParametersError(arg, "string")

                if not os.path.isfile(val):
                    raise ec.ParsePathDoesNotExistError(arg, val)
                
                tree_fp = val

            if arg == "string":
                if "file_path" in ct_fn_param_dict:
                    raise ec.ParseMutuallyExclusiveParametersError(arg, "file_path")

                tree_str = val

            try:
                # TODO: call pjr.read_nwk_tree_str()
                pass

            except:
                # TODO: capture dendropy's error messages and return custom ec.Exception...
                raise RuntimeError
            
        
    @classmethod
    def create_ct_fn_obj(
        cls,
        ct_fn_id: str,
        ct_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> ty.Optional[
                # ty.Union[
                    AnnotatedTree
                # ]
                ]:
        """
        Build and return constant function object.

        Args:
            ct_fn_id (str): Name of constant function to being called
            ct_fn_param_dict (dict): Dictionary containing constant
                function parameter names (str) as keys and lists (of either
                strings or NodePGMs) as values

        Returns:
            Object: one of a variety of objects to be stored within a clamped
            stochastic node
        """

        # validate input
        if ct_fn_id == "read_tree":
            # health checks inside
            return cls.init_return_ann_tr(ct_fn_param_dict)

        return None