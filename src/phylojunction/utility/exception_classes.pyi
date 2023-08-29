import typing as ty
from _typeshed import Incomplete

class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str
    def __init__(self, cmd_line: str, message: str) -> None: ...

class InexistentVariableError(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str) -> None: ...

class VariableAssignmentError(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str, message: str = ...) -> None: ...

class VariableMisspec(Exception):
    rv_name: str
    message: str
    def __init__(self, rv_name: str) -> None: ...

class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str
    def __init__(self, fn_name: str, message: str) -> None: ...

class NotBetweenZeroAndOneError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str="") -> None: ...
    def __str__(self) -> str: ...

class IncorrectDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str
    obj_name: Incomplete
    def __init__(self, container_name: str, obs_len: int, exp_len: int = ...) -> None: ...

class DimensionalityError(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str) -> None: ...

class NodeInferenceDimensionalityError(Exception):
    node_name: str
    message: str
    def __init__(self, node_name: str) -> None: ...

class NoPlatingAllowedError(Exception):
    det_name: str
    message: str
    node_pgm_name: str
    def __init__(self, det_name: str, problematic_node_pgm_name: str, message: str = ...) -> None: ...

class RequireSameParameterType(Exception):
    message: str
    obj_name: str
    n_diff_par: int
    def __init__(self, obj_name: str, n_diff_par: int) -> None: ...

class ReplicateNumberError(Exception):
    node_name: str
    message: str
    def __init__(self, node_name, message: str = ...) -> None: ...

class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str
    def __init__(self, rv_name: str, dn_name: str, message: str = ...) -> None: ...

class MissingStateDependentParameterError(Exception):
    epoch_missing_param: int
    symmetric_diff_set: ty.Set[ty.Any]
    message: str
    def __init__(self, epoch_missing_param: int, symmetric_diff_set: ty.Set[ty.Any], message: str="") -> None: ...

class RepeatedStateDependentParameterError(Exception):
    epoch_w_repeated_param: int
    repeated_state: ty.Union[int, ty.Tuple[int]]
    message: str
    def __init__(self, epoch_w_repeated_param: int, repeated_state: ty.Union[int, ty.Tuple[int]], message: str="") -> None: ...

class StateDependentParameterMisspec(Exception):
    message: str
    def __init__(self, message: str = ...) -> None: ...

class SSEStashMisspec(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

class InvalidMCMCChainLength(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

# Initialization errors
class ObjInitInvalidArgError(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str, par_name: str, message: str = "") -> None: ...
    def __str__(self) -> str: ...

class ObjInitMissingParameterError(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str, par_name: str, message: str = "") -> None: ...
    def __str__(self) -> str: ...

class ObjInitIncorrectDimensionError(Exception):
    obj_name: str
    message: str
    def __init__(self, obj_name: str, container_name: str, obs_len: int, exp_len: int = 0) -> None: ...
    def __str__(self) -> str: ...

# Syntax errors: Sampling distribution, det. functions exceptions #
class ParseFunctionArgError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str) -> None: ...

class ParseNotAParameterError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseRequireSingleValueError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ParseRequireIntegerError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ParseRequireNumericError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ParseMissingParameterError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseMissingArgumentError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str) -> None: ...

class ParseMissingSpecificationError(Exception):
    message: str
    def __init__(self, obj2spec_name: str) -> None: ...

class ParseDnInitFailError(Exception):
    dn_name: str
    def __init__(self, dn_name: str, message: str) -> None: ...
    
class ParseDetFnInitFailError(Exception):
    det_fn_name: str
    def __init__(self, det_fn_name: str, message: str) -> None: ...

# PGM exceptions #
class NodePGMNodeStatCantFloatError(Exception):
    message: str
    def __init__(self, node_name: str) -> None: ...
    def __str__(self) -> str: ...
    
# Tree exceptions #
class AnnotatedTreeMisspec(Exception):
    message: str
    def __init__(self, message) -> None: ...
    def __str__(self) -> str: ...

class AnnotatedTreeLineageMissannotation(Exception):
    message: str
    def __init__(self, message) -> None: ...
    def __str__(self) -> str: ...

# Generation exceptions #
class GenerateFailError(Exception):
    message: str
    def __init__(self, dn_name: str, message: str) -> None: ...
    def __str__(self) -> str: ...