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

class NoSpecificationError(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

class FunctionArgError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str) -> None: ...

class InvalidFunctionArgError(Exception):
    function_name: str
    message: str
    def __init__(self, function_name: str, message: str) -> None: ...

class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str
    def __init__(self, fn_name: str, message: str) -> None: ...

class NotAParameterError(Exception):
    par_name: str
    message: str
    def __init__(self, par_name: str, message: str = ...) -> None: ...

class WrongDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str
    obj_name: Incomplete
    def __init__(self, container_name: str, obs_len: int, exp_len: int = ..., message: str = ...) -> None: ...

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

class RequireScalarError(Exception):
    dn_name: str
    message: str
    arg: Incomplete
    def __init__(self, dn_name: str, arg: str) -> None: ...

class ReplicateNumberError(Exception):
    node_name: str
    message: str
    def __init__(self, node_name, message: str = ...) -> None: ...

class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str
    def __init__(self, rv_name: str, dn_name: str, message: str = ...) -> None: ...

class SSEAtomicRateMisspec(Exception):
    message: str
    def __init__(self, message: str = ...) -> None: ...

class SSEStashMisspec(Exception):
    message: str
    def __init__(self, message: str) -> None: ...

class AnnotatedTreeMisspec(Exception):
    message: str
    def __init__(self, message) -> None: ...

class DnInitMisspec(Exception):
    dn_name: str
    message: str
    def __init__(self, dn_name: str, message: str) -> None: ...

class InvalidMCMCChainLength(Exception):
    message: str
    def __init__(self, message: str) -> None: ...
