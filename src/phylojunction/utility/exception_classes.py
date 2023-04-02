import typing as ty

class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str

    def __init__(self, cmd_line: str, message: str) -> None:
        self.cmd_line = cmd_line
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: " + self.cmd_line + \
            "\n\nThe line above had a syntax problem. " + self.message


class InexistentVariableError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str) -> None:
        self.rv_name = rv_name
        self.message = ""
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: There was a problem when assigning \'" \
            + self.rv_name + "\' as an argument. It seems this " \
            + "variable has not been not initialized."


class VariableAssignmentError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str, message: str="") -> None:
        self.message = message
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Cannot assign values of different types to a " \
            + "variable (\'" + self.rv_name + "\'). Exiting..."


class VariableMisspec(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str) -> None:
        self.message = ""
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Could not parse value out of variable \'" \
            + self.rv_name + "\'."


class NoSpecificationError(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Function call had not arguments. "+ self.message


class FunctionArgError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str) -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Could not set value for function parameter " \
            + self.par_name + ". " + self.message


# TODO: coalesce with exception above, maybe grab arg name?
class InvalidFunctionArgError(Exception):
    function_name: str
    message: str

    def __init__(self, function_name: str, message: str) -> None:
        self.function_name = function_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Could not specify function " + self.function_name + ". " + self.message


class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str

    def __init__(self, fn_name: str, message: str) -> None:
        self.fn_name = fn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Could not call function " + self.fn_name \
            + " because argument values are incongruent in some way. " \
            + self.message


class NotAParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str="") -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: \'" + self.par_name + "\' does not seem to be " \
            + "a valid parameter." + self.message


class NotBetweenZeroAndOneError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str="") -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: \'" + self.par_name + "\' must be contained in " \
            + "[0,1], but one or more of its values was " + self.message \
            + "."


class WrongDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str

    def __init__(self, container_name: str, obs_len: int,
        exp_len: int=0, message: str="") -> None:
        self.obj_name = container_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: Container " + self.obj_name + " had a different " \
            + "dimension than expected. " + self.message + "."


class DimensionalityError(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str) -> None:
        self.dn_name = dn_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When specifying distribution " + self.dn_name + ":\n" + \
            "  (i) the number of values provided for each parameters differed, or\n" + \
            "  (ii) the number(s) was larger than the specified number of samples (i.e., draws), or\n" + \
            "  (iii) the number(s) was larger than 1, but smaller than the specified number of samples.\nExiting..."


class NodeInferenceDimensionalityError(Exception):
    node_name: str
    message: str

    def __init__(self, node_name: str) -> None:
        self.node_name = node_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When converting values of node " + self.node_name + " to inference script, the " + \
            "number of values as larger than 1, but smaller than the number of samples specified for " + \
            "other nodes. Unsure how to parse values. Exiting..."


class NoPlatingAllowedError(Exception):
    det_name: str
    message: str
    node_pgm_name: str

    def __init__(self, det_name: str, problematic_node_pgm_name: str, message: str="") -> None:
        self.det_name = det_name
        self.message = message
        self.node_pgm_name = problematic_node_pgm_name
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When executing " + self.det_name + "(), replicates were detected for argument " + \
            self.node_pgm_name + ". Plating is not supported for this deterministic function.\nExiting... "


class RequireScalarError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = ""
        
    def __str__(self) -> str:
        return "ERROR: When specifying object " + self.obj_name \
            + "'s parameter \'" + self.arg + "\', more than one value was " \
            + "provided. A scalar is required."


class RequireSameParameterType(Exception):
    message: str
    obj_name: str
    n_diff_par: int

    def __init__(self, obj_name: str, n_diff_par: int) -> None:
        self.obj_name = obj_name
        self.n_diff_par = n_diff_par
        self.message = ""
        
    def __str__(self) -> str:
        return "ERROR: When specifying object " + self.obj_name \
            + " only one type of parameter is allowed. Found " \
            + str(self.n_diff_par) + "."


class ReplicateNumberError(Exception):
    node_name: str
    message: str

    def __init__(self, node_name, message="") -> None:
        self.node_name = node_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When going through values in " + self.node_name + ", the number of replicates differed among simulations. One or more of the simulations failed to generate the specified number of replicates. Exiting..."


class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str

    def __init__(self, rv_name: str, dn_name: str, message: str="") -> None:
        self.rv_name = rv_name
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "\Warning: Distribution " + self.dn_name + " will be called to sample " + self.rv_name + " using vectorization."


class MissingStateDependentParameterError(Exception):
    epoch_missing_param: int
    symmetric_diff_set: ty.Set[ty.Any]
    symmetric_diff_str: str
    message: str

    def __init__(self, epoch_missing_param: int,
        symmetric_diff_set: ty.Set[ty.Any],
        message: str="") -> None:

        self.epoch_missing_param = epoch_missing_param
        self.symmetric_diff_str = "( "
        
        for i in symmetric_diff_set:
            ith_str = ""
            
            if isinstance(i, int):
                ith_str = str(i)
            
            elif isinstance(i, tuple):
                ith_str = "(" + ", ".join(str(j) for j in i) + ")"
            
            if self.symmetric_diff_str == "( ":
                self.symmetric_diff_str += ith_str
            
            else:
                self.symmetric_diff_str += ", " + ith_str
        
        self.symmetric_diff_str += " )"
        
        self.message = message
        
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: One or more state-dependent parameters were " \
            + "missing at least from time slice " \
            + str(self.epoch_missing_param) + ". The symmetric " \
            + "difference between time-slice state sets is " \
            + self.symmetric_diff_str


class StateDependentParameterMisspec(Exception):
    message: str

    def __init__(self, message: str="") -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Misspecified state-dependent parameter input. " \
            + self.message


class SSEStashMisspec(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Misspecified SSE stash parameter input. " + self.message


class DnInitMisspec(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Distribution " + self.dn_name + " was not properly initialized. " + self.message


class InvalidMCMCChainLength(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: " + self.message


# Tree exceptions #
class AnnotatedTreeMisspec(Exception):
    message: str

    def __init__(self, message) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "ERROR: Misspecified AnnotatedTree input. " + self.message


class AnnotatedTreeLineageMissannotation(Exception):
    message: str
    
    def __init__(self, message) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: AnnotatedTree cannot be annotated this way. " + self.message