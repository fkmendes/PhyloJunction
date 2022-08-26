class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str

    def __init__(self, cmd_line: str, message: str) -> None:
        self.cmd_line = cmd_line
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: " + self.cmd_line + \
            "\nThe line above had a syntax problem. " + self.message

class InexistentVariableError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str) -> None:
        self.rv_name = rv_name
        self.message = ""
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: There was a problem when assigning \'" + self.rv_name + \
            "\' as an argument. It seems this variable has not been not initialized."

class VariableAssignmentError(Exception):
    rv_name: str
    message: str

    def __init__(self, rv_name: str, message: str="") -> None:
        self.message = message
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Cannot assign values of different types to a variable (\'" + self.rv_name + "\'). Exiting..."

class VariableMisspec(Exception):
    rv_name: str

    def __init__(self, rv_name: str) -> None:
        self.message = ""
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Could not parse value out of variable \'" + self.rv_name + "\'."

class NoSpecificationError(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Function call had not arguments. "+ self.message

class FunctionArgError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str) -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Could not set value for function parameter " + self.par_name + ". " + self.message

class InvalidFunctionArgError(Exception):
    function_name: str
    message: str

    def __init__(self, function_name: str, message: str) -> None:
        self.function_name = function_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Could not specify function " + self.function_name + ". " + self.message

class FunctionArgsMismatchError(Exception):
    fn_name: str
    message: str

    def __init__(self, fn_name: str, message: str) -> None:
        self.fn_name = fn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Could not call function " + self.fn_name + " because argument values are incongruent in some way. " + self.message + ". Exiting... "

class NotAParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str="") -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: \'" + self.par_name + "\' does not seem to be a valid parameter." + self.message + " Exiting... "

class WrongDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str

    def __init__(self, container_name: str, obs_len: int, exp_len: int=0, message: str="") -> None:
        self.obj_name = container_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "\nERROR: Container " + self.obj_name + " had a different dimension than expected. " + self.message + ". Exiting... "

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
            "  (ii) the number(s) was larger than the specified number of draws, or\n" + \
            "  (iii) thet number(s) was larger than 1, but smaller than the specified number of draws.\nExiting..."

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
    dn_name: str
    message: str

    def __init__(self, dn_name: str, arg: str) -> None:
        self.dn_name = dn_name
        self.arg = arg
        self.message = ""
        
    def __str__(self) -> str:
        return "ERROR: When specifying distribution " + self.dn_name + "'s parameter " + self.arg + ", more than one value was provided. A scalar is required. Exiting..."

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

class SSEAtomicRateMisspec(Exception):
    message: str

    def __init__(self, message: str="") -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Misspecified SSE atomic rate parameter input. " + self.message

class SSEStashMisspec(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Misspecified SSE stash parameter input. " + self.message

class AnnotatedTreeMisspec(Exception):
    message: str

    def __init__(self, message) -> None:
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Misspecified AnnotatedTree input. " + self.message

class DnInitMisspec(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self) -> str:
        return "\nERROR: Distribution " + self.dn_name + " was not properly initialized. " + self.message

class InvalidMCMCChainLength(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "\nERROR: " + self.message