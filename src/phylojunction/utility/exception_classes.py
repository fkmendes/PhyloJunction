class ScriptSyntaxError(Exception):
    def __init__(self, cmd_line, message):
        self.cmd_line = cmd_line
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: " + self.cmd_line + \
            "\nThe line above had a syntax problem. " + self.message

class InexistentVariableError(Exception):
    def __init__(self, rv_name):
        self.rv_name = rv_name
        self.message = ""
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: There was a problem when assigning \'" + self.rv_name + \
            "\' as an argument. It seems this variable has not been not initialized."

class VariableAssignmentError(Exception):
    def __init__(self, rv_name):
        self.message = ""
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Cannot assign values of different types to a variable (\'" + self.rv_name + "\'). Exiting..."

class VariableMisspec(Exception):
    def __init__(self, rv_name):
        self.message = ""
        self.rv_name = rv_name
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Could not parse value out of variable \'" + self.rv_name + "\'."

class NoSpecificationError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Function call had not arguments. "+ self.message

class FunctionArgError(Exception):
    def __init__(self, par_name, message):
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Could not set value for function parameter " + self.par_name + ". " + self.message

class InvalidFunctionArgError(Exception):
    def __init__(self, function_name, message):
        self.function_name = function_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Could not specify function " + self.function_name + ". " + self.message

class FunctionArgsMismatchError(Exception):
    def __init__(self, fn_name, message):
        self.fn_name = fn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Could not call function " + self.fn_name + " because argument values are incongruent in some way. " + self.message + ". Exiting... "

class NotAParameterError(Exception):
    def __init__(self, par_name, message=""):
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: \'" + self.par_name + "\' does not seem to be a valid parameter." + self.message + " Exiting... "

class WrongDimensionError(Exception):
    def __init__(self, container_name, obs_len, exp_len=0, message=""):
        self.obj_name = container_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "\nERROR: Could not set value for function parameter " + self.par_name + ". " + self.message + ". Exiting... "

class DimensionalityError(Exception):
    def __init__(self, dn_name):
        self.dn_name = dn_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self):
        return "ERROR: When specifying distribution " + self.dn_name + ":\n" + \
            "  (i) the number of values provided for each parameters differed, or\n" + \
            "  (ii) the number(s) was larger than the specified number of draws, or\n" + \
            "  (iii) thet number(s) was larger than 1, but smaller than the specified number of draws.\nExiting..."

class NodeInferenceDimensionalityError(Exception):
    def __init__(self, node_name):
        self.node_name = node_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self):
        return "ERROR: When converting values of node " + self.node_name + " to inference script, the " + \
            "number of values as larger than 1, but smaller than the number of samples specified for " + \
            "other nodes. Unsure how to parse values. Exiting..."

class NoPlatingAllowedError(Exception):
    def __init__(self, det_name, problematic_node_pgm_name, message=""):
        self.det_name = det_name
        self.message = message
        self.node_pgm_name = problematic_node_pgm_name
        super().__init__(self.message)

    def __str__(self):
        return "ERROR: When executing " + self.det_name + "(), replicates were detected for argument " + \
            self.node_pgm_name + ". Plating is not supported for this deterministic function.\nExiting... "

class RequireScalarError(Exception):
    def __init__(self, dn_name, arg, message=""):
        self.dn_name = dn_name
        self.arg = arg
        self.message = message
        
    def __str__(self):
        return "ERROR: When specifying distribution " + self.dn_name + "'s parameter " + self.arg + ", more than one value was provided. A scalar is required. Exiting..."

class ReplicateNumberError(Exception):
    def __init__(self, node_pgm_name, message=""):
        self.node_pgm_name = node_pgm_name
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "ERROR: When going through values in " + self.node_pgm_name + ", the number of replicates differed among simulations. One or more of the simulations failed to generate the specified number of replicates. Exiting..."

class DimensionalityWarning(Exception):
    def __init__(self, rv_name, dn_name, message=""):
        self.rv_name = rv_name
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "\Warning: Distribution " + self.dn_name + " will be called to sample " + self.rv_name + " using vectorization."

class SSEAtomicRateMisspec(Exception):
    def __init__(self, message=""):
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Misspecified SSE atomic rate parameter input. " + self.message

class SSEStashMisspec(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Misspecified SSE stash parameter input. " + self.message

class AnnotatedTreeMisspec(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Misspecified AnnotatedTree input. " + self.message

class DnInitMisspec(Exception):
    def __init__(self, dn_name, message):
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)
    
    def __str__(self):
        return "\nERROR: Distribution " + self.dn_name + " was not properly initialized. " + self.message

class InvalidMCMCChainLength(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return "\nERROR: " + self.message