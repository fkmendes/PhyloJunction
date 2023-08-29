import typing as ty


class ScriptSyntaxError(Exception):
    cmd_line: str
    message: str

    def __init__(self, cmd_line: str, message: str) -> None:
        self.cmd_line = cmd_line
        self.message = "\n\nERROR: " \
            + self.cmd_line \
            + ("\n\nThe line above had a syntax problem and could not "
               "be tokenized. ") \
            + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


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

    def __init__(self, rv_name: str, message: str = "") -> None:
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


class NotBetweenZeroAndOneError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str = "") -> None:
        self.par_name = par_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: \'" + self.par_name + "\' must be contained in " \
            + "[0,1], but one or more of its values was " + self.message \
            + "."


class IncorrectDimensionError(Exception):
    container_name: str
    obs_len: int
    exp_len: int
    message: str

    def __init__(self,
                 container_name: str,
                 obs_len: int,
                 exp_len: int = 0) -> None:
        self.obj_name = container_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = \
            "Incorrect dimension of container " + self.obj_name \
            + ", which was of size " + str(obs_len)  + ". "
        
        if exp_len != 0:
            self.message += "The expected dimension was " + str(exp_len) \
            + "."
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class DimensionalityError(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str) -> None:
        self.dn_name = dn_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When specifying distribution " + self.dn_name \
            + ":\n(i) the number of values provided for each parameters " \
            + "differed, or\n (ii) the number(s) was larger than the " \
            + "specified number of samples (i.e., draws), or\n (iii) the " \
            + "number(s) was larger than 1, but smaller than the specified " \
            + "number of samples."


class NodeInferenceDimensionalityError(Exception):
    node_name: str
    message: str

    def __init__(self, node_name: str) -> None:
        self.node_name = node_name
        self.message = ""
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When converting values of node " + self.node_name \
            + " to inference script, the number of values as larger than " \
            + "1, but smaller than the number of samples specified for " \
            + "other nodes. Unsure how to parse values."


class NoPlatingAllowedError(Exception):
    det_name: str
    message: str
    node_pgm_name: str

    def __init__(self,
                 det_name: str,
                 problematic_node_pgm_name: str,
                 message: str = "") -> None:
        self.det_name = det_name
        self.message = message
        self.node_pgm_name = problematic_node_pgm_name
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: When executing " + self.det_name + "(), replicates " \
            + "were detected for argument " + self.node_pgm_name \
            + ". Plating is not supported for this deterministic function."


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
        return "ERROR: When going through values in " \
            + self.node_name + ", the number of replicates " \
            + "differed among simulations. One or more of the " \
            + "simulations failed to generate the specified " \
            + "number of replicates. Exiting..."


class DimensionalityWarning(Exception):
    rv_name: str
    dn_name: str
    message: str

    def __init__(self, rv_name: str, dn_name: str, message: str = "") -> None:
        self.rv_name = rv_name
        self.dn_name = dn_name
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "Warning: Distribution " + self.dn_name \
            + " will be called to sample " + self.rv_name \
            + " using vectorization."


class MissingStateDependentParameterError(Exception):
    epoch_missing_param: int
    symmetric_diff_set: ty.Set[ty.Any]
    symmetric_diff_str: str
    message: str

    def __init__(self,
                 epoch_missing_param: int,
                 symmetric_diff_set: ty.Set[ty.Any],
                 message: str = "") -> None:
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


class RepeatedStateDependentParameterError(Exception):
    epoch_w_repeated_param: int
    repeated_state: ty.Tuple[int]
    repeated_state_str: str
    message: str

    def __init__(self,
                 epoch_w_repeated_param: int,
                 repeated_state: ty.Union[int, ty.Tuple[int]],
                 message: str = "") -> None:
        self.epoch_w_repeated_param = epoch_w_repeated_param
        self.repeated_state = repeated_state
        self.repeated_state_str = ""
        self.message = message

        if isinstance(repeated_state, int):
            self.repeated_state_str = str(repeated_state)

        elif isinstance(repeated_state, tuple):
            self.repeated_state_str = "("

            for i in repeated_state:
                ith_str = str(i)

                if self.repeated_state_str == "(":
                    self.repeated_state_str += ith_str

                else:
                    self.repeated_state_str += ", " + ith_str

            self.repeated_state_str += ")"

        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: State-dependent parameter defined by states " \
            + self.repeated_state_str + " were repeated in epoch " \
            + str(self.epoch_w_repeated_param)


class StateDependentParameterMisspec(Exception):
    message: str

    def __init__(self, message: str = "") -> None:
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


class InvalidMCMCChainLength(Exception):
    message: str

    def __init__(self, message: str) -> None:
        self.message = message
        super().__init__(self.message)

    def __str__(self) -> str:
        return "ERROR: " + self.message


# Initialization errors
class ObjInitInvalidArgError(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str, par_name: str, message: str = "") -> None:
        self.dn_name = dn_name
        self.message = "Could not initialize the object of \'" + dn_name + \
            "\'. The argument provided to " + par_name + " was invalid."
        if message:
            self.message += " " + message
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitMissingParameterError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, par_name: str, message: str = "") -> None:
        self.obj_name = obj_name
        self.message = "Could not initialize the object of " + obj_name + \
            " . Parameter " + par_name + " was not provided an argument."
        if message:
            self.message += " " + message
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ObjInitIncorrectDimensionError(Exception):
    obj_name: str
    message: str

    def __init__(self,
                 obj_name: str,
                 container_name: str,
                 obs_len: int,
                 exp_len: int = 0) -> None:
        self.obj_name = obj_name
        self.obs_len = obs_len
        self.exp_len = exp_len
        self.message = "Could not initialize the object of " + obj_name \
            + ". Incorrect dimension of container " + container_name \
            + ", which was of size " + str(obs_len) + ". "
        
        if exp_len != 0:
            self.message += "The expected dimension was " + str(exp_len) \
            + "."
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message

# Syntax errors: sampling distribution, det. function exceptions #
class ParseFunctionArgError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, message: str) -> None:
        self.par_name = par_name
        self.message = "Could not set value for function parameter " \
            + self.par_name + ". " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return  self.message


class ParseNotAParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str) -> None:
        self.par_name = par_name
        self.message = "\'" + par_name + "\' is not a valid parameter."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message
    

class ParseRequireSingleValueError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
                + self.arg + ("\', more than one value was "
                              "provided. A scalar is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequireIntegerError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
                + self.arg + ("\', something other than an integer was "
                              "provided. An integer is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseRequireNumericError(Exception):
    obj_name: str
    message: str

    def __init__(self, obj_name: str, arg: str) -> None:
        self.obj_name = obj_name
        self.arg = arg
        self.message = \
            "When specifying object " + obj_name + "'s parameter \'" \
                + self.arg + ("\', something other than a number was "
                              "provided. An number is required.")
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingParameterError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str) -> None:
        self.par_name = par_name
        self.message = "Parameter \'" + par_name + "\' is missing."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingArgumentError(Exception):
    par_name: str
    message: str

    def __init__(self, par_name: str, arg_name: str = "") -> None:
        self.par_name = par_name
        self.message = "Argument "
        
        if arg_name:
            self.message += "\'" + arg_name + "\' "
            
        self.message += "for \'" + par_name + "\' is missing."
        
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseMissingSpecificationError(Exception):
    message: str

    def __init__(self, obj2spec_name: str) -> None:
        self.message = "Missing specification when calling \'" \
            + obj2spec_name + "\'."
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseDnInitFailError(Exception):
    dn_name: str
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.dn_name = dn_name
        self.message = "Parsing the specification of \'" + self.dn_name \
            + "\' failed. " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ParseDetFnInitFailError(Exception):
    det_fn_name: str
    message: str

    def __init__(self, det_fn_name: str, message: str) -> None:
        self.det_fn_name = det_fn_name
        self.message = "Deterministic output from \'" + det_fn_name \
            + "\' could not be instantiated. " + message
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


# PGM exceptions #
class NodePGMNodeStatCantFloatError(Exception):
    message: str

    def __init__(self, node_name: str) -> None:
        self.message = "ERROR: When summarizing " + node_name \
            + " into string, a node statistic could not be float-ified"

    def __str__(self) -> str:
        return self.message


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
        return "ERROR: AnnotatedTree cannot be annotated this way. " \
            + self.message

# Generation exceptions #
class GenerateFailError(Exception):
    message: str

    def __init__(self, dn_name: str, message: str) -> None:
        self.message = "ERROR: Could not generate a sample from " \
            + dn_name + ". " + message
        
    def __str__(self) -> str:
        return self.message