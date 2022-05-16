import sys
sys.path.extend(["../", "../phylojunction"]) # necessary to run it as standalone on command line (from phylojunction/ or phylojunction/interface/)
import typing as ty
import re

# pj imports
import pgm.pgm as pgm
import utility.exception_classes as ec

# general
cannot_start_with_this_regex = re.compile('[\(\[\{|^&+\-%*/=!><_~0-9]')
cannot_end_with_this_regex = re.compile('[\(\[\{|^&+\-%*/=!><_~]$')
whitespace_regex = re.compile("\s+")
int_or_float_regex = re.compile("(\d+(?:\.\d+)?)")
vector_value_regex = re.compile("[\[\]]")

# assignment
assign_regex = re.compile("\s*(<-)\s*") # TODO: make sure only one '-' makes into this regex
character_value_regex = re.compile("\s*([a-zA-Z]+)\s*")
quoted_character_value_regex = re.compile("\s*(\"[a-zA-Z]+\")\s*")

# sampling dn
sampled_as_regex = re.compile("\s*(~)\s*")
sampling_dn_spec_regex = re.compile("([a-zA-Z]+_*[a-zA-Z]*)\((.+)\)")

# det functions
deterministic_regex = re.compile("\s*(:=)\s*")


def val_or_obj(pgm_obj: pgm.ProbabilisticGraphicalModel, val: ty.List[str]) -> ty.List[ty.Union[pgm.NodePGM, str]]:
    """_summary_

    Args:
        pgm_obj (pgm.ProbabilisticGraphicalModel): Probabilistic graphical model object
        val (str): List of strings, one for each argument coming from a parsed command

    Raises:
        ec.InexistentVariableError: _description_

    Returns:
        (str): List of strings (representing values or quoted-enclosed strings) and/or stochastic node objects
    """
    val_or_obj_list: ty.List[ty.Union[pgm.NodePGM, str]] = []
    
    for v in val:
        if type(v) == str:
            # checking if string could potentially be a node object that has a name
            if re.match(character_value_regex, v):
                # if it does find a node with that name, we add the node object
                try:
                    val_or_obj_list.append(pgm_obj.node_name_val_dict[v]) # appending StochasticNodePGM
                except:
                    raise ec.InexistentVariableError(v)

            # if string is instead number, or a word enclosed in quotes, we add the string
            else:
                val_or_obj_list.append(v)
        
    return val_or_obj_list


def parse_spec(pgm_obj: pgm.ProbabilisticGraphicalModel, fn_spec_str: str, cmd_line: str) -> ty.Tuple[ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]], ty.List[pgm.NodePGM]]:
    spec_dict: ty.Dict[str, str] = tokenize_fn_spec(fn_spec_str, cmd_line)
    spec_dict_return: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]

    parent_pgm_nodes = list()
    for param_name, an_arg in spec_dict.items():

        #############
        # IMPORTANT #
        #############
        arg_list = list() # arguments of dn/det parameters will come out as lists

        # if argument is a list
        if re.match(vector_value_regex, an_arg):
            # print("\n\n parsing arg " + str(an_arg) + " as vector")
            arg_list = parse_val_vector(an_arg)
        # if scalar variable
        else:
            # print("\n\n parsing arg " + str(an_arg) + " as str")
            arg_list.append(an_arg)
        
        val_obj_list = val_or_obj(pgm_obj, arg_list)
        spec_dict_return[param_name] = val_obj_list
        # { param_name_str: [ number_or_quoted_str_str1, a_NodePGM1, number_or_quoted_str_str2, aNodePGM2 ] }

        for vo in val_obj_list:
            if type(vo) == pgm.NodePGM:
                parent_pgm_nodes.append(vo)

    # values in spec_dict will be lists of strings or lists of NodePGMs
    return spec_dict_return, parent_pgm_nodes


def parse_val_vector(vec_str: str) -> ty.List[str]:
    """Parse string specifying vector when a variable is assigned

    Args:
        vec_str (str): A string specifying a vector, e.g., "[1, 2, 3]"

    Raises:
        ec.ScriptSyntaxError: If there are more than two squared brackets, this is thrown

    Returns:
        (str): Vector of strings, each being one of the values inside the vector being specified by \'vec_str\'
    """
    if len(re.findall(vector_value_regex, vec_str)) > 2:
        raise ec.ScriptSyntaxError(vec_str, "Something went wrong during variable assignment. If a vector of values is specified, there can only be one left and one right squared brackets. Exiting...")

    # ok!
    else:
        scalar_values_str = re.sub(whitespace_regex, "", re.sub(vector_value_regex, "", vec_str))
        return scalar_values_str.split(",")


def tokenize_fn_spec(fn_spec_str: str, cmd_line: str) -> ty.Dict[str, str]:
    """Tokenize string containing function specification into dictionary

    Args:
        fn_spec (str): String containing the specification of function (i.e., what appears between \'(\' and \')\'

    Returns:
        (dict): Dictionary containing function parameter strings as keys, argument strings as values (e.g., {"mean": "0.0"})
    """

    stop_token_capture_chs = tuple(["]", "="])

    spec_dict = dict()
    token, last_key_token = str(), str()
    first_ch = True
    reading_vector = False

    for ch in fn_spec_str:
        if first_ch and re.match(cannot_start_with_this_regex, ch):
            raise ec.ScriptSyntaxError(cmd_line, "Cannot start specification of object with special character. Exiting...")

        # ignore whitespace
        if re.match(whitespace_regex, ch): continue

        # finished reading token (could be parameter or argument)
        elif ch in stop_token_capture_chs:
            # token was parameter, we make it a key
            if ch == "=":
                spec_dict[token] = ""
                last_key_token = token
                
            # character closed a vector, we make it a value
            elif ch == "]":
                token += ch
                spec_dict[last_key_token] = token
                reading_vector = False # matters only when ch == "]"

            token = str()
        
        # if not stopping, keep adding
        else:
            if ch == ",":
                # finished reading a parameter=argument pair
                if not reading_vector and token:
                    spec_dict[last_key_token] = token
                    token = str()
                    continue

                # reading vector, we care about the comma
                if reading_vector:
                    token += ch
                    continue
            
            # keep adding as long as it's not a comma
            else:
                if ch == "[":
                    reading_vector = True
                
                token += ch

            first_ch = False
    
    # gotta wrap up and add last one, if it wasn't a vector
    # (vectors add themselves up when they finish no matter what)
    if token:
        spec_dict[last_key_token] = token

    return spec_dict

if __name__ == "__main__":
    # can be called from interface/
    # $ python3 cmd_parse_utils.py
    # 
    # can also be called from phylojunction/
    # $ python3 interface/cmd_parse_utils.py
    # or
    # $ python3 -m interface.cmd_parse_utils
    #
    # can also be called from VS Code, if open folder is phylojuction/
    
    fn_spec1 = "name=\"lambda0\", value=1.0, event=\"w_speciation\", states=[0, 0, 0]"
    cmd_line1 = "a <- sse_rate(" + fn_spec1 + ")"
    token_dict1 = tokenize_fn_spec(fn_spec1, cmd_line1)
    
    print(token_dict1)
