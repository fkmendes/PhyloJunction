import typing as ty
import numpy as np
import pandas as pd # type: ignore
import pickle # type: itnore

# pj imports
import phylojunction.pgm.pgm as pgm

def read_text_file(fp_string: str) -> ty.List[str]:
    """Read and parse text file into list of strings (one per line)

    Args:
        fp_string (str): String containing file path to text file being read

    Returns:
        str: List of strings, each being a line of the input text file
    """
    cmd_line_list = list()
    
    with open(fp_string, "r") as infile:
        for line in infile:
            line = line.rstrip()

            if line: cmd_line_list.append(line)
    
    return cmd_line_list


def read_serialized_pgm(fp_string: str) -> pgm.ProbabilisticGraphicalModel:
    """Read binary file storing PGM from a previous PJ session

    Args:
        fp_string (str): String containing file path to binary file storing PGM from previous PJ session

    Returns:
        pgm.ProbabilisticGraphicalModel: Probabilistic graphical model object to be initialized and returned
    """
    pgm_obj: pgm.ProbabilisticGraphicalModel
    
    with open(fp_string, "rb") as picklefile:
        pgm_obj = pickle.load(picklefile)

    return pgm_obj