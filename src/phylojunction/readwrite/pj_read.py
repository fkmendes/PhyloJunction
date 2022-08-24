import typing as ty
import numpy as np
import pandas as pd # type: ignore
import pickle # type: ignore
import string
import csv

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


def read_csv_into_dataframe(fp_string: str) -> pd.DataFrame:
    """Read .csv file into a pandas DataFrame

    Args:
        fp_string (str): String containing file path to text file being read

    Returns:
        pd.DataFrame: pandas DataFrame object (empty DataFrame if not CSV)
    """
    if is_csv(fp_string):
        return pd.read_csv(fp_string)

    else:
        return pd.DataFrame()


def is_csv(fp_string: str) -> bool:
    """Check if file in provided path is in CSV format

    Args:
        fp_string (str): String containing file path to text file being read

    Returns:
        bool: If file is in CSV format
    """
    try:
        with open(fp_string, newline="") as csvfile:
            start = csvfile.read(4096) # reading file in 4 KiB chunks

            # isprintable does not allow newlines, printable does not allow umlauts...
            if not all([c in string.printable or c.isprintable() for c in start]):
                print("Inside readwrite.pj_read.is_csv():\n    ERROR: Attempted to read CSV file, but it does not seem to be in CSV format.")
                return False
            dialect = csv.Sniffer().sniff(start)
            
            return True
    
    except csv.Error:
        # could not get a csv dialect -> probably not a csv.
        print("Inside readwrite.pj_read.is_csv():\n    ERROR: Attempted to read CSV file, but it does not seem to be in CSV format.")
        return False 