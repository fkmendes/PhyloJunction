import typing as ty
import numpy as np
import pandas as pd  # type: ignore
import pickle  # type: ignore
import string
import csv

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec


def read_text_file(fp_string: str) -> ty.List[str]:
    """Read and parse text file into list of strings (one per line)

    Args:
        fp_string (str): String containing file path to text file being
            read

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
        fp_string (str): String containing file path to binary file storing
            PGM from previous PJ session

    Returns:
        pgm.ProbabilisticGraphicalModel: Probabilistic graphical model
            object to be initialized and returned
    """

    pgm_obj: pgm.ProbabilisticGraphicalModel

    with open(fp_string, "rb") as picklefile:
        pgm_obj = pickle.load(picklefile)

    return pgm_obj


def read_csv_into_dataframe(fp_string: str) -> pd.DataFrame:
    """Read .csv file into a pandas DataFrame

    Args:
        fp_string (str): String containing file path to text file being
            read

    Returns:
        pd.DataFrame: pandas DataFrame object (empty DataFrame if not
            CSV)
    """

    if is_csv(fp_string):
        return pd.read_csv(fp_string)

    else:
        return pd.DataFrame()


def is_csv(fp_string: str) -> bool:
    """Check if file in provided path is in CSV format

    Args:
        fp_string (str): String containing file path to text file being
            read

    Returns:
        bool: If file is in CSV format
    """

    try:
        with open(fp_string, newline="") as csvfile:
            start = csvfile.read(4096)  # reading file in 4 KiB chunks

            # isprintable does not allow newlines, printable does not allow
            # umlauts...
            if not all([c in string.printable or c.isprintable()
                        for c in start]):
                print(("Inside readwrite.pj_read.is_csv():\n    ERROR: "
                       "Attempted to read CSV file, but it does not seem "
                       "to be in CSV format."))

                return False

            dialect = csv.Sniffer().sniff(start)

            return True

    except csv.Error:
        # could not get a csv dialect -> probably not a csv.
        print(("Inside readwrite.pj_read.is_csv():\n    ERROR: Attempted to "
               "read CSV file, but it does not seem to be in CSV format."))

        return False


def parse_cli_str_write_fig(str_write_fig: str) \
        -> ty.Dict[str, ty.Tuple[int]]:
    """
    """

    node_range_dict: ty.Dict[str, ty.Tuple[int]] = dict()

    node_range_list = str_write_fig.split(";")

    # only one node name passed, no range
    # so we assume just the first repl of
    # first sample
    if len(node_range_list) == 1:
        if node_range_list[0].isnumeric():
            raise ec.PJCLIInvalidInput(
                "-f",
                ("If no range is provided, you must input a node name "
                 "(string) whose figure to plot, instead of an integer.")
            )
            # raise RuntimeError

        node_range_dict[node_range_list[0]] = tuple([0])

    else:
        node_names_str, node_ranges_str = node_range_list
        node_names_list = node_names_str.split(",")
        node_ranges_list = node_ranges_str.split(",")

        if len(node_names_list) != len(node_ranges_list):
            raise ec.PJCLIInvalidInput(
                "-f",
                ("If ranges are provided, the number of ranges must "
                 "match the number of node names to plot figures for.")
            )
            # raise RuntimeError

        for idx, node_name in enumerate(node_names_list):
            range_str = node_ranges_list[idx]

            try:
                range_tup = tuple(int(i) for i in range_str.split("-"))

            except ValueError:
                raise ec.PJCLIInvalidInput(
                    "-f",
                    ("Ranges must be defined by integers."))
                # raise RuntimeError

            if node_name in node_range_dict:
                raise ec.PJCLIInvalidInput(
                    "-f",
                    ("Node names appeared more than once. Nodes to plot "
                     "figures for must be unique.")
                )
                # raise RuntimeError

            node_range_dict[node_name] = range_tup

    return node_range_dict


if __name__ == "__main__":

    # should be ok
    print(parse_cli_str_write_fig("tr"))
    print(parse_cli_str_write_fig("tr;0"))
    print(parse_cli_str_write_fig("tr;0-10"))
    print(parse_cli_str_write_fig("rv1,tr;0,0"))
    print(parse_cli_str_write_fig("rv1,tr;0-1,0"))
    print(parse_cli_str_write_fig("rv1,tr;0-1,0-1"))

    # should raise hell (try one at a time)
    # parse_cli_str_write_fig("0")
    # parse_cli_str_write_fig("rv;tr")
    # parse_cli_str_write_fig("rv,tr;0")
    # parse_cli_str_write_fig("rv,tr;0,")