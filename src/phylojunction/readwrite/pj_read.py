import os
import typing as ty
import numpy as np
import pandas as pd  # type: ignore
import pickle  # type: ignore
import dendropy as dp  # type: ignore
import string
import csv
import collections

# pj imports
import phylojunction.pgm.pgm as pgm
import phylojunction.utility.exception_classes as ec
import phylojunction.utility.helper_functions as pjh
from phylojunction.data.tree import AnnotatedTree

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def read_text_file(fp_string: str) -> ty.List[str]:
    """Read and parse text file into list of strings (one per line).

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
    """Read binary file storing PGM from a previous PJ session.

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
    """Read .csv file into a pandas DataFrame.

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
    """Check if file in provided path is in CSV format.

    Args:
        fp_string (str): String containing file path to text file being
            read

    Returns:
        bool: True if file is in CSV format, False otherwise
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
    """Parse command-line string argument for generating node plots.

    Args:
        str_write_fig (str): User-provided string as argument to command-line
            interface -f parameter (e.g., 'tr;0-10')
    
    Returns:
        Dictionary with node names as keys (str), tuple with range start
            (int, int)
    """

    node_range_dict: ty.Dict[str, ty.Tuple[int]] = dict()

    node_range_list = str_write_fig.split(";")

    # only one node name passed, we assume all samples
    # are needed, range tuple is empty
    if len(node_range_list) == 1:
        if node_range_list[0].isnumeric():
            raise ec.PJCLIInvalidInput(
                "-f",
                ("If no range is provided, you must input a node name "
                 "(string) whose figure to plot, instead of an integer.")
            )
            # raise RuntimeError

        node_range_dict[node_range_list[0]] = tuple([])

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
                raise ec.PJCLIInvalidInputError(
                    "-f",
                    ("Ranges must be defined by non-negative integers."))
                # raise RuntimeError

            if node_name in node_range_dict:
                raise ec.PJCLIInvalidInputError(
                    "-f",
                    ("Node names appeared more than once. Nodes to plot "
                     "figures for must be unique.")
                )
                # raise RuntimeError

            node_range_dict[node_name] = range_tup

    return node_range_dict


def read_nwk_tree_str(nwk_tree_path_or_str: str,
                      fn_name: str = "read_tree",
                      in_file: bool = True,
                      node_names_attribute: str = "",
                      epsilon: float = 1e-12) \
        -> AnnotatedTree:
    """Read Newick tree string directly or in provided file.

    Args:
        nwk_tree_path_or_str (str): Tree Newick string, or path to file
            containing tree Newick string
        in_file (bool): If tree string is in a file being passed as argument
            (True) or if Newick string is being passed directly (False).
            Defaults to 'True'.
        node_names_attribute (str): Defaults to empty string "".

    Returns:
        AnnotatedTree with populated attributes dictionary
            (i.e., a dendropy.Tree that has been annotated)
    """

    tr_str: str = str()
    node_attr_dict = collections.defaultdict(
        pjh.create_str_defaultdict)
    node_id_to_name: dict[str, str] = dict()

    # if tree provided in separate file
    if in_file:
        if not os.path.isfile(nwk_tree_path_or_str):
            raise ec.PJIOFileDoesNotExistError(fn_name,
                                               nwk_tree_path_or_str)
        
        with open(nwk_tree_path_or_str, "r") as infile:
            tr_str = infile.readline()

        dp_tr = dp.Tree.get(data=tr_str, schema="newick")

    # if tree is provided as newick string
    else:
        dp_tr = dp.Tree.get(data=nwk_tree_path_or_str, schema="newick")

    # first we deal with having a root vs an origin
    is_origin = (True if len(dp_tr.seed_node.child_nodes()) == 1 \
                 else False)
    
    root_height = dp_tr.max_distance_from_root()
    seen_root = False

    # debugging
    print("root height = " + str(root_height))
    print("is_origin = " + str(is_origin))

    # now we annotate nodes
    for nd in dp_tr.preorder_node_iter():
        nd.is_sa_lineage = False
        nd.is_sa_dummy_parent = False

        # (1) origin and/or root
        if nd is dp_tr.seed_node:
            nd.is_sa_dummy_parent = False
            nd.is_sa = False
            nd.alive = False
            
            nd_name = ("origin" if is_origin else "root")
            nd.taxon = dp.Taxon(label=nd_name)
            nd.label = nd_name
            
            if nd_name == "root":
                seen_root = True
        
        # (2) not origin
        else:
            if not nd.is_leaf():
                # initial values
                nd.is_sa_dummy_parent = False
                nd.is_sa = False
                nd.alive = False

                # seeing if origin child is root
                # or the dummy parent of a sampled ancestor
                for ch_nd in nd.child_node_iter():
                    if abs(ch_nd.edge_length) <= epsilon:
                        nd.is_sa_dummy_parent = True
                        nd.is_sa = False
                        break

                # root!
                if not nd.is_sa_dummy_parent and \
                        not seen_root:
                    nd.taxon = dp.Taxon(label="root")
                    nd.label = "root"
                    seen_root = True

            # annotate node as sampled ancestor or not
            if abs(nd.edge_length) <= epsilon:
                nd.is_sa_dummy_parent = False
                nd.is_sa = True

            else:
                nd.is_sa = False
        
            # annotate node as alive or not
            if abs(nd.distance_from_root() - root_height) <= epsilon:
                nd.alive = True

            else:
                nd.alive = False

            # annotate as sa_lineage
            if nd.parent_node.is_sa_dummy_parent:
                nd.is_sa_lineage = True

            # now we deal with translating between node
            # ids and node names
            nd_id = nd.annotations[node_names_attribute].value
            nd_name = "nd" + nd_id

            if nd.is_leaf():
                nd.label = nd.taxon.label

            # internal and no name provided
            elif nd.taxon is None:                
                if nd.label == "None":
                    nd.label = nd_name
                    nd.taxon = dp.Taxon(label=nd_name)

                else:
                    nd.taxon = dp.Taxon(label=nd.label)

                dp_tr.taxon_namespace.add_taxon(nd.taxon)

            node_id_to_name[nd_id] = nd.label
            node_attr_dict[nd.label][node_names_attribute] = nd_id

    # debugging
    # print(dp_tr.taxon_namespace)

    # for nd in dp_tr.postorder_node_iter():
    #     print(nd is dp_tr.seed_node)
    #     print(nd.parent_node)

    return AnnotatedTree(
            dp_tr,
            3,
            start_at_origin=is_origin,
            tree_died=False,
            read_as_newick_string=True)

    # print(dp_tr.taxon_namespace)
    # print(node_attr_dict)
    # print(node_id_to_name)


if __name__ == "__main__":

    # should be ok
    # print(parse_cli_str_write_fig("tr"))
    # print(parse_cli_str_write_fig("tr;0"))
    # print(parse_cli_str_write_fig("tr;0-10"))
    # print(parse_cli_str_write_fig("rv1,tr;0,0"))
    # print(parse_cli_str_write_fig("rv1,tr;0-1,0"))
    # print(parse_cli_str_write_fig("rv1,tr;0-1,0-1"))

    # should raise hell (try one at a time)
    # parse_cli_str_write_fig("0")
    # parse_cli_str_write_fig("rv;tr")
    # parse_cli_str_write_fig("rv,tr;0")
    # parse_cli_str_write_fig("rv,tr;0,")

    # should all work!
    tr = read_nwk_tree_str("examples/trees_maps_files/turtle.tre",
    # tr = read_nwk_tree_str("examples/trees_maps_files/dummy_tree1.tre",
    # tr = read_nwk_tree_str("examples/trees_maps_files/dummy_tree2.tre",
    # tr = read_nwk_tree_str("examples/trees_maps_files/dummy_tree3.tre",
                      node_names_attribute="index")
    
    # dummy trees should be:
    # (((sp1:1.0,sp2:1.5)nd4:1.0,sp3:1.7)root:0.5)origin:0.5
    # ((sp1:1.0,sp2:1.5)nd4:1.0,sp3:1.7)root
    # (((sp1:1.0,sp2:1.5)root:1.0,sp3:0.0)nd2:0.5)origin:0.5