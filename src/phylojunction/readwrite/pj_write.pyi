import typing as ty

import pandas
import pandas as pd  # type: ignore
import matplotlib.pyplot as plt  # type: ignore

# pj imports
import phylojunction.pgm.pgm as pgm

def write_text_output(outfile_handle: ty.IO, content_string_list: ty.List[str]) -> None: ...
def write_data_df(outfile_handle: ty.IO, data_df: pd.DataFrame, format: str = ...) -> None: ...
def write_fig_to_file(outfile_path: str,fig_obj: plt.Figure) -> None: ...
def prep_data_df(sample_size: int, node_dag_list: ty.List[pgm.NodeDAG], write_nex_states: bool=False) -> ty.Tuple[ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]], ty.List[ty.Union[ty.Dict[str, pd.DataFrame], ty.Dict[str, str]]]]: ...
def prep_trees_rb_smap_dfs(dag_obj: pgm.DirectedAcyclicGraph, tree_dag_node_name_list: ty.List[str], attr_being_mapped: str) -> ty.Tuple[ty.Dict[str, ty.List[pandas.DataFrame]], ty.Dict[str, ty.List[str]]]: ...
def prep_data_filepaths_dfs(scalar_output_stash: ty.List[ty.Union[pd.DataFrame, ty.Dict[int, pd.DataFrame]]], tree_output_stash: ty.List[ty.Union[ty.Dict[str, pd.DataFrame], ty.Dict[str, str]]] = []) -> ty.Tuple[ty.List[str], ty.List[ty.Union[pd.DataFrame, str]]]: ...
def dump_pgm_data(dir_string: str, dag_obj: pgm.DirectedAcyclicGraph, prefix: str, write_nex_states: bool = ...) -> None: ...
def dump_trees_rb_smap_dfs(dir_string: str, dag_obj: pgm.DirectedAcyclicGraph, tr_dag_node_name_list: ty.List[str], mapped_attr_name: str, prefix: str = "") -> None: ...
def dump_serialized_pgm(dir_string: str, dag_obj: pgm.DirectedAcyclicGraph, prefix: str = ...) -> None: ...
def get_write_inference_rev_scripts(all_sims_model_spec_list: ty.List[str], all_sims_mcmc_logging_spec_list: ty.List[str], dir_list: ty.List[str], prefix: str = ..., write2file: bool = ...) -> ty.List[str]: ...
