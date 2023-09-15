import typing as ty
import pandas as pd

# pj imports
import phylojunction.pgm.pgm as pgm

def read_text_file(fp_string: str) -> ty.List[str]: ...
def read_serialized_pgm(fp_string: str) -> pgm.ProbabilisticGraphicalModel: ...
def read_csv_into_dataframe(fp_string: str) -> pd.DataFrame: ...
def is_csv(fp_string: str) -> bool: ...
def parse_cli_str_write_fig(str_write_fig: str) -> ty.Dict[str, ty.Tuple[int]]: ...