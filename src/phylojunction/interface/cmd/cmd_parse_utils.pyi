import pgm.pgm as pgm
import typing as ty
from _typeshed import Incomplete

cannot_start_with_this_regex: Incomplete
cannot_end_with_this_regex: Incomplete
whitespace_regex: Incomplete
int_or_float_regex: Incomplete
vector_value_regex: Incomplete
assign_regex: Incomplete
character_value_regex: Incomplete
quoted_character_value_regex: Incomplete
sampled_as_regex: Incomplete
sampling_dn_spec_regex: Incomplete
deterministic_regex: Incomplete

def val_or_obj(pgm_obj: pgm.ProbabilisticGraphicalModel, val: ty.List[str]) -> ty.List[ty.Union[pgm.NodePGM, str]]: ...
def parse_spec(pgm_obj: pgm.ProbabilisticGraphicalModel, fn_spec_str: str, cmd_line: str) -> ty.Tuple[ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]], ty.List[pgm.NodePGM]]: ...
def parse_val_vector(vec_str: str) -> ty.List[str]: ...
def tokenize_fn_spec(fn_spec_str: str, cmd_line: str) -> ty.Dict[str, str]: ...
