import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import typing as ty

class PJDetFnGrammar:
    det_fn_grammar_dict: ty.Dict[str, ty.Set[str]]
    @classmethod
    def grammar_check(cls, det_fn_id: str, fn_param: str) -> bool: ...
    @classmethod
    def init_return_sse_rate(cls, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) -> sseobj.MacroevolStateDependentRateParameter: ...
    @classmethod
    def init_return_macroevol_handler(cls, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) -> sseobj.MacroevolEventHandler: ...
    @classmethod
    def create_det_fn_obj(cls, det_fn_id: str, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) -> ty.Optional[ty.Union[sseobj.MacroevolStateDependentRateParameter, sseobj.MacroevolEventHandler]]: ...
