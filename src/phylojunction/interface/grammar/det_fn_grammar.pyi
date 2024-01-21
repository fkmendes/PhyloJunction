import phylojunction.pgm.pgm as pgm
import phylojunction.calculation.discrete_sse as sseobj
import typing as ty

class PJDetFnGrammar:
    det_fn_grammar_dict: ty.Dict[str, ty.Set[str]]
    @classmethod
    def grammar_check(cls, det_fn_id: str, fn_param: str) -> bool: ...
    @classmethod
    def init_return_state_dep_rate(cls, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) -> sseobj.DiscreteStateDependentRate: ...
    @classmethod
    def init_return_sse_stash(cls, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) -> ty.Tuple[sseobj.MacroevolEventHandler, sseobj.DiscreteStateDependentProbabilityHandler]: ...
    @classmethod
    def create_det_fn_obj(cls, det_fn_id: str, det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodeDAG]]]) -> ty.Optional[ty.Union[sseobj.DiscreteStateDependentRate, sseobj.MacroevolEventHandler]]: ...
