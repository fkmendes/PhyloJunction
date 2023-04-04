import typing as ty

# pj imports
import phylojunction.calculation.discrete_sse as sseobj
import phylojunction.utility.exception_classes as ec
import phylojunction.interface.grammar.det_fn_discrete_sse as detsse
import phylojunction.pgm.pgm as pgm

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"

class PJDetFnGrammar():

    det_fn_grammar_dict: ty.Dict[str, ty.Set[str]]

    ##########################################
    #  All available deterministic functions #
    ##########################################
    det_fn_grammar_dict = {
        "sse_prob": set(["name", "value", "state"]),
        "sse_rate": set(["name", "value", "event", "states"]),
        "sse_wrap": set(["flat_rate_mat", "flat_prob_mat", "n_states", "seed_age", "epoch_age_ends", "n_epochs"])
    }

    @classmethod
    def grammar_check(cls, det_fn_id: str, fn_param: str) -> bool:
        if fn_param in cls.det_fn_grammar_dict[det_fn_id]:
            return True
        return False

    @classmethod
    def init_return_state_dep_rate(cls,
        det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> sseobj.DiscreteStateDependentRate:
        
        return detsse.make_DiscreteStateDependentRate(det_fn_param_dict)

    @classmethod
    def init_return_state_dep_prob(cls,
        det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> sseobj.DiscreteStateDependentProbability:
        
        return detsse.make_DiscreteStateDependentProbability(det_fn_param_dict)

    @classmethod
    def init_return_sse_stash(cls,
        det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
        -> sseobj.SSEStash:
        
        return detsse.make_SSEStash(det_fn_param_dict)

    @classmethod
    def create_det_fn_obj(cls,
        det_fn_id: str,
        det_fn_param_dict: ty.Dict[str, ty.List[ty.Union[str, pgm.NodePGM]]]) \
            -> ty.Optional[
                ty.Union[
                    sseobj.DiscreteStateDependentRate,
                    sseobj.DiscreteStateDependentProbability,
                    sseobj.SSEStash
                ]]:

        # validate input
        if det_fn_id == "sse_rate":
            return cls.init_return_state_dep_rate(det_fn_param_dict)

        if det_fn_id == "sse_prob":
            return cls.init_return_state_dep_prob(det_fn_param_dict)

        if det_fn_id == "sse_wrap":
            for arg in det_fn_param_dict:
                if not cls.grammar_check(det_fn_id, arg):
                    raise ec.NotAParameterError(arg)
            return cls.init_return_sse_stash(det_fn_param_dict)

        return None

if __name__ == "__main__":
    # can be called from interface/grammar/
    # $ python3 det_fn_grammar.py
    # 
    # can also be called from phylojunction/
    # $ python3 interface/grammar/det_fn_grammar.py
    # or
    # $ python3 -m interface.grammar.det_fn_grammar
    #
    # can also be called from VS Code, if open folder is phylojuction/

    print(PJDetFnGrammar.grammar_check("sse_rate", "name")) # True
    print(PJDetFnGrammar.grammar_check("sse_rate", "cloud")) # False

    # TODO: add code to run init_return_state_dep_rate and init_return_sse_stash