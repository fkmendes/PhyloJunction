from enum import Enum, EnumMeta

def exp_root_height_yule_ntaxa(birth_rate: float, n_taxa: int) -> float: ...
def exp_extant_count_bd(birth_rate: float, death_rate: float, tree_age: float, speciation_k: int = ...) -> float: ...

class MetaEnum(EnumMeta):
    def __contains__(cls, item): ...

class BaseEnum(Enum, metaclass=MetaEnum): ...

class ParametricDistribution(str, BaseEnum):
    EXPONENTIAL: str
    LOGNORMAL: str
    NORMAL: str
    GAMMA: str
    UNIFORM: str
