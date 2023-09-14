import math
from enum import Enum, EnumMeta

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


def exp_root_height_yule_ntaxa(birth_rate: float, n_taxa: int) -> float:
    """
    Return expected root height of Yule tree simulated up to a specified
    number of tips

    Args:
        birth_rate (float): Birth-rate (lambda) of Yule process.
        n_taxa (int): Tree is simulated until this number of tips.

    Returns:
        float: Expected root height
    """

    expected_root_height: float = 0.0

    # Gernhardt, T. (2008), J. Theor. Biol.
    k = 1  # first speciation event (root!)
    for i in range(k + 1, n_taxa + 1):
        expected_root_height += (1.0 / (float(i) * birth_rate))

    return expected_root_height


def exp_extant_count_bd(birth_rate: float,
                        death_rate: float,
                        tree_age: float,
                        speciation_k: int = 1) -> float:
    """Return expected count of extant taxa in a birth-death tree

    Args:
        birth_rate (float): Birth-rate (lambda) of birth-death process.
        death_rate (float): Death-rate (mu) of birth-death process.
        tree_age (float): Age of birth-death tree (starting from origin).
        speciation_k (int): Extant taxa are counted conditioning on process
            starting on k-th speciation event. Default: 1.

    Returns:
        int: Expected count of extant taxa.
    """

    k: int = speciation_k  # default: first speciation event (root!)

    if speciation_k:
        k = speciation_k

    expected_extant_count = \
        k * math.exp((birth_rate - death_rate) * tree_age)

    return expected_extant_count


class MetaEnum(EnumMeta):
    # Note to self: later read up on Python's metaclasses
    def __contains__(cls, item):
        try:
            cls(item)

        except ValueError:
            return False

        return True


class BaseEnum(Enum, metaclass=MetaEnum):
    pass


class ParametricDistribution(str, BaseEnum):
    EXPONENTIAL = "expn"
    LOGNORMAL = "lnormal"
    NORMAL = "normal"
    GAMMA = "gamma"
    UNIFORM = "uniform"
