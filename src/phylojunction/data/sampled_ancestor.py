import typing as ty

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class SampledAncestor():

    label: str
    # label of node on whose subtending branch the SA is attached
    lineage_node_label: str
    # SA's (forward) occurrence time (i.e., origin or root = 0.0)
    global_time: float
    state: int

    # time to lineage node (node defining the branch on which the
    # SA is placed); this class member is filled out (if unfilled)
    # when the lineage node either undergoes an event, or if the
    # tree reaches a stop condition
    time_to_lineage_node: float
    str_representation: str

    def __init__(self,
                 label: str,
                 lineage_node_label: str,
                 global_time: float,
                 state: int = 0,
                 time_to_lineage_node: float = -1.0) -> None:

        self.label = label
        self.lineage_node_label = lineage_node_label
        self.global_time = global_time
        self.state = state
        self.time_to_lineage_node = time_to_lineage_node
        self.str_representation = \
            "Sampled ancestor \'" + label + "\':\n    Branch: \'" \
            + lineage_node_label + "\'\n    Time: " + str(global_time)

    def __str__(self) -> str:
        if self.time_to_lineage_node < 0:
            self.str_representation += \
                "\n    Above \'" + self.label + "\': " \
                + str(self.time_to_lineage_node)

        return self.str_representation

    def __repr__(self) -> str:
        if self.time_to_lineage_node < 0:
            self.str_representation += \
                "\n    Above \'" + self.label + "\': "\
                + str(self.time_to_lineage_node)

        return self.str_representation
