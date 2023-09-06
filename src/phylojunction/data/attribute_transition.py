import typing as ty

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class AttributeTransition():

    attr_label: str
    from_state: int
    to_state: int
    subtending_node_label: str  # parent can be dummy node
    # (forward) occurrence time of attribute change
    # (i.e., origin or root = 0.0)
    global_time: float
    str_representation: str

    def __init__(self,
                 attr_label: str,
                 subtending_node_label: str,
                 global_time: float,
                 from_state: int,
                 to_state: int) -> None:

        self.attr_label = attr_label
        self.subtending_node_label = subtending_node_label
        self.global_time = global_time
        self.from_state = from_state
        self.to_state = to_state

        self.str_representation = \
            "Attribute (\'" + attr_label + "\') transition:" + \
            "\n    Time: " + str(global_time) + \
            "\n    Subtending node: " + subtending_node_label + \
            "\n    Departing state: " + str(from_state) + \
            "\n    Arriving state: " + str(to_state)

    def update_daughter_members(self,
                                daughter_node_label: str,
                                daughter_node_time: float) -> None:

        if not self.daughter_node_label:
            # TODO: do stuff
            pass

        # we only need to set daughter members once
        else:
            pass

    def __str__(self) -> str:
        return self.str_representation
