import typing as ty

__author__ = "Fabio K. Mendes"
__email__ = "f.mendes@wustl.edu"


class AttributeTransition():

    attr_label: str
    from_state: int
    to_state: int  # if happening at speciation, child 1
    to_state2: int  # if happening at speciation, child 2
    at_speciation: bool

    # (i) if object stored in an AnnotatedTree's at_dict, this is
    # the node subtending a branch (note that this transition may be
    # one that happen right at speciation, but it is nonetheless
    # 'framed' as an anagenetic transition happening along a branch
    # for the purposes of painting branches when tree plotting)
    #
    # (ii) if object is stored in an AnnotatedTree's clado_at_dict,
    # this node is the speciating node; the clado_at_dict is not used
    # for plotting, but instead assists when tabulating and parsing
    # stochastic maps
    #
    # (parent can be dummy node)
    subtending_or_speciating_node_label: str

    # (forward) occurrence time of attribute change
    # (i.e., origin or root = 0.0)
    global_time: float
    str_representation: str

    def __init__(self,
                 attr_label: str,
                 subtending_node_label: str,
                 global_time: float,
                 from_state: int,
                 to_state: int,
                 to_state2: ty.Optional[int] = None,
                 at_speciation: bool = False) -> None:

        self.attr_label = attr_label
        self.subtending_or_speciating_node_label = subtending_node_label
        self.global_time = global_time
        self.from_state = from_state
        self.to_state = to_state
        self.at_speciation = at_speciation

        self.str_representation = \
            "Attribute (\'" + attr_label + "\') transition:" + \
            "\n    Time: " + str(global_time) + \
            "\n    Subtending node: " + subtending_node_label + \
            "\n    Departing state: " + str(from_state) + \
            "\n    Arriving state: " + str(to_state)

        # if second arriving state provided, this is a cladogenetic
        # attribute transition, we must adjust the string representation
        if to_state2 or at_speciation:
            self.str_representation = \
                self.str_representation.replace("Subtending", "Speciating")
            self.str_representation = \
                self.str_representation.replace("Departing", "Parent")
            self.str_representation = \
                self.str_representation.replace("Arriving", "Left child")
            self.str_representation += "\n    Right child state: " + \
                                       str(to_state2)

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
