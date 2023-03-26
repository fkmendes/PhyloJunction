import typing as ty

class AttributeTransition():

    attr_label: str
    from_state: int
    to_state: int
    
    # parent can be dummy node
    subtending_node_label: str
    # time_elapsed_since_parent_node: float
    
    # daughter can be dummy node
    # daughter_node_label: str
    # time_to_daughter_node: float

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
        
        self.str_representation = "Attribute (\'" + attr_label + "\') transition:" + \
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