import typing as ty

class SampledAncestor():

    label: str
    lineage_node_label: str # label of node on whose subtending branch the SA is attached
    global_time: float # SA's (forward) occurrence time (origin or root = 0.0)
    
    # time to lineage node (node defining the branch on which the SA is placed);
    # this class member is filled out (if unfilled) when the lineage node either undergoes an
    # event, or if the tree reaches a stop condition
    time_to_lineage_node: float

    def __init__(self, label: str, lineage_node_label: str, global_time: float, time_to_lineage_node: float=-1.0):
        self.label = label
        self.lineage_node_label = lineage_node_label
        self.global_time = global_time
        self.time_to_lineage_node = time_to_lineage_node
