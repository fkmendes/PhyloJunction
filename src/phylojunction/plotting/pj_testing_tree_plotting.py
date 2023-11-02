# Time-heterogeneous GeoSSE tree
#
# (((nd10:0.05585634266420793[&state=1],(nd12:0.01444225652017525[&state=0],nd13:0.01444225652017525[&state=2])nd11:0.04141408614403268[&state=2])nd2:1.225238932333245[&state=2],((nd6:0.15152443528855317[&state=1],nd7:0.15152443528855317[&state=1])nd4:0.16841576897032884[&state=1],(nd8:0.11388809491478649[&state=1],nd9:0.0982030908340934[&state=1])nd5:0.2060521093440955[&state=1])nd3:0.9611550707385709[&state=1])root:0.7189047250025469[&state=0])origin:0.0[&state=0];
# at_dict = 
# nd3    Attribute ('state') transition:
#     Time: 1.6028010851528143
#     Subtending node: nd3
#     Departing state: 0
#     Arriving state: 2,
#        Attribute ('state') transition:
#     Time: 1.6714902163855656
#     Subtending node: nd3
#     Departing state: 2
#     Arriving state: 1
# nd2    Attribute ('state') transition:
#     Time: 1.884598041236693
#     Subtending node: nd2
#     Departing state: 0
#     Arriving state: 2
# nd10    Attribute ('state') transition:
#     Time: 1.944143657335792
#     Subtending node: nd10
#     Departing state: 2
#     Arriving state: 1
# nd12    Attribute ('state') transition:
#     Time: 1.9855577434798246
#     Subtending node: nd12
#     Departing state: 2
#     Arriving state: 0

import phylojunction.readwrite.pj_read as pjr
import phylojunction.data.attribute_transition as pjat
import phylojunction.data.tree as pjt

import dendropy as dp
import matplotlib
import matplotlib.pyplot as plt  # type: ignore


if __name__ == "__main__":
    # tr_str = "(((nd10:0.05585634266420793[&state=1],(nd12:0.01444225652017525[&state=0],nd13:0.01444225652017525[&state=2])nd11:0.04141408614403268[&state=2])nd2:1.225238932333245[&state=2],((nd6:0.15152443528855317[&state=1],nd7:0.15152443528855317[&state=1])nd4:0.16841576897032884[&state=1],(nd8:0.11388809491478649[&state=1],nd9:0.0982030908340934[&state=1])nd5:0.2060521093440955[&state=1])nd3:0.9611550707385709[&state=1])root:0.7189047250025469[&state=0])origin:0.0[&state=0];"
    tr_str = "(((nd10[&state=1]:0.05585634266420793,(nd12[&state=0]:0.01444225652017525,nd13[&state=2]:0.01444225652017525)nd11[&state=2]:0.04141408614403268)nd2[&state=2]:1.225238932333245,((nd6[&state=1]:0.15152443528855317,nd7[&state=1]:0.15152443528855317)nd4[&state=1]:0.16841576897032884,(nd8[&state=1]:0.11388809491478649,nd9[&state=1]:0.0982030908340934)nd5[&state=1]:0.2060521093440955)nd3[&state=1]:0.9611550707385709)root[&state=0]:0.7189047250025469)origin[&state=0]:0.0;"

    ann_tr = pjr.read_nwk_tree_str(
        tr_str,
        "read_tree",
        in_file=False)
    
    at_dict = \
        {"nd3":
        [
            pjat.AttributeTransition(
                "state",
                "nd3",
                1.6028010851528143,
                0,
                2),
            pjat.AttributeTransition(
                "state",
                "nd3",
                1.6714902163855656,
                2,
                1)
        ],
        "nd2": [pjat.AttributeTransition("state", "nd2", 1.884598041236693, 0, 2)],
        "nd10": [pjat.AttributeTransition("state", "nd10", 1.944143657335792, 2, 1)],
        "nd12": [pjat.AttributeTransition("state", "nd12", 1.9855577434798246, 2, 0)]
        }

    ann_tr.at_dict = at_dict
    ann_tr.populate_nd_attr_dict(["state"], read_as_newick_str=True)
    
    fig = matplotlib.pyplot.figure()

    ax = fig.add_axes([0.25, 0.2, 0.5, 0.6])
    ax.patch.set_alpha(0.0)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    pjt.plot_ann_tree(
        ann_tr,
        ax,
        use_age=False,
        start_at_origin=True,
        sa_along_branches=False,
        attr_of_interest="state")
    # sa_along_branches=False

    plt.show()