library(phytools)

get.tr.h <- function(a.tr) {
    max(nodeHeights(a.tr))
}

get.specs <- function(a.tr, w.states=TRUE, w.sa=FALSE, w.root=FALSE) {
    extant.tip.labels = getExtant(a.tr)
    n.tips = length(extant.tip.labels) # extant!
    n.total = length(a.tr$tip.label)
    if (w.sa) {
        n.sa = get.n.sa(a.tr)
        n.total = n.total - n.sa
    }

    ## root age for all tips, fossil and extant
    tr.h = 0.0
    if (n.tips > 1) {
        tr.h = get.tr.h(a.tr) # if no SAs, just do this

    ## root age for all tips, fossil and extant, but we gotta remove the
    ## SAs because their implementation is hacky (they are tips
    ## and their parent nodes are not true speciation events, so those internal
    ## nodes cannot really be the root node; this is particularly important
    ## because SAs can also happen between the origin and the root, and thus be
    ## above the root
        if (w.sa) {
            if (n.tips > 1) {
                sa.labels = list.sa(a.tr) # get SA labels
                extant.fossil.tip.labels.no.sa = setdiff(a.tr$tip.label, sa.labels) # only non-SA tip labels
                root.node.nr = findMRCA(a.tr, extant.fossil.tip.labels.no.sa) # get MRCA node number of only non-SA tips
                node.heights.df = nodeHeights(a.tr) # all node heights

                ## weird case that seems ot happen occasionally
                ## if (root.node.nr > nrow(node.heights.df)) {
                ##     tr.h = max(node.heights.df[,2]) - nodeheight(a.tr, root.node.nr)
                ## } else {

                ## maximum distance between root and some node (has to be extant node)
                tr.h = max(dist.nodes(a.tr)[root.node.nr,])
            }
        }
    }

    else if (w.root) {
        tr.h = max(nodeHeights(a.tr)[,2])
    }

    ## else: just a single tip and starting at origin means root.age remains at  0.0

    if (w.states) {
        n.0 = sum(a.tr$tip.state==0 & names(a.tr$tip.state) %in% extant.taxa)
        n.1 = sum(a.tr$tip.state==1 & names(a.tr$tip.state) %in% extant.taxa)
        spec.df = data.frame(t(c(tr.h, n.total, n.tips, n.0, n.1)))
        names(spec.df) = c("tr.h", "n.total", "n.tips", "n.0", "n.1")
    }

    else {
        spec.df = data.frame(t(c(tr.h, n.total, n.tips)))
        names(spec.df) = c("tr.h", "n.total", "n.tips")
    }

    spec.df
}

list.sa <- function(a.tr) {

    epsilon = 1e-8

    ## first get the node numbers of the tips
    nodes = sapply(a.tr$tip.label, function(x,y) which(y==x), y=a.tr$tip.label)

    ## then get the edge lengths for those nodes
    edge.lengths = setNames(a.tr$edge.length[sapply(nodes,
                                                    function(x,y) which(y==x),y=a.tr$edge[,2])],names(nodes))

    names(edge.lengths[edge.lengths < 1e-8])
}

get.n.sa <- function(a.tr) {

    sa.names = list.sa(a.tr)

    length(sa.names) # n.sa
}
