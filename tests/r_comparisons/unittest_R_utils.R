library(phytools)

get.tr.h <- function(a.tr) {
    max(nodeHeights(a.tr))
}

get.specs <- function(a.tr, w.states=TRUE, w.sa=FALSE, w.root=FALSE, three.states=FALSE) {
    extant.tip.labels = getExtant(a.tr)
    n.tips = length(extant.tip.labels) # extant!
    n.total = length(a.tr$tip.label)
    if (w.sa) {
        n.sa = get.n.sa(a.tr)
        n.total = n.total - n.sa
    }

    ## root age for all tips, fossil and extant
    tr.h = 0.0

    # starting from root
    if (w.root) {
        tr.h = max(nodeHeights(a.tr)[,2])
    }
    # starting from origin
    else {
        tr.h = get.tr.h(a.tr) # we have extant and extinct tips (n.total doesn't include SAs)

        if (w.sa) {
            sa.labels = list.sa(a.tr) # get SA labels
            pruned.tr = drop.tip(a.tr, a.tr$tip.label[match(sa.labels, a.tr$tip.label)])
            tr.h = get.tr.h(pruned.tr)
        }
    }

    if (w.states) {
        if (!three.states) {
            n.0 = sum(a.tr$tip.state==0) # both extant and extinct
            n.1 = sum(a.tr$tip.state==1)
            spec.df = data.frame(t(c(tr.h, n.total, n.tips, n.0, n.1)))
            names(spec.df) = c("tr.h", "n.total", "n.tips", "n.0", "n.1")
        } else {
            n.0 = sum(a.tr$tip.state==1) # both extant and extinct
            n.1 = sum(a.tr$tip.state==2)
            n.2 = sum(a.tr$tip.state==3)
            spec.df = data.frame(t(c(tr.h, n.total, n.tips, n.0, n.1, n.2)))
            names(spec.df) = c("tr.h", "n.total", "n.tips", "n.0", "n.1", "n.2")
        }
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
