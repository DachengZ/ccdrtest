## all kinds of conversions

sFedL2graphedL <- function(sFedL) {
    ## edgeList of sF: children, then parent.
    ## e.g. sF$edges[[2]] = 1 means edge 1->2 exists.
    ## edgeL of graphNEL: parent, then children.
    ## e.g. edL(g)[[1]]$edges = 2 means edge 1->2 exists.
    pp <- length(sFedL)
    edL <- vector("list", pp)
    for(i in 1:pp) edL[[i]] <- list(edges = NULL)
    for(i in 1:pp) {
        if(length(sFedL[[i]]) > 0) {
            for(j in sFedL[[i]]) edL[[j]]$edges <- c(edL[[j]]$edges, i)
        }
    }
    names(edL) <- as.character(1:pp)
    return(edL)
}

sparsebnFit2graphedL <- function(sF) {
    ## edgelist of sF: children, then parent.
    ## e.g. sF$edges[[2]] = 1 means edge 1->2 exists.
    ## edgelist of graphNEL: parent, then children.
    ## e.g. edL(g)[[1]]$edges = 2 means edge 1->2 exists.
    return(sFedL2graphedL(sF$edges))
}

sparsebnFit2graph <- function(sF) {
    ## edgelist of sF: children, then parent
    ## edgelist of graphNEL: parent, then children
    edL <- sparsebnFit2graphedL(sF)
    return(graph::graphNEL(nodes = names(edL), edgeL = edL, edgemode = 'directed'))
}

sparsebnFit2matrix <- function(sF) {
    ## to matrix
    ## m_{ij} = 1 <--> edge i->j exists
    ## or t(get.adjacency.matrix()) for a sparse matrix
    ## or wgtMatrix(, FALSE)
    return(as(sparsebnFit2graph(sF), "matrix") != 0)
}

graphedL2sFedL <- function(gedL) {
    ## edgelist of sF: children, then parent.
    ## e.g. sF$edges[[2]] = 1 means edge 1->2 exists.
    ## edgelist of graphNEL: parent, then children.
    ## e.g. edL(g)[[1]]$edges = 2 means edge 1->2 exists.
    pp <- length(gedL)
    sFedL <- vector("list", pp)
    for(j in 1:pp) sFedL[[j]] <- integer(0)
    for(i in 1:pp) {
        if(length(gedL[[i]]$edges) > 0) {
            for(j in gedL[[i]]$edges) sFedL[[j]] <- c(sFedL[[j]], i)
        }
    }
    structure(sFedL, class = c("edgeList", "list"))
}

## adjacency matrix to graph
## m_{ij} = 1 <--> edge i->j exists
## no longer need to transpose here
adj2graph <- function(m, newname = NULL) {
    g <- as(graphAM(adj = m, edgemode = 'directed'), "graphNEL")
    if(!is.null(newname)) nodes(g) <- newname
    return(g)
}

## permute nodes (and edges) of a graph
permutenodes.edgeL <- function(edL, o) {
    ## edgeL from graphNEL
    ## so that when plot, they are the same
    ## but look at the matrix, the rows and columns are permuted
    ## o is the new order of nodes
    ## i is at position q[i] in vector o
    ## so edge i->j in the new graph is q[i] -> q[j] in the original graph
    q <- order(o)
    nodes0 <- names(edL)
    nn <- length(nodes0)

    edL1 <- vector("list", length = nn)
    names(edL1) <- nodes0[o]
    for(i in 1:nn) {
        edL1[[q[i]]]$edges <- q[edL[[i]]$edges]
    }
    return(edL1)
}

permutenodes.edgeList <- function(edL, o) {
    ## edgeList from sparsebnUtils
    ## so that when plot, they are the same
    ## but look at the matrix, the rows and columns are permuted
    ## o is the new order of nodes
    ## i is at position q[i] in vector o
    ## so edge i->j in the new graph is q[i] -> q[j] in the original graph
    q <- order(o)
    ## nodes0 <- names(edL) ## names(edL) is NULL
    nn <- length(edL)

    edL1 <- vector("list", length = nn)
    ## names(edL1) <- nodes0[o]
    for(i in 1:nn) {
        edL1[[q[i]]] <- q[edL[[i]]]
    }
    structure(edL1, class = c("edgeList", "list"))
}

permutenodes.graph <- function(g, o) {
    ## so that when plot, they are the same
    ## but look at the matrix, the rows and columns are permuted
    ## o is the new order of nodes
    ## i is at position q[i] in vector o
    ## so edge i->j in the new graph is q[i] -> q[j] in the original graph
    q <- order(o)
    nodes0 <- nodes(g)
    edges0 <- edgeL(g)

    nn <- length(nodes0)
    nodes1 <- nodes0[o]

    edges1 <- vector("list", length = nn)
    names(edges1) <- nodes1
    for(i in 1:nn) {
        edges1[[q[i]]]$edges <- q[edges0[[i]]$edges]
    }
    return(graphNEL(nodes = nodes1, edgeL = edges1, edgemode = 'directed'))
}
