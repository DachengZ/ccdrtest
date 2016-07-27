## compare graphs and calculate SHD, fdr, ...

compare.edgeL <- function(edL.est, edL.true, o = NULL) {
    ## edgeL from graphNEL
    ## i %in% edgeL[[j]]$edges means there is an edge j -> i
    if(!is.null(o)) edL.true <- permutenodes.edgeL(edL.true, o)
    estL <- sapply(edL.est, getElement, "edges")
    trueL <- sapply(edL.true, getElement, "edges")
    nn <- length(trueL)
    ## do not use edgeL. it refers to node names

    P <- TP <- R <- FP <- 0
    for(j in 1:nn) {
        lj <- length(estL[[j]])
        P <- P + lj ## P: number of estimated (predicted) edges
        if(lj != 0) for(i in estL[[j]]) {
            ## so j -> i in estimated graph
            if(i %in% trueL[[j]]) TP <- TP + 1 ## TP: number of true positives
            else {
                if(j %in% trueL[[i]]) R <- R + 1 ## R: number of reversed edges
                else FP <- FP + 1 ## FP: number of false positives
            }
        }
    }
    Tedge <- sum(sapply(trueL, length))
    if(P == 0) fdr = 0 else fdr = (R + FP) / P
    ### Fedge <- pp * (pp - 1) / 2 - Tedge
    ### fpr = (R + FP) / Fedge
    ## SHD <- Tedge - TP + FP
    return(c(p = P, tp = TP, r = R, fp = FP, tpr = TP / Tedge, fdr = fdr, shd = Tedge - TP + FP))
}

compare.graph <- function(g.est, g.true, o = NULL) {
    return(compare.edgeL(edgeL(g.est), edgeL(g.true), o))
}

compare.edgeList <- function(estL, trueL, o = NULL) {
    ## edgeList from sparsebnUtils
    ## i %in% edL[[j]] means there is an edge i -> j
    if(!is.null(o)) trueL <- permutenodes.edgeList(trueL, o)
    nn <- length(trueL)

    P <- TP <- R <- FP <- 0
    for(j in 1:nn) {
        lj <- length(estL[[j]])
        P <- P + lj ## P: number of estimated (predicted) edges
        if(lj != 0) for(i in estL[[j]]) {
            ## so i -> j in the estimated graph
            if(i %in% trueL[[j]]) TP <- TP + 1 ## TP: number of true positives
            else {
                if(j %in% trueL[[i]]) R <- R + 1 ## R: number of reversed edges
                else FP <- FP + 1 ## FP: number of false positives
            }
        }
    }
    Tedge <- sum(sapply(trueL, length))
    if(P == 0) fdr = 0 else fdr = (R + FP) / P
    ### Fedge <- pp * (pp - 1) / 2 - Tedge
    ### fpr = (R + FP) / Fedge
    ## SHD <- Tedge - TP + FP
    return(c(p = P, tp = TP, r = R, fp = FP, tpr = TP / Tedge, fdr = fdr, shd = Tedge - TP + FP))
}

compare.sFg <- function(sFedL, gedL, o = NULL) {
    ## sFedL is edgeList from the estimated graph
    ## gedL is edgeL from the true graph
    if(!is.null(o)) gedL <- permutenodes.edgeL(gedL, o)
    return(compare.edgeList(sFedL, graphedL2sFedL(gedL)))
}
