source("sparsebnconvert.R")
### compare with the original
### compareGraphs is wrong?
### compareGraphs(ccdr2graph(ccdr.path[[3]]), g1)

# compare.graph <- function(g.est, g.true, o = NULL) {
#   m.est <- t(as(g.est, "matrix") != 0)
#   m.true <- t(as(g.true, "matrix") != 0)
#   return(compare.adj.matrix(m.est, m.true, o))
# }

compare.graph <- function(g.est, g.true, o = NULL) {
    if(!is.null(o)) g.true <- permutenodes(g.true, o)
    trueL <- sapply(edgeL(g.true), getElement, "edges")
    estL <- sapply(edgeL(g.est), getElement, "edges")
    nn <- length(trueL)
    ## do not use edgeL. it refers to node names

    P <- TP <- R <- FP <- 0
    for(j in 1:nn) {
        lj <- length(estL[[j]])
        P <- P + lj ## P: number of estimated (predicted) edges
        if(lj != 0) for(i in estL[[j]]) {
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

compare.cFg <- function(cF, g, o = NULL) {
    gcF <- ccdrFit2graph(cF)
    return(compare.graph(cF, g, o))
}

### do we need to order back?
### compare.adj.matrix(mcF, m0)
### compare.adj.matrix(mcF, m0[o, o])
### compare.adj.matrix(mcF[q, q], m0)
### and why different from
### compareGraphs(ccdrFit2graph(ccdr.path[[14]]), g)
