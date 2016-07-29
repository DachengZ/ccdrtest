## compare graphs and calculate SHD, fdr, ...

compare.edgeList <- function(estL, trueL, o = NULL) {
    ## edgeList from sparsebnUtils
    ## i %in% edL[[j]] means there is an edge i -> j
    nn <- length(trueL)
    if(length(estL) != nn) stop("Two edgeLists have different number of nodes.")
    if(!is.null(o)) {
        if(any(sort(o) != 1:nn)) stop("Wrong permute order!")
        else {
            ## o is the new order of nodes
            ## i is at position q[i] in vector o
            ## so edge i->j in the new graph is q[i] -> q[j] in the original graph
            q <- order(o)
            ## nodes0 <- names(edL) ## names(edL) is NULL
            tempL <- vector("list", length = nn)
            for(i in 1:nn) {
                tempL[[q[i]]] <- q[trueL[[i]]]
            }
            structure(tempL, class = c("edgeList", "list"))
            trueL <- tempL
        }
    }

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