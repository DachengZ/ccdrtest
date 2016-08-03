rmvDAG.fix <- function(dag, n, ivnlist = NULL, m = 0, s = 1) {
    ## function for generating random data based on DAG
    ## Currently only supports single-node intervention in each sample
    ## To do: to support ivnlist

    ## Based on rmvDAG
    ## Original Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler
    ## ----------------------------------------------------------------------

    ## check input &  initialize variables
    stopifnot(is(dag, "graph"),
              (p <- length(dag@nodes)) >= 2)

    if(!is.null(ivnlist)) if(length(ivnlist) != n) stop("wrong length of ivnlist!")

    ## as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
    ## weightMatrix <- t(as(dag,"matrix"))
    ## weightMatrix <- if(back.compatible) wgtMatrix.0(dag) else wgtMatrix(dag)
    weightMatrix <- wgtMatrix(dag)

    ## check if top. sorted
    nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
    if (nrow(nonZeros) > 0) {
        if (any(nonZeros[,1] - nonZeros[,2] < 0) || any(diag(weightMatrix) != 0))
            stop("Input DAG must be topologically ordered!")
    }

    ## compute X matrix X_i
    errMat <- matrix(rnorm(n * p), nrow = n)
    colnames(errMat) <- 1:p
    if(is.null(ivnlist)) ivnlist <- vector("list", n) ## todo: change to 0?
    if (sum(weightMatrix) > 0) {
        X <- errMat
        for(i in 1:n) {
            v <- ivnlist[[i]]
            if(1 %in% v) X[i, 1] <- X[i, 1] * s + m
            for(j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
                ## in case input vfix is not ordered
                if(j %in% v) X[i, j] <- X[i, j] * s + m
                else {
                    ij <- 1:(j-1)
                    X[i, j] <- X[i, j] + X[i, ij] %*% weightMatrix[j, ij]
                }
            }
        }
        return(X)
    }
    else return(errMat)
} # END rmvDAG.fix
