rmvDAG.fix <- function(dag, n, vfix = NULL, m = 0, s = 1) {
    ## function for generating random data based on DAG
    ## Currently only supports single-node intervention in each sample
    ## To do: to support ivnlist

    ## Based on rmvDAG
    ## Original Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler
    ## ----------------------------------------------------------------------

    ## check input &  initialize variables
    stopifnot(is(dag, "graph"),
              (p <- length(nodes(dag))) >= 2)

    if(!is.null(vfix)) if(length(vfix) != n) stop("wrong length of vfix!")

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
    if(is.null(vfix)) vfix <- rep(p + 1, n) ## todo: change to 0?
    if (sum(weightMatrix) > 0) {
        X <- errMat
        for(i in 1:n) {
            v <- vfix[i]
            if(v == 1) X[i, 1] <- X[i, 1] * s + m
            for(j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
                ## in case input vfix is not ordered
                if(v == j) X[i, j] <- X[i, j] * s + m
                else {
                    ij <- 1:(j-1)
                    X[i, j] <- X[i, j] + X[i, ij] %*% weightMatrix[j, ij]
                }
            }
        }
        return(X)
    }
    else return(errMat)
} # END RMVDAG.FIX
