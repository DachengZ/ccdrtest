library(ccdrAlgorithm)
setwd("~/ccdrtest") ## change this accordingly

library(bnlearn)
library(graph)
## library(RBGL)
library(pcalg)

source("ccdr-rmvDAG.R") ## generate random data with optional intervention
source("ccdr-convert.R") ## consider use functions from sparsebnUtils packages
source("ccdr-compare.R") ## compare graphs. different from `compareGraphs`

### Test function wrapper
ccdrtest <- function(nn, pp, num.edges, itvtimes = NULL, vfix = NULL, gamma = 2, lambdas.length = 20) {
    # Goal: Based on a same graph, we generata two sets of data,
    #       one observational, one with intervention
    #       And test CCDr, CCDri, PC on them.

    ## Parameters for the DAG
    # nn: number of samples to draw for observational data
    # pp: number of many nodes in the DAG
    # num.edges: number of *expected* edges in the DAG
    ss <- num.edges / pp # the expected number of parents *per node*

    ## To record performance
    metric <- matrix(0, 4, 8) ## to record performance metrics

    ### Generate a random DAG using the pcalg method randomDAG
    beta.min <- 0.5
    beta.max <- 1
    edge.pr <- 2 * ss / (pp - 1)
    g <- randomDAG(n = pp, prob = edge.pr, lB = beta.min, uB = beta.max) # Note that the edge weights are selected at random here!

    ### Observational data
    vfix0 <- as.integer(rep(pp + 1, nn))
    dat <- rmvDAG.fix(dag = g, n = nn, vfix = vfix0)

    ## Permute columns
    o <- sample(1:pp)
    o1 <- c(o, pp + 1)
    q <- order(o) ## o[q] == q[o] == 1:pp
    q1 <- order(o1)
    dat1 <- dat[, o] # permute the columns to randomize node ordering
    g1 <- permutenodes.graph(g, o)
    edL1 <- edgeL(g1)
    ## permutations actually affects estimates

    ## pack as sparsebnData
    suppressMessages(sData.null <- sparsebnUtils::sparsebnData(dat1, type = "continuous"))
    ## CCDrA
    ptm <- proc.time()
    ccdrA.path <- ccdr.run(data = sData.null, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdri.path)
    compare.path <- sapply(lapply(ccdrA.path, getElement, "edges"), compare.sFg, edL1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    metric[1, 1:7] <- compare.path[, z]
    metric[1, 8] <- (proc.time() - ptm)[3]
    # metric[1, 8:10] <- (proc.time() - ptm)[1:3]

    ## MMHC from `bnlearn`
    # PC from `pcalg` gives bi-directed DAG # Why?
    ptm <- proc.time()
    # pc.fit <- pc(suffStat = list(C = cor(dat1), n = nn1), indepTest = gaussCItest, alpha = 0.01, p = pp, verbose = F)
    mmhc.fit <- mmhc(as.data.frame(dat1))
    mmhc.graph <- as.graphNEL(mmhc.fit)
    metric[2, 1:7] <- compare.graph(mmhc.graph, g1)
    metric[2, 8] <- (proc.time() - ptm)[3]
    # metric[2, 8:10] <- (proc.time() - ptm)[1:3]

    true.edges <- metric[1, 2] - metric[1, 4] + metric[1, 7]

    ### Intervention data
    if(is.null(vfix)) {
        if(is.null(itvtimes)) itvtimes <- floor(nn / pp)
        vfix <- rep(sample(pp), itvtimes)
    } else itvtimes <- paste("~", length(vfix) / pp, sep = "")
    # it seems that, here itvtimes can be non-integer;
    # R will replace itvtimes by floor(itvtimes)
    nn1 <- length(vfix)
    dat <- rmvDAG.fix(dag = g, n = nn1, vfix = vfix)
    dat1 <- dat[, o] # permute the columns to randomize node ordering
    vfix1 <- as.list(q1[vfix])
    g1 <- permutenodes.graph(g, o)
    edL1 <- edgeL(g1)
    ## permutations actually affects estimates

    ## pack as sparsebnData
    sData.vfix <- sparsebnUtils::sparsebnData(dat1, type = "continuous", ivn = vfix1)

    ## CCDrA
    ptm <- proc.time()
    ccdrA.path <- ccdr.run(data = sData.vfix, gamma = gamma, lambdas.length = lambdas.length, alpha = 10, verbose = FALSE)
    # print(ccdri.path)
    compare.path <- sapply(lapply(ccdrA.path, getElement, "edges"), compare.sFg, edL1)
    shd.val <- compare.path[7, ]
    z <- which(shd.val == min(shd.val))
    z <- z[length(z)]
    metric[3, 1:7] <- compare.path[, z]
    metric[3, 8] <- (proc.time() - ptm)[3]
    # metric[3, 8:10] <- (proc.time() - ptm)[1:3]

    ## MMHC from `bnlearn`
    # PC from `pcalg` gives bi-directed DAG # Why?
    ptm <- proc.time()
    # pc.fit <- pc(suffStat = list(C = cor(dat1), n = nn1), indepTest = gaussCItest, alpha = 0.01, p = pp, verbose = F)
    mmhc.fit <- mmhc(as.data.frame(dat1))
    mmhc.graph <- as.graphNEL(mmhc.fit)
    metric[4, 1:7] <- compare.graph(mmhc.graph, g1)
    metric[4, 8] <- (proc.time() - ptm)[3]
    # metric[4, 8:10] <- (proc.time() - ptm)[1:3]

    namelist <- rep(c("CCDrA", "MMHC"), 2)
    nnlist <- c(rep(nn, 2), rep(length(vfix), 2))
    colnames(metric) <- c("P", "TP", "R", "FP", "TPR", "FDR", "SHD", "elapsed")
    itvlist <- c(rep("obs", 2), rep("itv", 2))
    itvtimes <- c(rep(NA, 2), rep(itvtimes, 2)) # use 0 instead of NA?
    return(data.frame(name = namelist, nn = nnlist, pp, num.edges,
                      itv = itvlist, itvtimes = itvtimes, true.edges,
                      metric))
}

m <- vector("list", 50 * 15)

# specify either
# 1) itvtimes: sample(pp) and replicate, so intervention times per node is the same
# 2) vfix: a vector of node labels, so that each row has intervention on one node.

for(i in    1:50) m[[i]] <- ccdrtest(nn = 20, pp = 10, num.edges =  5, itvtimes = 2)
for(i in  51:100) m[[i]] <- ccdrtest(nn = 20, pp = 10, num.edges = 10, itvtimes = 2)
for(i in 101:150) m[[i]] <- ccdrtest(nn = 20, pp = 10, num.edges = 20, itvtimes = 2)
for(i in 151:200) m[[i]] <- ccdrtest(nn = 50, pp = 20, num.edges = 10, itvtimes = 2)
for(i in 201:250) m[[i]] <- ccdrtest(nn = 50, pp = 20, num.edges = 20, itvtimes = 2)
for(i in 251:300) m[[i]] <- ccdrtest(nn = 50, pp = 20, num.edges = 40, itvtimes = 2)
for(i in 301:350) m[[i]] <- ccdrtest(nn = 100, pp = 50, num.edges = 25, itvtimes = 2)
for(i in 351:400) m[[i]] <- ccdrtest(nn = 100, pp = 50, num.edges = 50, itvtimes = 2)
for(i in 401:450) m[[i]] <- ccdrtest(nn = 100, pp = 50, num.edges = 100, itvtimes = 2)
for(i in 451:500) m[[i]] <- ccdrtest(nn = 200, pp = 100, num.edges = 50, itvtimes = 2)
for(i in 501:550) m[[i]] <- ccdrtest(nn = 200, pp = 100, num.edges = 100, itvtimes = 2)
for(i in 551:600) m[[i]] <- ccdrtest(nn = 200, pp = 100, num.edges = 200, itvtimes = 2)
for(i in 601:650) m[[i]] <- ccdrtest(nn = 500, pp = 200, num.edges = 100, itvtimes = 2)
for(i in 651:700) m[[i]] <- ccdrtest(nn = 500, pp = 200, num.edges = 200, itvtimes = 2)
for(i in 701:750) m[[i]] <- ccdrtest(nn = 500, pp = 200, num.edges = 400, itvtimes = 2)

saveRDS(m, file = "ccdrtest.rds")

# To average over one case:
m7 <- cbind(m[[301]][, 1:6], Reduce('+', lapply(m[301:350], '[', , 7:15)) / 50); m7
m8 <- cbind(m[[351]][, 1:6], Reduce('+', lapply(m[351:400], '[', , 7:15)) / 50); m8
m9 <- cbind(m[[401]][, 1:6], Reduce('+', lapply(m[401:450], '[', , 7:15)) / 50); m9

# getm <- function(row, column) {
#     return(sapply(m, "[", row, column))
# }
#
# ## Check one simulation:
# # m[[5]]
#
# pdf(file = "obs.pdf")
#     par(oma = rep(0, 4), mar = c(4, 4, 2, 1))
#     plot(getm(1, 15), type = "p", pch = 1, main = "Timing Observational Data", xlab = "test index", ylab = "Time elapsed (seconds)", ylim = c(0, 5))
#     points(getm(1, 15), type = "l", lty = 1)
#     points(getm(2, 15), type = "p", pch = 20)
#     points(getm(2, 15), type = "l", lty = 2)
#     points(getm(3, 15), type = "p", pch = 24)
#     points(getm(3, 15), type = "l", lty = 3)
#     legend(1, 5, legend = c("CCDr", "CCDri", "MMHC"), lty = c(1, 2, 3), pch = c(1, 20, 24))
# dev.off()
#
# pdf(file = "itv.pdf")
#     par(oma = rep(0, 4), mar = c(4, 4, 2, 1))
#     plot(getm(4, 15), type = "p", pch = 1, main = "Timing Interventional Data", xlab = "test index", ylab = "Time elapsed (seconds)", ylim = c(0, 5))
#     points(getm(4, 15), type = "l", lty = 1)
#     points(getm(5, 15), type = "p", pch = 20)
#     points(getm(5, 15), type = "l", lty = 2)
#     points(getm(6, 15), type = "p", pch = 24)
#     points(getm(6, 15), type = "l", lty = 3)
#     legend(1, 5, legend = c("CCDr", "CCDri", "MMHC"), lty = c(1, 2, 3), pch = c(1, 20, 24))
# dev.off()
