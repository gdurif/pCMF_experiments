### log2 transform: log2(count+1)
# element-wise
log2_coord <- function(x) {
    return(log2(x+1))
}
# matrix-wise
log2m <- function(X, ncore=1) {
    require(parallel)
    cl <- makeCluster(ncore)
    res <- parApply(cl, X, 2, log2_coord)
    stopCluster(cl)
    return(res)
}

# tranform to integer matrix-wise
float2int <- function(X, ncore=1) {
    require(parallel)
    cl <- makeCluster(ncore)
    res <- parApply(cl, X, 2, as.integer)
    stopCluster(cl)
    return(res)
}
