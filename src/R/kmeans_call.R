#### functions for CLUSTERING

# MODEL: k-means
# ALGO: k-means

kmeans_call <- function(data, data_labels, k, iterMax=100, nstart=1, mySeed=NULL) {
    
    ## LIBRAIRIES
    require(mclust)
    
    ## INIT
    set.seed(mySeed)
    
    ## RUNNNING
    comp_time <- system.time(
        cluster_out <- kmeans(data, centers=k, 
                              iter.max=iterMax, nstart=nstart))[3]
    
    ## RandIndex
    adjRI <- adjustedRandIndex(cluster_out$cluster, as.integer(data_labels))
    
    ## within sum of square
    wSS <- cluster_out$tot.withinss
    
    ## BIC
    # bic <- BIC(cluster_out) # from the stat package, use logLik.kmeans
    bic <- BIC2(cluster_out)
    
    ## return
    return(list(cluster_out=cluster_out, adjRI=adjRI, wSS=wSS, bic=bic, Ncluster=k, comp_time=comp_time))
}


BIC2 <- function(fit) {
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    # return(data.frame(AIC = D + 2*m*k,
    #                   BIC = D + log(n)*m*k))
    return(D + log(n)*m*k)
}

logLik.kmeans <- function(object) structure(
    object$tot.withinss,
    df = nrow(object$centers)*ncol(object$centers),
    nobs = length(object$cluster)
)
