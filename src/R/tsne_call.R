### wrapper for Rtsne

tsne_call <- function(X, ncomp, verbose, perplexity, ncore=NULL, ...) {
    
    if(!is.null(ncore)) {
        require(RhpcBLASctl)
        blas_set_num_threads(ncore)
        omp_set_num_threads(ncore)
    }
    
    require(Rtsne)
    
    X <- as.matrix(X)
    X <- apply(X, c(1,2), as.double)
    
    tmp <- Rtsne(X=X, dims=ncomp, 
                 perplexity=perplexity, 
                 verbose=verbose, ...)
    
    res <- list(U=tmp$Y, 
                perplexity=perplexity)
    
    #### RETURN
    return(res)
}
