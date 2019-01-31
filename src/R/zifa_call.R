## wrapper for ZIFA python package

zifa_call <- function(X, ncomp, ncore=NULL) {
    
    if(!is.null(ncore)) {
        require(RhpcBLASctl)
        blas_set_num_threads(ncore)
        omp_set_num_threads(ncore)
    }
    
    require(reticulate)
    
    warning("If you use raw counts, you should use a log2 transform")
    
    X <- as.matrix(X)
    X <- apply(X, c(1,2), as.double)
    
    ####
    use_virtualenv(file.path(PYENV, "zifa"))
    ZIFA <- import("ZIFA")
    
    tmp <- ZIFA$block_ZIFA$fitModel(X, as.integer(ncomp))
    names(tmp[[1]]) <- c("U", "params")
    
    ####
    res <- list(U=tmp[[1]]$U, 
                V=tmp[[1]]$params$A,
                kept_genes=tmp[[2]])
    
    #### RETURN
    return(res)
}
