### wrapper for spca

spca_call <- function(X, ncomp, centerX=TRUE, orthog=FALSE, 
                      verbose=TRUE, cv=TRUE, fit=TRUE, nfolds=5, 
                      lambda=NULL, nlambda=100, niter=40, niter.cv=20,
                      full_result=FALSE, ncore=NULL) {
    
    if(!is.null(ncore)) {
        require(RhpcBLASctl)
        blas_set_num_threads(ncore)
        omp_set_num_threads(ncore)
    }
    
    require(PMA)
    
    X <- as.matrix(X)
    p <- 1:ncol(X)
    
    if((cv==FALSE) && is.null(lambda)) {
        stop("Message from spca: error, if cross-validation is not processed, user should provide a real value for lambda")
    }
    
    if((!cv) && (!fit)) {
        warning("Message from spca: nor cross-validation neither fitting will be done")
    }
    
    ## tuning
    if(cv) {
        cv_out <- SPC.cv(X, 
                         sumabsvs=seq(1.2, sqrt(ncol(X)), len=nlambda), 
                         nfolds=5, 
                         niter=niter.cv, 
                         center=centerX,
                         orth=orthog, 
                         trace=verbose)
        
        lambda <- cv_out$bestsumabsv
    } else {
        cv_out <- NULL
        lambda <- NULL
    }
    
    ## fitting
    if(fit) {
        spca_out <- SPC(X, 
                        sumabsv=lambda, 
                        niter=niter, 
                        K=ncomp, 
                        cnames=colnames(X), 
                        center=centerX,
                        orth=orthog, 
                        trace=verbose, 
                        compute.pve=TRUE)
        
        U <- t(t(as.matrix(spca_out$u)) * spca_out$d)
        V <- t(t(as.matrix(spca_out$v)) * spca_out$d)
        
        exp_var <- spca_out$prop.var.explained
    } else {
        spca_out <- NULL
        U <- NULL
        V <- NULL
        exp_var <- NULL
    }
    
    res <- NULL
    if(full_result) {
        res <- list(U=U, V=V, lambda=lambda, exp_var=exp_var, 
                    spca_out=spca_out, cv_out=cv_out,
                    centerX=centerX, ncomp=ncol(U))
    } else {
        res <- list(U=U, V=V, lambda=lambda, exp_dev=exp_var,
                    centerX=centerX, ncomp=ncol(U))
    }
    
    return(res)
}
