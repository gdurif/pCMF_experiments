### wrapper for pca function

pca_call <- function(X, ncomp, centerX=TRUE, scaleX=FALSE) {
    
    X <- as.matrix(X)
    X <- apply(X, c(1,2), as.double)
    
    res <- prcomp(X, center=centerX, scale.=scaleX)
    V <- as.matrix(res$rotation)
    U <- as.matrix(res$x)
    
    exp_dev <- cumsum(res$sdev^2)/sum(res$sdev^2)
    
    return(list(U=U, V=V, exp_dev=exp_dev, centerX=centerX, scaleX=scaleX, ncomp=ncol(U)))
    
}
