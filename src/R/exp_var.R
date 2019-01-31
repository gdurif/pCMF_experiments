### version 1
## parameters
# X : matrix n x p of observations
# U : matrix n x K of observation coordinates in the latent basis
## output
# explained variance
exp_var <- function(X, U) {
    require(MASS)
    X.cent <- scale(X, center=TRUE, scale=FALSE)
    X.proj <- U %*% ginv(t(U) %*% U) %*% t(U) %*% X
    X.proj.cent <- scale(X.proj, center=TRUE, scale=FALSE)
    return(sum(diag(t(X.proj.cent) %*% X.proj.cent)) / sum(diag(t(X.cent) %*% X.cent)))
}