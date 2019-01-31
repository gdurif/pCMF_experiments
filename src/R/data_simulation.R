#### generate data for simulations


data_generation <- function(n, p, K, ZI=0, nblockU, nblockV, epsilonU, espilonV, p0, reorder=TRUE) {
    
    ### zero-inflation
    if(ZI==0) {
        prob1 <- rep(1,p)
    } else if(ZI==1) {
        prob1 <- round(runif(p, 0.4, 0.6), digits=2)
    } else if(ZI==2) {
        prob1 <- round(runif(p, 0.2, 0.4), digits=2)
    }
    
    ## organization of U
    signalBlockU <- matrix(rep(4, each=nblockU*nblockU), nrow=nblockU, ncol=nblockU)
    epsilonBlockU <- diag(1, nblockU, nblockU)
    signalBlockU <- signalBlockU + epsilonU * epsilonBlockU
    blockAlpha1 <- blockMatrix(nrow=n, ncol=K, nRowBlock=nblockU, nColBlock=nblockU, signalBlock=signalBlockU)
    alpha1 <- blockAlpha1$mat
    alpha2 <- matrix(1, nrow=n, ncol=K)
    
    ## organization of V
    signalBlockV <- matrix(rep(4, each=nblockV*nblockV), nrow=nblockV, ncol=nblockV)
    epsilonBlockV <-  diag(1, nblockV, nblockV)
    signalBlockV <- signalBlockV + epsilonV * epsilonBlockV
    blockBeta1 <- blockMatrix(nrow=p0, ncol=K, nRowBlock=nblockV, nColBlock=nblockV, signalBlock=signalBlockV)
    beta1 <- NULL
    if(p-p0 > 0) {
        beta1 <- rbind(blockBeta1$mat, matrix(0.7, nrow=p-p0, ncol=K))
    } else {
        beta1 <- blockBeta1$mat
    }
    beta2 <- matrix(1, nrow=p, ncol=K)
    
    ## generating the data
    data <- dataGeneration(n=n, p=p, K=K, alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2, 
                           ZI=(ZI!=0), prob1=prob1, reorder=reorder)
    
    ## return
    res <- list(data=data, 
                indivLabels=blockAlpha1$idRows[data$orderInd],
                varLabels=blockBeta1$idRows[data$orderVar],
                blockAlpha1=blockAlpha1, blockBeta1=blockBeta1)
    return(res)
}

shift <- function(x) {return(c(tail(x, length(x)-1), head(x,1)))}
compos <- function(x, n, fun) {
    res <- matrix(NA, nrow=length(x), ncol=n)
    res[,1] = x
    if(n>=2) {
        for(i in 2:n) {
            res[,i] = fun(res[,i-1])
        }
    }
    return(res)
}

custom_generate_factor_matrix <- function(nfeature, K, ngroup=1, average_signal=1,
                                          group_separation=0.5,
                                          distribution="uniform",
                                          shuffle_feature=TRUE,
                                          prop_noise_feature=0,
                                          noise_level=0.5,
                                          seed=NULL) {
    # input check
    if(!is.numeric(nfeature) | !(nfeature>0) | (floor(nfeature) - nfeature !=0)) {
        stop("'nfeature' input parameter should be a positive integer")
    }
    if(!is.numeric(K) | !(K>0) | (floor(K) - K != 0)) {
        stop("'K' input parameter should be a positive integer")
    }
    if(!is.numeric(ngroup) | !(ngroup>0) | (ngroup>K)) {
        stop("'ngroup' input parameter should be a positive integer <= 'K'")
    }
    if(!is.numeric(average_signal) | (average_signal<0)) {
        stop("'average_signal' input parameter should be a real positive value")
    }
    if(!is.numeric(group_separation) | (group_separation<0) | (group_separation>1)) {
        stop("'group_separation' input parameter should be a real value in [0,1]")
    }
    if(!(distribution %in% c("uniform", "gamma"))) {
        stop("'distribution' input parameter should be \"uniform\" or \"gamma\"")
    }
    if(!is.logical(shuffle_feature)) {
        stop("'shuffle_feature' input parameter should be a boolean")
    }
    if(!is.numeric(prop_noise_feature) | (prop_noise_feature<0) | (prop_noise_feature>1)) {
        stop("'prop_noise_feature' input parameter should be a real value in [0,1]")
    }
    if(!is.numeric(noise_level) | (noise_level<0)) {
        stop("'noise_level' input parameter should be a real positive value")
    }
    if(!is.null(seed)) {
        if(!is.integer(seed) | (seed<0) ) {
            stop("'seed' input parameter should be NULL or a positive integer")
        }
    }
    
    # set seed if necessary
    if(!is.null(seed)) {
        set.seed(seed)
    }
    
    # generate matrix of parameter for each block, dim ngroup x ngroup
    param_block_matrix <- diag(average_signal, nrow=ngroup, ncol=ngroup)
    if(ngroup>1) {
        param_block_matrix <- param_block_matrix + average_signal * (1 - group_separation) * (1 - diag(1, nrow=ngroup, ncol=ngroup))
    }
    # shuffle groups (columns of 'param_block_matrix')
    col_order = sample.int(n=ngroup, size=ngroup, replace=FALSE)
    param_block_matrix <- param_block_matrix[,col_order]
    
    # block id
    id_block <- 1:ngroup
    
    ### structuring block matrix
    nfeature0 <- floor( (1 - prop_noise_feature) * nfeature)
    mat0 <- NULL
    id_rows <- NULL
    
    if(nfeature0 > 0) {
        
        param_matrix <- NULL
        
        if(ngroup > 1) {
            # assign features (=rows) to a block
            id_rows <- sort(rep(id_block, length=nfeature0))
            
            # assign components (=cols) to a block
            id_cols <- sort(rep(id_block, length=K))
            
            # construction of the parameter matrix by block, dim nfeature x K
            param_matrix <- matrix(NA, nrow=nfeature0, ncol=K)
            for(row_block in id_block) {
                for(col_block in id_block) {
                    rows_in_block <- (1:nfeature0)[id_rows == row_block]
                    nrow_block <- length(rows_in_block)
                    cols_in_block <- (1:K)[id_cols == col_block]
                    ncol_block <- length(cols_in_block)
                    
                    param_matrix[rows_in_block, cols_in_block] <- matrix(param_block_matrix[row_block, col_block], nrow=nrow_block, ncol=ncol_block)
                }
            }
        } else {
            id_rows <- rep(1, nfeature0)
            param_matrix <- matrix(average_signal, nrow=nfeature0, ncol=K)
        }
        
        # generation of the factor sub-matrix
        if(distribution == "uniform") {
            mat0 <- sapply(1:K, function(k) return(runif(nfeature0, min=0, max=2*param_matrix[,k]/sqrt(K)))) # matrix nfeature0 x K
        } else if(distribution == "gamma") {
            mat0 <- sapply(1:K, function(k) return(rgamma(nfeature0, shape=param_matrix[,k]/sqrt(K), rate=1))) # matrix nfeature0 x K
        }
        
    }
    
    ### noise block
    mat_noise <- NULL
    nfeature_noise <- nfeature - nfeature0
    if(nfeature_noise > 0) {
        # generation of the factor sub-matrix
        if(distribution == "uniform") {
            mat_noise <- sapply(1:K, function(k) return(runif(nfeature_noise, min=0, max=2*noise_level/sqrt(K)))) # matrix nfeature_noise x K
        } else if(distribution == "gamma") {
            mat_noise <- sapply(1:K, function(k) return(rgamma(nfeature_noise, shape=noise_level/sqrt(K), rate=1))) # matrix nfeature_noise x K
        }
    }
    
    ### concatenate structuring block matrix and noise block matrix
    mat <- rbind(mat0, mat_noise)
    
    # shuffle individuals
    feature_order <- 1:nfeature
    if(shuffle_feature) {
        feature_order <- sample.int(n=nfeature, size=nfeature, replace=FALSE)
    }
    mat <- mat[feature_order,]
    feature_label <- c(id_rows, rep(0, nfeature_noise))[feature_order]
    
    ### output
    res <- list(factor_matrix=mat, feature_order=feature_order, feature_label=feature_label)
    return(res)
    
}