### function to run experiment

## experiment 1: clustering and visualization
experiment1 <- function(input) {
    
    if(!is.data.frame(input)) {
        stop("pb with input")
    } else {
        if(nrow(input)>1) stop("pb with input")
    }
    
    print(input)
    
    rep <- input$rep
    ncomp <- input$ncomp
    drop_out <- input$drop_out
    group_separation <- input$group_separation
    prop_noise_feature <- input$prop_noise_feature
    nrep <- input$nrep
    n <- input$n
    p <- input$p
    K <- input$K
    ngroup <- input$ngroup
    
    if(!is.null(SEED)) {
        set.seed(SEED)
    }
    
    factorU <- NULL
    factorV <- NULL
    count_data <- NULL
    X <- NULL
    keep_genes <- NULL
    for(count in 1:50) {
        average_signal <- c(250,100,250)
        factorU <- generate_factor_matrix(n, K, ngroup=ngroup, 
                                          average_signal=average_signal,
                                          group_separation=group_separation,
                                          distribution="exponential",
                                          shuffle_feature=FALSE,
                                          tag="cell")
        
        tmp1 <- 20*(1-prop_noise_feature)
        tmp2 <- 20*prop_noise_feature
        prop_noise_feature <- round(rbeta(1, 
                                          shape1 = runif(1, min=tmp2-2, max=tmp2+2), 
                                          shape2 = runif(1, min=tmp1-2, max=tmp1+2)), 
                                    digits = 3)
        
        factorV <- generate_factor_matrix(p, K, ngroup=2, average_signal=80,
                                          group_separation=0.8,
                                          distribution="exponential",
                                          shuffle_feature=TRUE,
                                          prop_noise_feature=prop_noise_feature,
                                          tag="gene")
        U <- factorU$factor_matrix
        V <- factorV$factor_matrix
        
        tmp1 <- 20*(1-drop_out)
        tmp2 <- 20*drop_out
        prob1 <- round(rbeta(p, 
                             shape1 = runif(1, min=tmp1-2, max=tmp1+2), 
                             shape2 = runif(1, min=tmp2-2, max=tmp2+2)), 
                       digits = 3)
        
        count_data <- generate_count_matrix(n, p, K, U, V,
                                            ZI=TRUE, prob1=prob1)
        X <- count_data$X
        Xmax <- apply(X,2,max)
        keep_genes <- (Xmax>0) & (Xmax<=quantile(Xmax, probs=0.95)) #& (apply(X,2,function(x) sum(x!=0))>=0.05*nrow(X))
        
        if(sum(keep_genes) > 0.5 * p) break
    }
    
    X <- X[,keep_genes]
    Xlog <- log2m(X)
    indiv_label <- as.factor(factorU$feature_label)
    feature_label <- as.factor(factorV$feature_label[keep_genes])
    
    additional_factor <- cbind(input[1,],
                               data.frame(prop0 = 1-sum(X!=0)/prod(dim(X))),
                               data.frame(ngene = ncol(X)),
                               data.frame(prop_noise = prop_noise_feature))
    
    ## run
    print(paste0("dim X = ", paste0(dim(X), collapse=", ")))
    
    tmp_data_to_plot_U <- NULL
    tmp_data_to_plot_V <- NULL
    tmp_result <- NULL
    
    if(!is.null(SEED)) {
        set.seed(NULL)
    }
    
    ### sparse pCMF
    filt_spcmf <- rep(TRUE, ncol(X))
    print("### pCMF (sparse)")
    
    # prior_S <- 1-exp(-apply(X,2,sd)/mean(X[X!=0]))
    # filt_spcmf <- prior_S > 0.2
    Xtmp <- X[,filt_spcmf]
    # print(paste0("dim X = ", paste0(dim(Xtmp), collapse=", "), " after pCMF pre-filtering"))
    
    prior_S <- 1-exp(-apply(Xtmp,2,sd)/mean(Xtmp[Xtmp!=0]))
    prob_S <- matrix(rep(prior_S, ncomp), ncol=ncomp)
    
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_sparse_gap_factor(Xtmp, K = ncomp, U = NULL, V = NULL,
                                                     verbose = FALSE, monitor = FALSE,
                                                     iter_max = 500, iter_min = 100,
                                                     init_mode = "random",
                                                     epsilon = 1e-3, additional_iter = 10,
                                                     conv_mode = 1, ninit = 1,
                                                     iter_init = 50, ncores = NCORE,
                                                     reorder_factor = FALSE, seed = NULL,
                                                     a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                     alpha1 = NULL, alpha2 = NULL,
                                                     beta1 = NULL, beta2 = NULL,
                                                     prob_S = prob_S, prior_S = prior_S,
                                                     prob_D = NULL, prior_D = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE) * tmp_res$sparse_param$prob_S
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE) * tmp_res$sparse_param$prob_S
        ## standard representation
        method <- "spcmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "spcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## selection
        filt_spcmf[filt_spcmf] <- tmp_res$sparse_param$prior_prob_S > 0.5
    }
    
    ### pCMF restim
    Xtmp <- X[,filt_spcmf]
    print("### pCMF (restim)")
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_gap_factor(Xtmp, K = ncomp, U = NULL, V = NULL,
                                              verbose = FALSE, monitor = FALSE,
                                              iter_max = 2000, iter_min = 1000,
                                              init_mode = "random",
                                              epsilon = 1e-3, additional_iter = 10,
                                              conv_mode = 1, ninit = 10,
                                              iter_init = 50, ncores = NCORE,
                                              reorder_factor = FALSE, seed = NULL,
                                              a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                              alpha1 = NULL, alpha2 = NULL,
                                              beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE)
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE)
        ## standard representation
        method <- "pcmf_restim"
        tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "pcmf_restim_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### pCMF
    print("### pCMF")
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_gap_factor(X, K = ncomp, U = NULL, V = NULL,
                                              verbose = FALSE, monitor = FALSE,
                                              iter_max = 2000, iter_min = 1000,
                                              init_mode = "random",
                                              epsilon = 1e-3, additional_iter = 10,
                                              conv_mode = 1, ninit = 10,
                                              iter_init = 50, ncores = NCORE,
                                              reorder_factor = FALSE, seed = NULL,
                                              a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                              alpha1 = NULL, alpha2 = NULL,
                                              beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE)
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE)
        ## standard representation
        method <- "pcmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "pcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### NMF
    print("### NMF")
    run_time <- system.time(
        tmp_res <- tryCatch(run_poisson_nmf(X, K = ncomp, U = NULL, V = NULL,
                                            verbose = FALSE, monitor = FALSE,
                                            iter_max = 500, iter_min = 100,
                                            epsilon = 1e-3, additional_iter = 10,
                                            conv_mode = 1, ninit = 10,
                                            iter_init = 50, ncores = NCORE,
                                            reorder_factor = FALSE, seed = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE)
        ## standard representation
        method <- "nmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label,
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### SPCA
    filt_spca <- rep(TRUE, ncol(Xlog))
    print("### SPCA")
    run_time <- system.time(
        tmp_res <- tryCatch(spca_call(Xlog, ncomp=ncomp, ncore=NCORE, verbose=FALSE),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev[ncomp]
        U <- tmp_res$U
        V <- tmp_res$V
        ## standard representation
        method <- "spca"
        tmp <- format_result(method, U, V, indiv_label, feature_label,
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## selection
        filt_spca <- apply(tmp_res$V, 1, function(x) return(any(x!=0)))
    }
    
    ### PCA (restim)
    Xlog_tmp <- Xlog[,filt_spca]
    print("### PCA (restim)")
    run_time <- system.time(
        tmp_res <- tryCatch(pca_call(Xlog_tmp, ncomp=ncomp),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = sum(tmp_res$exp_dev[1:ncomp])
        U <- tmp_res$U[,1:ncomp]
        V <- tmp_res$V[,1:ncomp]
        ## standard representation
        method <- "pca_restim"
        tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spca],
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### PCA
    print("### PCA")
    run_time <- system.time(
        tmp_res <- tryCatch(pca_call(Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = sum(tmp_res$exp_dev[1:ncomp])
        U <- tmp_res$U[,1:ncomp]
        V <- tmp_res$V[,1:ncomp]
        ## standard representation
        method <- "pca"
        tmp <- format_result(method, U, V, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### ZIFA
    print("### ZIFA")
    run_time <- system.time(
        tmp_res <- tryCatch(zifa_call(X=Xlog, ncomp=ncomp, ncore=NCORE),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        U <- tmp_res$U
        V <- tmp_res$V
        exp_dev <- tryCatch(exp_var(Xlog, U),
                            error=function(e) custom_error3(e))
        ## standard representation
        method <- "zifa"
        tmp <- format_result(method, U, V, indiv_label, feature_label[tmp_res$kept_genes], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### t-SNE
    if(ncomp==2) {
        print("### t-SNE")
        tmp_res <- tsne_analysis(Xlog, ncomp=2, verbose=FALSE, 
                                 ncore=NCORE, choose_best=TRUE,
                                 indiv_label=indiv_label,
                                 additional_factor=additional_factor)
        
        if(!is.null(tmp_res)) {
            tmp <- tmp_res$result
            if(!is.null(tmp)) {
                tmp <- tmp[,-ncol(tmp)]
                tmp_result <- rbind(tmp_result, tmp)
            }
            tmp <- tmp_res$data_to_plot_U
            if(!is.null(tmp)) {
                tmp <- tmp[,-ncol(tmp)]
                tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp)
            }
        }
    }
    
    res <- list(result=tmp_result, 
                data_to_plot_U=tmp_data_to_plot_U, 
                data_to_plot_V=tmp_data_to_plot_V)
    return(res)
}

## experiment 1bis: clustering and visualization (no NMF and no restim)
experiment1b <- function(input) {
    
    if(!is.data.frame(input)) {
        stop("pb with input")
    } else {
        if(nrow(input)>1) stop("pb with input")
    }
    
    print(input)
    
    rep <- input$rep
    ncomp <- input$ncomp
    drop_out <- input$drop_out
    group_separation <- input$group_separation
    prop_noise_feature <- input$prop_noise_feature
    nrep <- input$nrep
    n <- input$n
    p <- input$p
    K <- input$K
    ngroup <- input$ngroup
    
    if(!is.null(SEED)) {
        set.seed(SEED)
    }
    
    factorU <- NULL
    factorV <- NULL
    count_data <- NULL
    X <- NULL
    keep_genes <- NULL
    for(count in 1:50) {
        average_signal <- c(250,100,250)
        factorU <- generate_factor_matrix(n, K, ngroup=ngroup, 
                                          average_signal=average_signal,
                                          group_separation=group_separation,
                                          distribution="exponential",
                                          shuffle_feature=FALSE,
                                          tag="cell")
        
        tmp1 <- 20*(1-prop_noise_feature)
        tmp2 <- 20*prop_noise_feature
        prop_noise_feature <- round(rbeta(1, 
                                          shape1 = runif(1, min=tmp2-2, max=tmp2+2), 
                                          shape2 = runif(1, min=tmp1-2, max=tmp1+2)), 
                                    digits = 3)
        
        factorV <- generate_factor_matrix(p, K, ngroup=2, average_signal=80,
                                          group_separation=0.8,
                                          distribution="exponential",
                                          shuffle_feature=TRUE,
                                          prop_noise_feature=prop_noise_feature,
                                          tag="gene")
        U <- factorU$factor_matrix
        V <- factorV$factor_matrix
        
        tmp1 <- 20*(1-drop_out)
        tmp2 <- 20*drop_out
        prob1 <- round(rbeta(p, 
                             shape1 = runif(1, min=tmp1-2, max=tmp1+2), 
                             shape2 = runif(1, min=tmp2-2, max=tmp2+2)), 
                       digits = 3)
        
        count_data <- generate_count_matrix(n, p, K, U, V,
                                            ZI=TRUE, prob1=prob1)
        X <- count_data$X
        Xmax <- apply(X,2,max)
        keep_genes <- (Xmax>0) & (Xmax<=quantile(Xmax, probs=0.95)) #& (apply(X,2,function(x) sum(x!=0))>=0.05*nrow(X))
        
        if(sum(keep_genes) > 0.5 * p) break
    }
    
    X <- X[,keep_genes]
    Xlog <- log2m(X)
    indiv_label <- as.factor(factorU$feature_label)
    feature_label <- as.factor(factorV$feature_label[keep_genes])
    
    additional_factor <- cbind(input[1,],
                               data.frame(prop0 = 1-sum(X!=0)/prod(dim(X))),
                               data.frame(ngene = ncol(X)),
                               data.frame(prop_noise = prop_noise_feature))
    
    ## run
    print(paste0("dim X = ", paste0(dim(X), collapse=", ")))
    
    tmp_data_to_plot_U <- NULL
    tmp_data_to_plot_V <- NULL
    tmp_result <- NULL
    
    if(!is.null(SEED)) {
        set.seed(NULL)
    }
    
    ### sparse pCMF
    filt_spcmf <- rep(TRUE, ncol(X))
    print("### pCMF (sparse)")
    
    # prior_S <- 1-exp(-apply(X,2,sd)/mean(X[X!=0]))
    # filt_spcmf <- prior_S > 0.2
    Xtmp <- X[,filt_spcmf]
    # print(paste0("dim X = ", paste0(dim(Xtmp), collapse=", "), " after pCMF pre-filtering"))
    
    prior_S <- 1-exp(-apply(Xtmp,2,sd)/mean(Xtmp[Xtmp!=0]))
    prob_S <- matrix(rep(prior_S, ncomp), ncol=ncomp)
    
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_sparse_gap_factor(Xtmp, K = ncomp, U = NULL, V = NULL,
                                                     verbose = FALSE, monitor = FALSE,
                                                     iter_max = 500, iter_min = 100,
                                                     init_mode = "random",
                                                     epsilon = 1e-3, additional_iter = 10,
                                                     conv_mode = 1, ninit = 1,
                                                     iter_init = 50, ncores = NCORE,
                                                     reorder_factor = FALSE, seed = NULL,
                                                     a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                     alpha1 = NULL, alpha2 = NULL,
                                                     beta1 = NULL, beta2 = NULL,
                                                     prob_S = prob_S, prior_S = prior_S,
                                                     prob_D = NULL, prior_D = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE) * tmp_res$sparse_param$prob_S
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE) * tmp_res$sparse_param$prob_S
        ## standard representation
        method <- "spcmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "spcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## selection
        filt_spcmf[filt_spcmf] <- tmp_res$sparse_param$prior_prob_S > 0.5
    }
    
    # ### pCMF restim
    # Xtmp <- X[,filt_spcmf]
    # print("### pCMF (restim)")
    # run_time <- system.time(
    #     tmp_res <- tryCatch(run_zi_gap_factor(Xtmp, K = ncomp, U = NULL, V = NULL,
    #                                           verbose = FALSE, monitor = FALSE,
    #                                           iter_max = 2000, iter_min = 1000,
    #                                           init_mode = "random",
    #                                           epsilon = 1e-3, additional_iter = 10,
    #                                           conv_mode = 1, ninit = 10,
    #                                           iter_init = 50, ncores = NCORE,
    #                                           reorder_factor = FALSE, seed = NULL,
    #                                           a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
    #                                           alpha1 = NULL, alpha2 = NULL,
    #                                           beta1 = NULL, beta2 = NULL),
    #                         error=function(e) custom_error2(e)))[3]
    # 
    # 
    # if(!is.null(tmp_res)) {
    #     exp_dev = tmp_res$exp_dev
    #     U <- getU(tmp_res, log_representation=FALSE)
    #     V <- getV(tmp_res, log_representation=FALSE)
    #     Ulog <- getU(tmp_res, log_representation=TRUE)
    #     Vlog <- getV(tmp_res, log_representation=TRUE)
    #     ## standard representation
    #     method <- "pcmf_restim"
    #     tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spcmf], 
    #                          exp_dev, run_time,
    #                          additional_factor=additional_factor)
    #     
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    #     ## log representation
    #     method <- "pcmf_restim_log"
    #     tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
    #                          exp_dev, run_time,
    #                          additional_factor=additional_factor)
    #     
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    # }
    
    ### pCMF
    print("### pCMF")
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_gap_factor(X, K = ncomp, U = NULL, V = NULL,
                                              verbose = FALSE, monitor = FALSE,
                                              iter_max = 2000, iter_min = 1000,
                                              init_mode = "random",
                                              epsilon = 1e-3, additional_iter = 10,
                                              conv_mode = 1, ninit = 10,
                                              iter_init = 50, ncores = NCORE,
                                              reorder_factor = FALSE, seed = NULL,
                                              a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                              alpha1 = NULL, alpha2 = NULL,
                                              beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE)
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE)
        ## standard representation
        method <- "pcmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "pcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    # ### NMF
    # print("### NMF")
    # run_time <- system.time(
    #     tmp_res <- tryCatch(run_poisson_nmf(X, K = ncomp, U = NULL, V = NULL,
    #                                         verbose = FALSE, monitor = FALSE,
    #                                         iter_max = 500, iter_min = 100,
    #                                         epsilon = 1e-3, additional_iter = 10,
    #                                         conv_mode = 1, ninit = 10,
    #                                         iter_init = 50, ncores = NCORE,
    #                                         reorder_factor = FALSE, seed = NULL),
    #                         error=function(e) custom_error2(e)))[3]
    # 
    # if(!is.null(tmp_res)) {
    #     exp_dev = tmp_res$exp_dev
    #     U <- getU(tmp_res, log_representation=FALSE)
    #     V <- getV(tmp_res, log_representation=FALSE)
    #     ## standard representation
    #     method <- "nmf"
    #     tmp <- format_result(method, U, V, indiv_label, feature_label,
    #                          exp_dev, run_time,
    #                          additional_factor=additional_factor)
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    # }
    
    ### SPCA
    filt_spca <- rep(TRUE, ncol(Xlog))
    print("### SPCA")
    run_time <- system.time(
        tmp_res <- tryCatch(spca_call(Xlog, ncomp=ncomp, ncore=NCORE, verbose=FALSE),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev[ncomp]
        U <- tmp_res$U
        V <- tmp_res$V
        ## standard representation
        method <- "spca"
        tmp <- format_result(method, U, V, indiv_label, feature_label,
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## selection
        filt_spca <- apply(tmp_res$V, 1, function(x) return(any(x!=0)))
    }
    
    # ### PCA (restim)
    # Xlog_tmp <- Xlog[,filt_spca]
    # print("### PCA (restim)")
    # run_time <- system.time(
    #     tmp_res <- tryCatch(pca_call(Xlog_tmp, ncomp=ncomp),
    #                         error=function(e) custom_error2(e)))[3]
    # 
    # if(!is.null(tmp_res)) {
    #     exp_dev = sum(tmp_res$exp_dev[1:ncomp])
    #     U <- tmp_res$U[,1:ncomp]
    #     V <- tmp_res$V[,1:ncomp]
    #     ## standard representation
    #     method <- "pca_restim"
    #     tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spca],
    #                          exp_dev, run_time,
    #                          additional_factor=additional_factor)
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    # }
    
    ### PCA
    print("### PCA")
    run_time <- system.time(
        tmp_res <- tryCatch(pca_call(Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = sum(tmp_res$exp_dev[1:ncomp])
        U <- tmp_res$U[,1:ncomp]
        V <- tmp_res$V[,1:ncomp]
        ## standard representation
        method <- "pca"
        tmp <- format_result(method, U, V, indiv_label, feature_label, 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### ZIFA
    print("### ZIFA")
    run_time <- system.time(
        tmp_res <- tryCatch(zifa_call(X=Xlog, ncomp=ncomp, ncore=NCORE),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        U <- tmp_res$U
        V <- tmp_res$V
        exp_dev <- tryCatch(exp_var(Xlog, U),
                            error=function(e) custom_error3(e))
        ## standard representation
        method <- "zifa"
        tmp <- format_result(method, U, V, indiv_label, feature_label[tmp_res$kept_genes], 
                             exp_dev, run_time,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### t-SNE
    if(ncomp==2) {
        print("### t-SNE")
        tmp_res <- tsne_analysis(Xlog, ncomp=2, verbose=FALSE, 
                                 ncore=NCORE, choose_best=TRUE,
                                 indiv_label=indiv_label,
                                 additional_factor=additional_factor)
        
        if(!is.null(tmp_res)) {
            tmp <- tmp_res$result
            if(!is.null(tmp)) {
                tmp <- tmp[,-ncol(tmp)]
                tmp_result <- rbind(tmp_result, tmp)
            }
            tmp <- tmp_res$data_to_plot_U
            if(!is.null(tmp)) {
                tmp <- tmp[,-ncol(tmp)]
                tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp)
            }
        }
    }
    
    res <- list(result=tmp_result, 
                data_to_plot_U=tmp_data_to_plot_U, 
                data_to_plot_V=tmp_data_to_plot_V)
    return(res)
}

## experiment 2: gene selection
experiment2 <- function(input) {
    
    if(!is.data.frame(input)) {
        stop("pb with input")
    } else {
        if(nrow(input)>1) stop("pb with input")
    }
    
    print(input)
    
    rep <- input$rep
    ncomp <- input$ncomp
    drop_out <- input$drop_out
    group_separation <- input$group_separation
    prop_noise_feature <- input$prop_noise_feature
    nrep <- input$nrep
    n <- input$n
    p <- input$p
    K <- input$K
    ngroup <- input$ngroup
    
    if(!is.null(SEED)) {
        set.seed(SEED)
    }
    
    factorU <- NULL
    factorV <- NULL
    count_data <- NULL
    X <- NULL
    keep_genes <- NULL
    for(count in 1:50) {
        average_signal <- c(250,100,250)
        factorU <- generate_factor_matrix(n, K, ngroup=ngroup, 
                                          average_signal=average_signal,
                                          group_separation=group_separation,
                                          distribution="exponential",
                                          shuffle_feature=FALSE,
                                          tag="cell")
        
        tmp1 <- 20*(1-prop_noise_feature)
        tmp2 <- 20*prop_noise_feature
        prop_noise_feature <- round(rbeta(1, 
                                          shape1 = runif(1, min=tmp2-2, max=tmp2+2), 
                                          shape2 = runif(1, min=tmp1-2, max=tmp1+2)), 
                                    digits = 3)
        
        factorV <- generate_factor_matrix(p, K, ngroup=2, average_signal=80,
                                          group_separation=0.8,
                                          distribution="exponential",
                                          shuffle_feature=TRUE,
                                          prop_noise_feature=prop_noise_feature,
                                          tag="gene")
        U <- factorU$factor_matrix
        V <- factorV$factor_matrix
        
        tmp1 <- 20*(1-drop_out)
        tmp2 <- 20*drop_out
        prob1 <- round(rbeta(p, 
                             shape1 = runif(1, min=tmp1-2, max=tmp1+2), 
                             shape2 = runif(1, min=tmp2-2, max=tmp2+2)), 
                       digits = 3)
        
        count_data <- generate_count_matrix(n, p, K, U, V,
                                            ZI=TRUE, prob1=prob1)
        X <- count_data$X
        Xmax <- apply(X,2,max)
        keep_genes <- (Xmax>0) & (Xmax<=quantile(Xmax, probs=0.95)) #& (apply(X,2,function(x) sum(x!=0))>=0.05*nrow(X))
        
        if(sum(keep_genes) > 0.5 * p) break
    }
    
    X <- X[,keep_genes]
    Xlog <- log2m(X)
    indiv_label <- as.factor(factorU$feature_label)
    feature_label <- as.factor(factorV$feature_label[keep_genes])
    
    additional_factor <- cbind(input[1,],
                               data.frame(prop0 = 1-sum(X!=0)/prod(dim(X))),
                               data.frame(ngene = ncol(X)),
                               data.frame(prop_noise = prop_noise_feature))
    
    ## run
    print(paste0("dim X = ", paste0(dim(X), collapse=", ")))
    
    tmp_data_to_plot_U <- NULL
    tmp_data_to_plot_V <- NULL
    tmp_result <- NULL
    tmp_sel_spcmf <- NULL
    tmp_sel_spca <- NULL
    
    if(!is.null(SEED)) {
        set.seed(NULL)
    }
    
    ### sparse pCMF
    filt_spcmf <- rep(TRUE, ncol(X))
    print("### pCMF (sparse)")
    
    # prior_S <- 1-exp(-apply(X,2,sd)/mean(X[X!=0]))
    # filt_spcmf <- prior_S > 0.2
    Xtmp <- X[,filt_spcmf]
    # print(paste0("dim X = ", paste0(dim(Xtmp), collapse=", "), " after pCMF pre-filtering"))
    
    prior_S <- 1-exp(-apply(Xtmp,2,sd)/mean(Xtmp[Xtmp!=0]))
    prob_S <- matrix(rep(prior_S, ncomp), ncol=ncomp)
    
    run_time <- system.time(
        tmp_res <- tryCatch(run_zi_sparse_gap_factor(Xtmp, K = ncomp, U = NULL, V = NULL,
                                                     verbose = FALSE, monitor = FALSE,
                                                     iter_max = 500, iter_min = 100,
                                                     init_mode = "random",
                                                     epsilon = 1e-3, additional_iter = 10,
                                                     conv_mode = 1, ninit = 1,
                                                     iter_init = 50, ncores = NCORE,
                                                     reorder_factor = FALSE, seed = NULL,
                                                     a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                     alpha1 = NULL, alpha2 = NULL,
                                                     beta1 = NULL, beta2 = NULL,
                                                     prob_S = prob_S, prior_S = prior_S,
                                                     prob_D = NULL, prior_D = NULL),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev
        U <- getU(tmp_res, log_representation=FALSE)
        V <- getV(tmp_res, log_representation=FALSE) * tmp_res$sparse_param$prob_S
        Ulog <- getU(tmp_res, log_representation=TRUE)
        Vlog <- getV(tmp_res, log_representation=TRUE) * tmp_res$sparse_param$prob_S
        ## selection
        sel_gene <- tmp_res$sparse_param$prior_prob_S > 0
        accuracy <- 1 - (sum(as.integer(sel_gene) != as.integer(factorV$feature_label[keep_genes]!=0)) / ncol(X))
        tmp_additional_factor <- cbind(additional_factor,
                                       data.frame(accuracy=accuracy))
        tmp_sel_spcmf <- list(config=tmp_additional_factor,
                              truth=colnames(X)[factorV$feature_label[keep_genes]!=0],
                              selection=colnames(X)[sel_gene])
        ## standard representation
        method <- "spcmf"
        tmp <- format_result(method, U, V, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=tmp_additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "spcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time,
                             additional_factor=tmp_additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    ### SPCA
    print("### SPCA")
    run_time <- system.time(
        tmp_res <- tryCatch(spca_call(Xlog, ncomp=ncomp, ncore=NCORE, verbose=FALSE),
                            error=function(e) custom_error2(e)))[3]
    
    if(!is.null(tmp_res)) {
        exp_dev = tmp_res$exp_dev[ncomp]
        U <- tmp_res$U
        V <- tmp_res$V
        ## selection
        sel_gene <- apply(tmp_res$V, 1, function(x) return(any(x!=0)))
        accuracy <- 1 - (sum(as.integer(sel_gene) != as.integer(factorV$feature_label[keep_genes]!=0)) / ncol(X))
        tmp_additional_factor <- cbind(additional_factor,
                                       data.frame(accuracy=accuracy))
        tmp_sel_spca <- list(config=tmp_additional_factor,
                             truth=colnames(X)[factorV$feature_label[keep_genes]!=0],
                             selection=colnames(X)[sel_gene])
        ## standard representation
        method <- "spca"
        tmp <- format_result(method, U, V, indiv_label, feature_label,
                             exp_dev, run_time,
                             additional_factor=tmp_additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
    }
    
    res <- list(result=tmp_result, 
                data_to_plot_U=tmp_data_to_plot_U, 
                data_to_plot_V=tmp_data_to_plot_V,
                sel_spcmf=tmp_sel_spcmf,
                sel_spca=tmp_sel_spca)
    return(res)
}
