### function to run experiment

## experiment: clustering, visualization and selection
experiment <- function(input, dataset, X, indiv_label, feature_label) {
    
    if(!is.data.frame(input)) {
        stop("pb with input")
    } else {
        if(nrow(input)>1) stop("pb with input")
    }
    
    print(input)
    
    rep <- input$rep
    ncomp <- input$ncomp
    presel <- input$presel
    
    ## pre-selection
    if(presel) {
        prior_variance <- 1-exp(-apply(X,2,sd)/mean(X[X!=0]))
        keep_genes <- prior_variance > 0.2
        X <- X[,keep_genes]
        feature_label <- as.factor(feature_label[keep_genes])
        # sum(colnames(X) %in% dea_result$Gene) / length(dea_result$Gene)
        print(paste0("dim X = ", paste0(dim(X), collapse=", "), " after re-filtering"))
    }
    
    ## log count
    Xlog <- log2m(X)
    
    ## factors
    additional_factor <- cbind(data.frame(dataset = dataset),
                               data.frame(ngroup = length(unique(indiv_label))),
                               data.frame(prop0 = 1-sum(X!=0)/prod(dim(X))),
                               data.frame(ngene = ncol(X)))
    
    ## run
    print(paste0("dim X = ", paste0(dim(X), collapse=", ")))
    
    tmp_data_to_plot_U <- NULL
    tmp_data_to_plot_V <- NULL
    tmp_result <- NULL
    tmp_sel_spcmf <- NULL
    tmp_res_dimred <- list()
    
    if(!is.null(SEED)) {
        set.seed(NULL)
    }
    
    ### sparse pCMF
    filt_spcmf <- rep(TRUE, ncol(X))
    print("### pCMF (sparse)")

    Xtmp <- X[,filt_spcmf]
    
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
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "spcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## selection
        filt_spcmf[filt_spcmf] <- tmp_res$sparse_param$prior_prob_S > 0.5
        tmp_sel_spcmf <- list(ncomp=ncomp, genes=colnames(X)[filt_spcmf])
        ## res U and V
        rownames(U) <- rownames(Xtmp)
        rownames(V) <- colnames(Xtmp)
        rownames(Ulog) <- rownames(Xtmp)
        rownames(Vlog) <- colnames(Xtmp)
        tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
                                                           input=input, U=U, Ulog=Ulog,
                                                           V=V, Vlog=Vlog)
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
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "pcmf_restim_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label[filt_spcmf], 
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## res U and V
        rownames(U) <- rownames(Xtmp)
        rownames(V) <- colnames(Xtmp)
        rownames(Ulog) <- rownames(Xtmp)
        rownames(Vlog) <- colnames(Xtmp)
        tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
                                                           input=input, U=U, Ulog=Ulog,
                                                           V=V, Vlog=Vlog)
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
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## log representation
        method <- "pcmf_log"
        tmp <- format_result(method, Ulog, Vlog, indiv_label, feature_label, 
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## res U and V
        rownames(U) <- rownames(X)
        rownames(V) <- colnames(X)
        rownames(Ulog) <- rownames(X)
        rownames(Vlog) <- colnames(X)
        tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
                                                           input=input, U=U, Ulog=Ulog,
                                                           V=V, Vlog=Vlog)
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
    #                          exp_dev, run_time, replicate=rep,
    #                          additional_factor=additional_factor)
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     ## res U and V
    #     rownames(U) <- rownames(X)
    #     rownames(V) <- colnames(X)
    #     tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
    #                                                        input=input, U=U, V=V)
    # }

    # ### SPCA
    # filt_spca <- rep(TRUE, ncol(Xlog))
    # print("### SPCA")
    # run_time <- system.time(
    #     tmp_res <- tryCatch(spca_call(Xlog, ncomp=ncomp, ncore=NCORE, verbose=FALSE),
    #                         error=function(e) custom_error2(e)))[3]
    # 
    # if(!is.null(tmp_res)) {
    #     exp_dev = tmp_res$exp_dev[ncomp]
    #     U <- tmp_res$U
    #     V <- tmp_res$V
    #     ## standard representation
    #     method <- "spca"
    #     tmp <- format_result(method, U, V, indiv_label, feature_label,
    #                          exp_dev, run_time, replicate=rep,
    #                          additional_factor=additional_factor)
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     ## selection
    #     filt_spca <- apply(tmp_res$V, 1, function(x) return(any(x!=0)))
    #     ## res U and V
    #     rownames(U) <- rownames(Xlog)
    #     rownames(V) <- colnames(Xlog)
    #     tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
    #                                                        input=input, U=U, V=V)
    # }
    # 
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
    #                          exp_dev, run_time, replicate=rep,
    #                          additional_factor=additional_factor)
    #     tmp_result <- rbind(tmp_result, tmp$result)
    #     tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
    #     ## res U and V
    #     rownames(U) <- rownames(Xlog_tmp)
    #     rownames(V) <- colnames(Xlog_tmp)
    #     tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
    #                                                        input=input, U=U, V=V)
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
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## res U and V
        rownames(U) <- rownames(Xlog)
        rownames(V) <- colnames(Xlog)
        tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
                                                           input=input, U=U, V=V)
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
                             exp_dev, run_time, replicate=rep,
                             additional_factor=additional_factor)
        tmp_result <- rbind(tmp_result, tmp$result)
        tmp_data_to_plot_U <- rbind(tmp_data_to_plot_U, tmp$data_to_plot_U)
        tmp_data_to_plot_V <- rbind(tmp_data_to_plot_V, tmp$data_to_plot_V)
        ## res U and V
        rownames(U) <- rownames(Xlog)
        rownames(V) <- colnames(Xlog)[tmp_res$kept_genes]
        tmp_res_dimred[[length(tmp_res_dimred)+1]] <- list(method=method,  dataset=dataset,
                                                           input=input, U=U, V=V)
    }
    
    ### t-SNE
    if(ncomp==2) {
        print("### t-SNE")
        tmp_res <- tsne_analysis(Xlog, ncomp=2, verbose=FALSE, 
                                 ncore=NCORE, choose_best=TRUE,
                                 indiv_label=indiv_label,
                                 replicate=rep,
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
                data_to_plot_V=tmp_data_to_plot_V,
                sel_spcmf=tmp_sel_spcmf,
                res_dimred=tmp_res_dimred)
    return(res)
}

## run DE analysis with M3Drop
m3drop_dea <- function(X, indiv_label) {
    require(M3Drop)
    Normalized_data <- M3DropCleanData(t(X), 
                                       labels = indiv_label, 
                                       is.counts=TRUE, min_detected_genes=2000)
    DE_genes <- M3DropDifferentialExpression(Normalized_data$data, 
                                             mt_method="fdr", mt_threshold=0.01,
                                             suppress.plot=TRUE)
    return(DE_genes)
}
