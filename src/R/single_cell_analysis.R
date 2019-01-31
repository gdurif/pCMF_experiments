### run comparison on single_cell data set

run_analysis <- function(dataset) {
    
    print(paste0("################### ", dataset))
    
    tmp <- load_dataset(dataset, del_zero=TRUE, cluster=cl)
    X <- tmp$X
    Xlog <- log2m(X)
    
    indiv_labels <- tmp$indiv_labels
    id_indiv_labels <- as.numeric(factor(indiv_labels))
    nb_group <- tmp$nb_group
    
    n <- nrow(X)
    p <- ncol(X)
    print(paste0("dim = ", paste0(dim(X), collapse=", ")))
    
    perplex_value <- c(1:80)
    
    dataToPlot <- NULL
    res_pcmf <- NULL
    res_pcmf_restim <- NULL
    res_nmf <- NULL
    res_pcmf_nmf <- NULL
    res_pcmf_zifa <- NULL
    res_pcmf_10 <- NULL
    res_pcmf_restim_10 <- NULL
    res_pca <- NULL
    res_zifa <- NULL
    sapply(paste0("res_tsne", perplex_value), function(x) assign(x, NULL, envir=.GlobalEnv))
    sapply(paste0("res_tsne_pcmf", perplex_value), function(x) assign(x, NULL, envir=.GlobalEnv))
    
    print("### pCMF")
    require(pCMF)
    time_run <- system.time(
        res_tmp <- tryCatch(run_zi_sparse_gap_factor(X, K = ncomp, sel_bound = 0.5, U = NULL, V = NULL,
                                                     verbose = TRUE, monitor = FALSE,
                                                     iter_max = 2000, iter_min = 1000,
                                                     init_mode = "random",
                                                     epsilon = 1e-8, additional_iter = 10,
                                                     conv_mode = 1, ninit = 10,
                                                     iter_init = 50, ncores = ncores,
                                                     reorder_factor = FALSE, seed = NULL,
                                                     a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                     alpha1 = NULL, alpha2 = NULL,
                                                     beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))
    
    
    if(!is.null(res_tmp)) {
        
        exp_dev = res_tmp$exp_dev
        U <- getU(res_tmp, log_representation=FALSE)
        V <- getV(res_tmp, log_representation=FALSE)
        Ulog <- getU(res_tmp, log_representation=TRUE)
        Vlog <- getV(res_tmp, log_representation=TRUE)
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        res_cluster2 <- tryCatch(kmeans_call(data=Ulog, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                 error=function(e) custom_error2(e))
        
        res_pcmf <- list(time = time_run,
                         res = res_tmp, 
                         exp_dev = exp_dev,
                         U = U,
                         Ulog = Ulog,
                         V = V,
                         Vlog = Vlog,
                         res_cluster = res_cluster,
                         res_cluster2 = res_cluster2)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pcmf", nrow(X)),
                                       dataset=rep(dataset, nrow(X))),
                            data.frame(comp1=Ulog[,1],
                                       comp2=Ulog[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pcmf_logRep", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    print("### pCMF restim")
    require(pCMF)
    if(!is.null(res_pcmf)) {
        time_run <- system.time(
            res_tmp <- tryCatch(run_zi_gap_factor(X[,res_pcmf$res$sparse_param$prior_prob_S!=0], 
                                                  K = ncomp, U = NULL, V = NULL,
                                                  verbose = TRUE, monitor = FALSE,
                                                  iter_max = 2000, iter_min = 1000,
                                                  init_mode = "random",
                                                  epsilon = 1e-8, additional_iter = 10,
                                                  conv_mode = 1, ninit = 10,
                                                  iter_init = 50, ncores = ncores,
                                                  reorder_factor = FALSE, seed = NULL,
                                                  a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                  alpha1 = NULL, alpha2 = NULL,
                                                  beta1 = NULL, beta2 = NULL),
                                error=function(e) custom_error2(e)))
        
        
        if(!is.null(res_tmp)) {
            
            exp_dev = res_tmp$exp_dev
            U <- getU(res_tmp, log_representation=FALSE)
            V <- getV(res_tmp, log_representation=FALSE)
            Ulog <- getU(res_tmp, log_representation=TRUE)
            Vlog <- getV(res_tmp, log_representation=TRUE)
            
            res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                    error=function(e) custom_error2(e))
            res_cluster2 <- tryCatch(kmeans_call(data=Ulog, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                     error=function(e) custom_error2(e))
            
            res_pcmf_restim <- list(time = time_run,
                                    res = res_tmp, 
                                    exp_dev = exp_dev,
                                    U = U,
                                    Ulog = Ulog,
                                    V = V,
                                    Vlog = Vlog,
                                    res_cluster = res_cluster,
                                    res_cluster2 = res_cluster2)
            
            dataToPlot <- rbind(dataToPlot,
                                data.frame(comp1=U[,1],
                                           comp2=U[,2],
                                           cell_type=indiv_labels,
                                           perplexity=rep(NA, nrow(X)),
                                           method=rep("pcmf_restim", nrow(X)),
                                           dataset=rep(dataset, nrow(X))),
                                data.frame(comp1=Ulog[,1],
                                           comp2=Ulog[,2],
                                           cell_type=indiv_labels,
                                           perplexity=rep(NA, nrow(X)),
                                           method=rep("pcmf_restim_logRep", nrow(X)),
                                           dataset=rep(dataset, nrow(X))))
        }
    }
    
    # print("### NMF")
    # require(pCMF)
    # time_run <- system.time(
    #     res_tmp <- tryCatch(run_poisson_nmf(X, K=ncomp, U = NULL, V = NULL,
    #                                         verbose = TRUE, monitor = FALSE,
    #                                         iter_max = 2000, iter_min = 1700,
    #                                         epsilon = 1e-8, additional_iter = 10,
    #                                         conv_mode = 1, ninit = 10,
    #                                         iter_init = 50, ncores = ncores,
    #                                         reorder_factor = FALSE,
    #                                         seed = NULL),
    #                         error=function(e) custom_error2(e)))
    # 
    # if(!is.null(res_tmp)) {
    #     
    #     exp_dev = res_tmp$exp_dev
    #     U <- getU(res_tmp)
    #     V <- getV(res_tmp)
    #     
    #     res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
    #                             error=function(e) custom_error2(e))
    #     
    #     res_nmf <- list(time = time_run,
    #                     res = res_tmp, 
    #                     exp_dev = exp_dev,
    #                     U = U,
    #                     V = V,
    #                     res_cluster = res_cluster)
    #     
    #     dataToPlot <- rbind(dataToPlot,
    #                         data.frame(comp1=U[,1],
    #                                    comp2=U[,2],
    #                                    cell_type=indiv_labels,
    #                                    perplexity=rep(NA, nrow(X)),
    #                                    method=rep("nmf", nrow(X)),
    #                                    dataset=rep(dataset, nrow(X))))
    # }
    
    # print("### pCMF after NMF")
    # if(!is.null(res_nmf)) {
    #     require(pCMF)
    #     time_run <- system.time(
    #         res_tmp <- tryCatch(run_zi_gap_factor(X, K=ncomp, U = res_nmf$U, V = res_nmf$V,
    #                                               verbose = TRUE, monitor = FALSE,
    #                                               iter_max = 2000, iter_min = 1000,
    #                                               init_mode = "random",
    #                                               epsilon = 1e-8, additional_iter = 10,
    #                                               conv_mode = 1, ninit = 1,
    #                                               iter_init = 50, ncores = ncores,
    #                                               reorder_factor = FALSE, seed = NULL,
    #                                               a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
    #                                               alpha1 = NULL, alpha2 = NULL,
    #                                               beta1 = NULL, beta2 = NULL),
    #                             error=function(e) custom_error2(e)))
    #     
    #     if(!is.null(res_tmp)) {
    #         
    #         exp_dev = res_tmp$exp_dev
    #         U <- getU(res_tmp, log_representation=FALSE)
    #         V <- getV(res_tmp, log_representation=FALSE)
    #         Ulog <- getU(res_tmp, log_representation=TRUE)
    #         Vlog <- getV(res_tmp, log_representation=TRUE)
    #         
    #         res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
    #                                 error=function(e) custom_error2(e))
    #         res_cluster2 <- tryCatch(kmeans_call(data=Ulog, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
    #                                  error=function(e) custom_error2(e))
    #         
    #         res_pcmf_nmf <- list(time = time_run,
    #                              res = res_tmp, 
    #                              exp_dev = exp_dev,
    #                              U = U,
    #                              Ulog = Ulog,
    #                              V = V,
    #                              Vlog = Vlog,
    #                              res_cluster = res_cluster,
    #                              res_cluster2 = res_cluster2)
    #         
    #         dataToPlot <- rbind(dataToPlot,
    #                             data.frame(comp1=U[,1],
    #                                        comp2=U[,2],
    #                                        cell_type=indiv_labels,
    #                                        perplexity=rep(NA, nrow(X)),
    #                                        method=rep("pcmf_nmf", nrow(X)),
    #                                        dataset=rep(dataset, nrow(X))),
    #                             data.frame(comp1=Ulog[,1],
    #                                        comp2=Ulog[,2],
    #                                        cell_type=indiv_labels,
    #                                        perplexity=rep(NA, nrow(X)),
    #                                        method=rep("pcmf_nmf_logRep", nrow(X)),
    #                                        dataset=rep(dataset, nrow(X))))
    #     }
    # }
    
    print("### PCA")
    time_run <- system.time(
        res_tmp <- tryCatch(pca_call(Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))
    
    if(!is.null(res_tmp)) {
        
        exp_dev = sum(res_tmp$exp_dev[1:ncomp])
        U <- res_tmp$U[,1:ncomp]
        V <- res_tmp$V[,1:ncomp]
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        
        res_pca <- list(time = time_run,
                        res = res_tmp, 
                        exp_dev = exp_dev,
                        U = U,
                        V = V,
                        res_cluster = res_cluster)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pca", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    print("### ZIFA")
    time_run <- system.time(
        res_tmp <- tryCatch(zifa_call(X=Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))
    
    if(!is.null(res_tmp)) {
        
        U <- res_tmp$U
        V <- res_tmp$V
        exp_dev <- tryCatch(exp_var(Xlog, U),
                            error=function(e) custom_error2(e))
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        
        res_zifa <- list(time = time_run,
                         res = res_tmp,
                         exp_dev = exp_dev,
                         U = U,
                         V = V, 
                         res_cluster = res_cluster)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("zifa", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    # print("### pCMF after ZIFA")
    # if(!is.null(res_zifa)) {
    #     U0 <- res_zifa$U
    #     V0 <- res_zifa$V
    #     
    #     ## translate into the positive quadrant
    #     dir_vector_U <- apply(U0, 2, min)
    #     dir_vector_V <- apply(V0, 2, min)
    #     
    #     U0 <- U0 - matrix(rep(dir_vector_U, nrow(U0)), nrow=nrow(U0), byrow=TRUE) + 0.1
    #     V0 <- V0 - matrix(rep(dir_vector_V, nrow(V0)), nrow=nrow(V0), byrow=TRUE) + 0.1
    #     
    #     require(pCMF)
    #     time_run <- system.time(
    #         res_tmp <- tryCatch(run_zi_gap_factor(X, K=ncomp, U = U0, V = V0,
    #                                               verbose = TRUE, monitor = FALSE,
    #                                               iter_max = 2000, iter_min = 1000,
    #                                               init_mode = "random",
    #                                               epsilon = 1e-8, additional_iter = 10,
    #                                               conv_mode = 1, ninit = 1,
    #                                               iter_init = 50, ncores = ncores,
    #                                               reorder_factor = FALSE, seed = NULL,
    #                                               a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
    #                                               alpha1 = NULL, alpha2 = NULL,
    #                                               beta1 = NULL, beta2 = NULL),
    #                             error=function(e) custom_error2(e)))
    #     
    #     if(!is.null(res_tmp)) {
    #         
    #         exp_dev = res_tmp$exp_dev
    #         U <- getU(res_tmp, log_representation=FALSE)
    #         V <- getV(res_tmp, log_representation=FALSE)
    #         Ulog <- getU(res_tmp, log_representation=TRUE)
    #         Vlog <- getV(res_tmp, log_representation=TRUE)
    #         
    #         res_cluster <- kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10)
    #         res_cluster2 <- kmeans_call(data=Ulog, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10)
    #         
    #         res_pcmf_zifa <- list(time = time_run,
    #                               res = res_tmp, 
    #                               exp_dev = exp_dev,
    #                               U = U,
    #                               Ulog = Ulog,
    #                               V = V,
    #                               Vlog = Vlog,
    #                               res_cluster = res_cluster,
    #                               res_cluster2 = res_cluster2)
    #         
    #         dataToPlot <- rbind(dataToPlot,
    #                             data.frame(comp1=U[,1],
    #                                        comp2=U[,2],
    #                                        cell_type=indiv_labels,
    #                                        perplexity=rep(NA, nrow(X)),
    #                                        method=rep("pcmf_zifa", nrow(X)),
    #                                        dataset=rep(dataset, nrow(X))),
    #                             data.frame(comp1=Ulog[,1],
    #                                        comp2=Ulog[,2],
    #                                        cell_type=indiv_labels,
    #                                        perplexity=rep(NA, nrow(X)),
    #                                        method=rep("pcmf_zifa_logRep", nrow(X)),
    #                                        dataset=rep(dataset, nrow(X))))
    #     }
    # }
    
    print("### t-sne")
    for(perplex in perplex_value) {
        
        time_run <- system.time(
            res_tmp <- tryCatch(tsne_call(Xlog, ncomp=2, verbose=TRUE, perplexity=perplex),
                                error=function(e) custom_error2(e)))
        
        if(!is.null(res_tmp)) {
            
            U <- res_tmp$U
            exp_dev <- tryCatch(exp_var(Xlog, U),
                                error=function(e) custom_error2(e))
            
            res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                    error=function(e) custom_error2(e))
            
            
            assign(paste0("res_tsne", perplex), 
                   list(time = time_run,
                        res = res_tmp,
                        exp_dev = exp_dev,
                        U = U,
                        res_cluster = res_cluster), 
                   envir=.GlobalEnv)
            
            dataToPlot <- rbind(dataToPlot,
                                data.frame(comp1=U[,1],
                                           comp2=U[,2],
                                           cell_type=indiv_labels,
                                           perplexity=rep(perplex, nrow(X)),
                                           method=rep("tsne", nrow(X)),
                                           dataset=rep(dataset, nrow(X))))
        }
    }
    
    print("### pCMF ncomp = 10")
    require(pCMF)
    time_run <- system.time(
        res_tmp <- tryCatch(run_zi_sparse_gap_factor(X, K = 10, sel_bound = 0.5, U = NULL, V = NULL,
                                                     verbose = TRUE, monitor = FALSE,
                                                     iter_max = 2000, iter_min = 1000,
                                                     init_mode = "random",
                                                     epsilon = 1e-8, additional_iter = 10,
                                                     conv_mode = 1, ninit = 10,
                                                     iter_init = 50, ncores = ncores,
                                                     reorder_factor = FALSE, seed = NULL,
                                                     a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                     alpha1 = NULL, alpha2 = NULL,
                                                     beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))
    
    if(!is.null(res_tmp)) {
        
        exp_dev = res_tmp$exp_dev
        U <- getU(res_tmp, log_representation=FALSE)
        V <- getV(res_tmp, log_representation=FALSE)
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        
        res_pcmf_10 <- list(time = time_run,
                            res = res_tmp, 
                            exp_dev = exp_dev,
                            U = U,
                            V = V,
                            res_cluster = res_cluster)
    }
    
    print("### pCMF restim ncomp = 10")
    require(pCMF)
    if(!is.null(res_pcmf_10)) {
        time_run <- system.time(
            res_tmp <- tryCatch(run_zi_gap_factor(X[,res_pcmf_10$res$sparse_param$prior_prob_S!=0], K = 10, U = NULL, V = NULL,
                                                  verbose = TRUE, monitor = FALSE,
                                                  iter_max = 2000, iter_min = 1000,
                                                  init_mode = "random",
                                                  epsilon = 1e-8, additional_iter = 10,
                                                  conv_mode = 1, ninit = 10,
                                                  iter_init = 50, ncores = ncores,
                                                  reorder_factor = FALSE, seed = NULL,
                                                  a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                                  alpha1 = NULL, alpha2 = NULL,
                                                  beta1 = NULL, beta2 = NULL),
                                error=function(e) custom_error2(e)))
        
        
        if(!is.null(res_tmp)) {
            
            exp_dev = res_tmp$exp_dev
            U <- getU(res_tmp, log_representation=FALSE)
            V <- getV(res_tmp, log_representation=FALSE)
            
            res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                    error=function(e) custom_error2(e))
            
            res_pcmf_restim_10 <- list(time = time_run,
                                       res = res_tmp, 
                                       exp_dev = exp_dev,
                                       U = U,
                                       V = V,
                                       res_cluster = res_cluster)
        }
    }
    
    print("### t-sne on pCMF output")
    if(!is.null(res_pcmf_restim_10)) {
        for(perplex in perplex_value) {
            
            time_run <- system.time(
                res_tmp <- tryCatch(tsne_call(res_pcmf_restim_10$U, ncomp=2, verbose=TRUE, perplexity=perplex, pca=FALSE),
                                    error=function(e) custom_error2(e)))
            
            if(!is.null(res_tmp)) {
                
                U <- res_tmp$U
                exp_dev <- tryCatch(exp_var(Xlog, U),
                                    error=function(e) custom_error2(e))
                
                res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                        error=function(e) custom_error2(e))
                
                
                assign(paste0("res_tsne_pcmf", perplex), 
                       list(time = time_run,
                            res = res_tmp,
                            exp_dev = exp_dev,
                            U = U,
                            res_cluster = res_cluster), 
                       envir=.GlobalEnv)
                
                dataToPlot <- rbind(dataToPlot,
                                    data.frame(comp1=U[,1],
                                               comp2=U[,2],
                                               cell_type=indiv_labels,
                                               perplexity=rep(perplex, nrow(X)),
                                               method=rep("tsne_pcmf", nrow(X)),
                                               dataset=rep(dataset, nrow(X))))
            }
        }
    }
    
    #### save results
    save(list=c("res_pcmf", "res_pca", "res_zifa", "res_nmf", "res_pcmf_nmf", "res_pcmf_zifa", "res_pcmf_10",
                "res_pcmf_restim", "res_pcmf_restim_10",
                paste0("res_tsne", perplex_value), 
                paste0("res_tsne_pcmf", perplex_value), 
                "dataToPlot"), 
         file=file.path(LOCRESDIR, paste0(dataset, "_pCMF_sparse_restim_results_ncomp_", ncomp, ".RData")))
    
}

## version 2

run_analysis2 <- function(dataset, ncomp) {
    
    print(paste0("################### ", dataset))
    print(paste0("## ncomp = ", ncomp))
    
    tmp <- load_dataset(dataset, del_zero=TRUE, cluster=cl)
    X <- tmp$X
    Xlog <- tmp$Xlog
    
    indiv_labels <- tmp$cell_label
    id_indiv_labels <- tmp$cell_id
    nb_group <- tmp$nb_group
    
    n <- nrow(X)
    p <- ncol(X)
    print(paste0("dim = ", paste0(dim(X), collapse=", ")))
    
    perplex_value <- c(1,seq(5,100,5))
    
    dataToPlot <- NULL
    res_pcmf <- NULL
    res_pcmf_nmf <- NULL
    res_pcmf_zifa <- NULL
    res_pca <- NULL
    res_zifa <- NULL
    sapply(paste0("res_tsne", perplex_value), function(x) assign(x, NULL, envir=.GlobalEnv))
    sapply(paste0("res_tsne_pcmf", perplex_value), function(x) assign(x, NULL, envir=.GlobalEnv))
    sapply(paste0("res_tsne_zifa", perplex_value), function(x) assign(x, NULL, envir=.GlobalEnv))
    
    print("### pCMF")
    require(pCMF)
    time_run <- system.time(
        res_tmp <- tryCatch(run_zi_gap_factor(X, K = ncomp, U = NULL, V = NULL,
                                              verbose = TRUE, monitor = FALSE,
                                              iter_max = 2000, iter_min = 1000,
                                              init_mode = "random",
                                              epsilon = 1e-8, additional_iter = 10,
                                              conv_mode = 1, ninit = 10,
                                              iter_init = 50, ncores = ncores,
                                              reorder_factor = FALSE, seed = NULL,
                                              a1 = NULL, a2 = NULL, b1 = NULL, b2 = NULL,
                                              alpha1 = NULL, alpha2 = NULL,
                                              beta1 = NULL, beta2 = NULL),
                            error=function(e) custom_error2(e)))
    
    
    if(!is.null(res_tmp)) {
        
        exp_dev = res_tmp$exp_dev
        U <- getU(res_tmp, log_representation=FALSE)
        V <- getV(res_tmp, log_representation=FALSE)
        Ulog <- getU(res_tmp, log_representation=TRUE)
        Vlog <- getV(res_tmp, log_representation=TRUE)
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        res_cluster2 <- tryCatch(kmeans_call(data=Ulog, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                 error=function(e) custom_error2(e))
        
        res_pcmf <- list(time = time_run,
                         res = res_tmp, 
                         exp_dev = exp_dev,
                         U = U,
                         Ulog = Ulog,
                         V = V,
                         Vlog = Vlog,
                         res_cluster = res_cluster,
                         res_cluster2 = res_cluster2)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pcmf", nrow(X)),
                                       dataset=rep(dataset, nrow(X))),
                            data.frame(comp1=Ulog[,1],
                                       comp2=Ulog[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pcmf_logRep", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    print("### PCA")
    time_run <- system.time(
        res_tmp <- tryCatch(pca_call(Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))
    
    if(!is.null(res_tmp)) {
        
        exp_dev = sum(res_tmp$exp_dev[1:ncomp])
        U <- res_tmp$U[,1:ncomp]
        V <- res_tmp$V[,1:ncomp]
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        
        res_pca <- list(time = time_run,
                        res = res_tmp, 
                        exp_dev = exp_dev,
                        U = U,
                        V = V,
                        res_cluster = res_cluster)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("pca", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    print("### ZIFA")
    time_run <- system.time(
        res_tmp <- tryCatch(zifa_call(X=Xlog, ncomp=ncomp),
                            error=function(e) custom_error2(e)))
    
    if(!is.null(res_tmp)) {
        
        U <- res_tmp$U
        V <- res_tmp$V
        exp_dev <- tryCatch(exp_var(Xlog, U),
                            error=function(e) custom_error2(e))
        
        res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                error=function(e) custom_error2(e))
        
        res_zifa <- list(time = time_run,
                         res = res_tmp,
                         exp_dev = exp_dev,
                         U = U,
                         V = V, 
                         res_cluster = res_cluster)
        
        dataToPlot <- rbind(dataToPlot,
                            data.frame(comp1=U[,1],
                                       comp2=U[,2],
                                       cell_type=indiv_labels,
                                       perplexity=rep(NA, nrow(X)),
                                       method=rep("zifa", nrow(X)),
                                       dataset=rep(dataset, nrow(X))))
    }
    
    print("### t-sne")
    for(perplex in perplex_value) {
        
        time_run <- system.time(
            res_tmp <- tryCatch(tsne_call(Xlog, ncomp=2, verbose=TRUE, perplexity=perplex, initial_dims=ncomp),
                                error=function(e) custom_error2(e)))
        
        if(!is.null(res_tmp)) {
            
            U <- res_tmp$U
            exp_dev <- tryCatch(exp_var(Xlog, U),
                                error=function(e) custom_error2(e))
            
            res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                    error=function(e) custom_error2(e))
            
            assign(paste0("res_tsne", perplex), 
                   list(time = time_run,
                        res = res_tmp,
                        exp_dev = exp_dev,
                        U = U,
                        res_cluster = res_cluster), 
                   envir=.GlobalEnv)
            
            dataToPlot <- rbind(dataToPlot,
                                data.frame(comp1=U[,1],
                                           comp2=U[,2],
                                           cell_type=indiv_labels,
                                           perplexity=rep(perplex, nrow(X)),
                                           method=rep("tsne", nrow(X)),
                                           dataset=rep(dataset, nrow(X))))
        }
    }
    
    print("### t-sne on pCMF output")
    if(!is.null(res_pcmf)) {
        for(perplex in perplex_value) {
            
            time_run <- system.time(
                res_tmp <- tryCatch(tsne_call(res_pcmf$U, ncomp=2, verbose=TRUE, perplexity=perplex, pca=FALSE),
                                    error=function(e) custom_error2(e)))
            
            if(!is.null(res_tmp)) {
                
                U <- res_tmp$U
                exp_dev <- tryCatch(exp_var(Xlog, U),
                                    error=function(e) custom_error2(e))
                
                res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                        error=function(e) custom_error2(e))
                
                
                assign(paste0("res_tsne_pcmf", perplex), 
                       list(time = time_run,
                            res = res_tmp,
                            exp_dev = exp_dev,
                            U = U,
                            res_cluster = res_cluster), 
                       envir=.GlobalEnv)
                
                dataToPlot <- rbind(dataToPlot,
                                    data.frame(comp1=U[,1],
                                               comp2=U[,2],
                                               cell_type=indiv_labels,
                                               perplexity=rep(perplex, nrow(X)),
                                               method=rep("tsne_pcmf", nrow(X)),
                                               dataset=rep(dataset, nrow(X))))
            }
        }
    }
    
    print("### t-sne on ZIFA output")
    if(!is.null(res_zifa)) {
        for(perplex in perplex_value) {
            
            time_run <- system.time(
                res_tmp <- tryCatch(tsne_call(res_zifa$U, ncomp=2, verbose=TRUE, perplexity=perplex, pca=FALSE),
                                    error=function(e) custom_error2(e)))
            
            if(!is.null(res_tmp)) {
                
                U <- res_tmp$U
                exp_dev <- tryCatch(exp_var(Xlog, U),
                                    error=function(e) custom_error2(e))
                
                res_cluster <- tryCatch(kmeans_call(data=U, data_labels=indiv_labels, k=nb_group, iterMax=200, nstart=10),
                                        error=function(e) custom_error2(e))
                
                
                assign(paste0("res_tsne_zifa", perplex), 
                       list(time = time_run,
                            res = res_tmp,
                            exp_dev = exp_dev,
                            U = U,
                            res_cluster = res_cluster), 
                       envir=.GlobalEnv)
                
                dataToPlot <- rbind(dataToPlot,
                                    data.frame(comp1=U[,1],
                                               comp2=U[,2],
                                               cell_type=indiv_labels,
                                               perplexity=rep(perplex, nrow(X)),
                                               method=rep("tsne_zifa", nrow(X)),
                                               dataset=rep(dataset, nrow(X))))
            }
        }
    }
    
    #### save results
    save(list=c("res_pcmf", "res_pca", "res_zifa", 
                paste0("res_tsne", perplex_value), 
                paste0("res_tsne_pcmf", perplex_value), 
                paste0("res_tsne_zifa", perplex_value), 
                "dataToPlot"), 
         file=file.path(LOCRESDIR, paste0(dataset, "_results_", ncomp, "comp.RData")))
    
}
