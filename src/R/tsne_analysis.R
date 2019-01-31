## running tsne_call for different value of perplexity and chosing the best one

tsne_analysis <- function(X, ncomp, verbose=FALSE, 
                          perplex_value=c(1,seq(5,100,5)), ncore=NULL,
                          choose_best=FALSE,
                          indiv_label=NULL,
                          replicate=NULL, 
                          additional_factor=NULL, ...) {
    
    method <- "tsne"
    
    result <- NULL
    data_to_plot_U <- NULL
    
    for(perplex in perplex_value) {
        
        print(paste0("perplexity = ", perplex))
        
        if(is.null(additional_factor)) {
            tmp_additional_factor <- data.frame(perplexity = perplex)
        } else {
            tmp_additional_factor <- cbind(additional_factor, 
                                       data.frame(perplexity = perplex))
        }
        
        ## run tsne
        run_time <- system.time(
            tmp_res <- tryCatch(tsne_call(X, ncomp=ncomp, verbose=verbose, 
                                          perplexity=perplex, ncore=ncore, 
                                          initial_dims=ncomp),
                                error=function(e) custom_error2(e)))[3]
        
        if(!is.null(tmp_res)) {
            
            U <- tmp_res$U
            exp_dev <- tryCatch(exp_var(Xlog, U),
                                error=function(e) custom_error3(e))
            
            tmp <- format_result(method, U, V=NULL, indiv_label, feature_label=NULL, 
                                 exp_dev, run_time, replicate=replicate,
                                 additional_factor=tmp_additional_factor)
            
            result <- rbind(result, tmp$result)
            data_to_plot_U <- rbind(data_to_plot_U, tmp$data_to_plot_U)
        } else {
            break
        }
    }
    
    if(choose_best) {
        if(!is.null(result) & !is.null(data_to_plot_U)) {
            best_perplex <- which.max(result$adjRI_U)
            data_to_plot_U <- subset(data_to_plot_U, 
                                     perplexity == result$perplexity[best_perplex])
            result <- result[best_perplex,]
        }
    }
    
    res <- list(result=result, data_to_plot_U=data_to_plot_U)
    return(res)
    
}
