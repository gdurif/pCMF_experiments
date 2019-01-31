
## clustering in latent sub-space after dimension reduction
# Input:
#    U: matrix of latent representation of individuals (cells)
#    V: matrix of latent representation of features (genes)
#    method: string, name of the method
#    indiv_label: vector of cell labels
#    feature_label: vector of gene clusters
clustering <- function(U, V, method, indiv_label, feature_label) {
    
    cluster_U <- NULL
    cluster_V <- NULL
    
    ## k-means clustering on cell sub-space
    if(!is.null(U) & !is.null(indiv_label)) {
        cluster_U <- tryCatch(kmeans_call(data=U, 
                                          data_labels=indiv_label, 
                                          k=length(unique(indiv_label)), 
                                          iterMax=200, 
                                          nstart=10),
                              error=function(e) custom_error2(e))
    }
    
    ## k-means clustering on gene sub-space
    if(!is.null(V) & !is.null(feature_label)) {
        cluster_V <- tryCatch(kmeans_call(data=V, 
                                          data_labels=feature_label, 
                                          k=length(unique(feature_label)), 
                                          iterMax=200, 
                                          nstart=10),
                              error=function(e) custom_error2(e))
    }
    
    res <- list(cluster_U=cluster_U,
                cluster_V=cluster_V)
    
    return(res)
}


## format latent factors for graphical representation
# Input:
#    U: matrix of latent representation of individuals (cells)
#    V: matrix of latent representation of features (genes)
#    method: string, name of the method
#    indiv_label: vector of cell labels
#    feature_label: vector of gene clusters
#    replicate: a tag to identify the replicate if necessary
dimred_result <- function(U, V, method, indiv_label, feature_label, replicate=NULL) {
    
    cell_plot <- NULL
    gene_plot <- NULL
    
    ## format results
    if(!is.null(U) & !is.null(indiv_label)) {
        cell_plot <- data.frame(comp1=U[,1],
                                comp2=U[,2],
                                cell_cluster=indiv_label,
                                method=rep(method, nrow(U)))
        
        if(!is.null(replicate)) {
            cell_plot$replicate <- rep(replicate, nrow(cell_plot))
        }
    }
    
    if(!is.null(V) & !is.null(feature_label)) {
        gene_plot <- data.frame(comp1=V[,1],
                                comp2=V[,2],
                                gene_cluster=feature_label,
                                method=rep(method, nrow(V)))
        
        if(!is.null(replicate)) {
            gene_plot$replicate <- rep(replicate, nrow(gene_plot))
        }
    }
    
    res <- list(cell_plot=cell_plot,
                gene_plot=gene_plot)
    
    return(res)
}

## format quantitative results after dimension reduction
# Input:
#    method: string, name of the method
#    exp_dev: numeric, explained deviance
#    run_time: numeric, computation time
#    clustering_res: results from function `clustering`
#    replicate: a tag to identify the replicate if necessary
quant_result <- function(method, exp_dev, run_time, clustering_res, replicate=NULL) {
    
    ## format results
    tmp_adjRI_U <- ifelse(!is.null(clustering_res$cluster_U), 
                          clustering_res$cluster_U$adjRI, NA)
    tmp_adjRI_V <- ifelse(!is.null(clustering_res$cluster_V), 
                          clustering_res$cluster_V$adjRI, NA)
    
    res <- data.frame(method=method, exp_dev=exp_dev, run_time=run_time,
                      adjRI_U=tmp_adjRI_U,
                      adjRI_V=tmp_adjRI_V,
                      row.names=NULL)
    
    if(!is.null(replicate)) {
        res$replicate <- rep(replicate, nrow(res))
    }
    
    return(res)
}


## formating result
# run clustering + quant_result + dimred_result
# Input:
#    method: string, name of the method
#    U: matrix of latent representation of individuals (cells)
#    V: matrix of latent representation of features (genes)
#    indiv_label: vector of cell labels
#    feature_label: vector of gene clusters
#    exp_dev: numeric, explained deviance
#    run_time: numeric, computation time
#    replicate: a tag to identify the replicate if necessary
#    additional_factor: data.frame with additional factor characterizing the run
format_result <- function(method, U, V, indiv_label, feature_label, 
                          exp_dev, run_time, replicate=NULL, 
                          additional_factor=NULL) {
    
    if(!is.null(additional_factor) & !is.data.frame(additional_factor)) {
        stop("'additional_factor' is not a data.frame")
    }
    
    result <- NULL
    data_to_plot_U <- NULL
    data_to_plot_V <- NULL
    
    ## clustering
    tmp_clust <- clustering(U, V, method, indiv_label, feature_label)
    
    ## quantitative results
    tmp <- quant_result(method, exp_dev, run_time, tmp_clust, replicate)
    tmp_additional_factor <- additional_factor[rep(1, nrow(tmp)),]
    rownames(tmp_additional_factor) <- NULL
    result <- cbind(tmp, tmp_additional_factor)
    
    ## dimension reduction
    tmp <- dimred_result(U, V, method, indiv_label, feature_label, replicate)
    if(!is.null(tmp$cell_plot)) {
        tmp_additional_factor <- additional_factor[rep(1, nrow(tmp$cell_plot)),]
        rownames(tmp_additional_factor) <- paste0(1:nrow(tmp$cell_plot))
        data_to_plot_U <- cbind(tmp$cell_plot, tmp_additional_factor)
    }
    if(!is.null(tmp$gene_plot)) {
        tmp_additional_factor <- additional_factor[rep(1, nrow(tmp$gene_plot)),]
        rownames(tmp_additional_factor) <- NULL
        data_to_plot_V <- cbind(tmp$gene_plot, tmp_additional_factor)
    }
    
    ## return
    res <- list(result=result, data_to_plot_U=data_to_plot_U, data_to_plot_V=data_to_plot_V)
    return(res)
}
