### t-SNE impact of perplexity

## src
source_dir(RENV)
library(pCMF)

## parameters
ncomp <- 2
drop_out <- c(0.3,0.5,0.7,0.9)
nrep <- 20

n <- 100
p <- 800
K <- 40
ngroup <- 3

## experiment

experiment <- function(input) {
    
    if(!is.data.frame(input)) if(nrow(input)>1) stop("pb with input")
    
    rep <- input$rep
    ncomp <- input$ncomp
    drop_out <- input$drop_out
    nrep <- input$nrep
    n <- input$n
    p <- input$p
    K <- input$K
    ngroup <- input$ngroup
    
    average_signal <- c(250,100,250)
    factorU <- generate_factor_matrix(n, K, ngroup=ngroup, 
                                      average_signal=average_signal,
                                      group_separation=0.8,
                                      distribution="exponential",
                                      shuffle_feature=FALSE,
                                      tag="cell")
    
    factorV <- generate_factor_matrix(p, K, ngroup=2, average_signal=80,
                                      group_separation=0.8,
                                      distribution="exponential",
                                      shuffle_feature=TRUE,
                                      prop_noise_feature=0.6,
                                      tag="gene")
    U <- factorU$factor_matrix
    V <- factorV$factor_matrix
    
    tmp <- 10*(1-drop_out)
    prob1 <- round(rbeta(p, 
                         shape1 = runif(tmp, min=tmp-1, max=tmp+2), 
                         shape2 = runif(1, min=12, max=22)), 
                   digits = 3)
    
    count_data <- generate_count_matrix(n, p, K, U, V,
                                        ZI=TRUE, prob1=prob1)
    X <- count_data$X
    Xmax <- apply(X,2,max)
    keep_genes <- Xmax>0 & (Xmax<quantile(Xmax, probs=0.95))
    X <- X[,keep_genes]
    Xlog <- log2m(X)
    indiv_label <- as.factor(factorU$feature_label)
    feature_label <- as.factor(factorV$feature_label[keep_genes])
    
    ## run
    print(paste0("dim = ", paste0(dim(X), collapse=", ")))
    
    data_to_plot_U <- NULL
    result <- NULL
    
    tmp_res <- tsne_analysis(Xlog, ncomp, verbose=FALSE, 
                             ncore=NCORE,
                             choose_best=FALSE,
                             indiv_label=indiv_label,
                             replicate=rep)
    
    if(!is.null(tmp_res)) {
        
        tmp <- tmp_res$result
        tmp$p0 <- rep(sum(X!=0)/prod(dim(X)), nrow(tmp))
        result <- rbind(result, tmp)
        
        tmp <- tmp_res$data_to_plot_U
        tmp$p0 <- rep(sum(X!=0)/prod(dim(X)), nrow(tmp))
        data_to_plot_U <- rbind(data_to_plot_U, tmp)
    }
    
    res <- list(result=result, data_to_plot_U=data_to_plot_U)
    return(res)
}

## run experiment
parameters <- expand.grid(ncomp=ncomp, drop_out=drop_out,
                          rep=1:nrep, n=n, p=p, K=K, ngroup=ngroup, 
                          KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)

data_to_plot_U <- NULL
result <- NULL
for(rep in 1:nrow(parameters)) {
    
    tmp_res <- experiment(parameters[rep,])
    
    if(!is.null(tmp_res)) {
        result <- rbind(result, tmp_res=result)
        data_to_plot_U <- rbind(data_to_plot_U, tmp_res$data_to_plot_U)
    }
}

## printing result
g1 <- ggplot(data_to_plot_U, aes(x=comp1, y=comp2)) 
g1 <- g1 + geom_point(aes(col=cell_cluster))
g1 <- g1 + facet_wrap(~method, scales="free")
g1 <- g1 + xlab("") + ylab("") + ggtitle("")
g1 <- g1 + theme_bw()
# g1

task <- paste0("visu_data_Baron_mouse_cell")
ggsave(filename=file.path(PROJDIR, "figure", paste0(task, ".eps")), 
       plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
