### selection (accuracy, stability)

## src
source_dir(RENV) # util functions
library(pCMF)
source(file.path(PROJDIR, "experiment", 
                 "simulated_data", "src_experiment.R")) # code for experiment

################################################################
## run experiment or directly load result (c.f. below)

## parameters
ncomp <- c(2,4,6,8,10)
drop_out <- 0.7
group_separation <- 0.8
prop_noise_feature <- c(0.2,0.4,0.6,0.8)
nrep <- 50

n <- 100
p <- 800
K <- 40
ngroup <- 3

## run experiment
parameters <- expand.grid(ncomp=ncomp, drop_out=drop_out, 
                          group_separation=group_separation,
                          prop_noise_feature=prop_noise_feature,
                          rep=1:nrep, n=n, p=p, K=K, ngroup=ngroup, 
                          KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)

print(paste0("number of experiment = ", nrow(parameters)))

data_to_plot_U <- NULL
data_to_plot_V <- NULL
result <- NULL
sel_spcmf <- list()
sel_spca <- list()
for(rep in 1:nrow(parameters)) {
    
    print(paste0("##### replicate ", rep))
    
    tmp_res <- experiment2(parameters[rep,])
    
    if(!is.null(tmp_res)) {
        result <- rbind(result, tmp_res$result)
        data_to_plot_U <- rbind(data_to_plot_U, tmp_res$data_to_plot_U)
        data_to_plot_V <- rbind(data_to_plot_V, tmp_res$data_to_plot_V)
        sel_spcmf[[length(sel_spcmf)+1]] <- tmp_res$sel_spcmf
        sel_spca[[length(sel_spca)+1]] <- tmp_res$sel_spca
    }
}

## averaging result
library(plyr)
library(dplyr)

library(doParallel)
registerDoParallel(cores=as.integer(NCORE))

property <- c("exp_dev", "run_time", "adjRI_U", "adjRI_V", "prop0", "ngene", "prop_noise", "accuracy")
feature <- colnames(result)[! colnames(result) %in% c(property, "rep", "replicate")]

average_result <- ddply(result, feature, 
                        function(x) data.frame(myMean(x[property], na.rm=TRUE),
                                               myVar(x[property], na.rm=TRUE)),
                        .parallel=TRUE)

## saving
DIR <- file.path(RESDIR, "experiment", "simulated_data")
ensure_dir(DIR)
save(list=c("data_to_plot_U", "data_to_plot_V", 
            "result", "average_result", 
            "sel_spcmf", "sel_spca"), 
     file=file.path(DIR, "gene_selection.RData"))

################################################################
if(II) {
    ## loading result
    DIR <- file.path(RESDIR, "experiment", "simulated_data")
    tmp <- load(file=file.path(DIR, "gene_selection.RData"))
    
    ## printing result
    library(ggplot2)
    DIR <- file.path(PROJDIR, "figure", "experiment", "simulated_data")
    ensure_dir(DIR)
    
    ## accuracy result (supp)
    method_to_plot <- c("spcmf", "spca")
    data_to_plot <- subset(average_result, method %in% method_to_plot & prop_noise_feature == 0.6)
    
    data_to_plot$method <- factor(data_to_plot$method,
                                  levels=c("spcmf", "spca"),
                                  labels=c("sparse pCMF", "sparse PCA"))
    
    g1 <- ggplot(data_to_plot, 
                 aes(x=ncomp, y=accuracy.m, 
                     ymin=accuracy.m - sqrt(accuracy.sd2) * 2/sqrt(50), 
                     ymax=accuracy.m + sqrt(accuracy.sd2) * 2/sqrt(50)))
    g1 <- g1 + geom_point(aes(col=method, shape=method))
    g1 <- g1 + geom_errorbar(width=0.3, size=0.5, aes(col=method))
    g1 <- g1 + geom_line(aes(col=method))
    g1 <- g1 + xlab("K") + ylab("accuracy") + ggtitle("")
    g1 <- g1 + ylim(c(-0.05,1))
    g1 <- g1 + theme_bw() + theme(legend.position="none")
    g1
    
    task <- paste0("supp_gene_selection_accuracy_ncomp")
    ggsave(filename=file.path(DIR, paste0(task, ".eps")),
           plot=g1, width=60, height=60, dpi=1200, scale=1, units="mm")
    
    ## accuracy result (supp)
    method_to_plot <- c("spcmf", "spca")
    data_to_plot <- subset(average_result, method %in% method_to_plot & ncomp == 2)
    
    data_to_plot$method <- factor(data_to_plot$method,
                                  levels=c("spcmf", "spca"),
                                  labels=c("sparse pCMF", "sparse PCA"))
    
    g1 <- ggplot(data_to_plot, 
                 aes(x=prop_noise_feature, y=accuracy.m, 
                     ymin=accuracy.m - sqrt(accuracy.sd2) * 2/sqrt(50), 
                     ymax=accuracy.m + sqrt(accuracy.sd2) * 2/sqrt(50)))
    g1 <- g1 + geom_point(aes(col=method, shape=method))
    g1 <- g1 + geom_errorbar(width=0.07, size=0.5, aes(col=method))
    g1 <- g1 + geom_line(aes(col=method))
    g1 <- g1 + xlab("proportion noisy genes") + ylab("accuracy") + ggtitle("")
    g1 <- g1 + ylim(c(-0.05,1))
    g1 <- g1 + theme_bw() #+ theme(legend.position="none")
    g1
    
    task <- paste0("supp_gene_selection_accuracy_noise")
    ggsave(filename=file.path(DIR, paste0(task, ".eps")),
           plot=g1, width=95, height=60, dpi=1200, scale=1, units="mm")
    
    
    
    ## cell plot (paper)
    method_to_plot <- c("pcmf_restim_log", "pca_restim", "zifa", "tsne")
    data_to_plot <- subset(data_to_plot_U, rep==1 & group_separation == 0.8 & drop_out == 0.9 & method %in% method_to_plot)
    
    data_to_plot$method <- factor(data_to_plot$method,
                                  levels=c("nmf", "pcmf_restim_log", "pca_restim", "zifa", "tsne"),
                                  labels=c("NMF", "pCMF", "PCA", "ZIFA", "t-SNE"))
    g1 <- ggplot(data_to_plot, 
                 aes(x=comp1, y=comp2))
    g1 <- g1 + geom_point(aes(col=cell_cluster))
    g1 <- g1 + facet_wrap(~method, scales="free", nrow=1)
    g1 <- g1 + xlab("") + ylab("") + ggtitle("")
    g1 <- g1 + theme_bw() #+ theme(legend.position = 'bottom')
    g1
    
    task <- paste0("paper_visu_effect_dropout_cells")
    ggsave(filename=file.path(DIR, paste0(task, ".eps")),
           plot=g1, width=310, height=80, dpi=1200, scale=1, units="mm")
    
    
    ## gene plot (paper)
    method_to_plot <- c("spcmf_log", "pca_restim", "zifa", "tsne")
    data_to_plot <- subset(data_to_plot_V, replicate==1 & group_separation == 0.8 & drop_out == 0.9 & method %in% method_to_plot)
    
    data_to_plot$method <- factor(data_to_plot$method,
                                  levels=c("nmf", "pcmf_restim_log", "pca_restim", "zifa", "tsne"),
                                  labels=c("NMF", "pCMF", "PCA", "ZIFA", "t-SNE"))
    g1 <- ggplot(data_to_plot, 
                 aes(x=comp1, y=comp2))
    g1 <- g1 + geom_point(aes(col=gene_cluster))
    g1 <- g1 + facet_wrap(~method, scales="free", nrow=1)
    g1 <- g1 + xlab("") + ylab("") + ggtitle("")
    g1 <- g1 + theme_bw() + theme(legend.position = 'bottom')
    g1
    
    task <- paste0("paper_visu_effect_dropout_genes")
    ggsave(filename=file.path(DIR, paste0(task, ".eps")),
           plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
}

# ## printing result
# library(ggplot2)
# 
# DIR <- file.path(RESDIR, "experiment", "simulated_data")
# load(file=file.path(DIR, "visu_effect_noisy_genes.RData"))
# 
# DIR <- file.path(PROJDIR, "figure", "experiment", "simulated_data")
# ensure_dir(DIR)
# 
# ## cell plot
# g1 <- ggplot(subset(data_to_plot_U, replicate==1), aes(x=comp1, y=comp2)) 
# g1 <- g1 + geom_point(aes(col=cell_cluster))
# g1 <- g1 + facet_wrap(drop_out~method, scales="free")
# g1 <- g1 + xlab("") + ylab("") + ggtitle("")
# g1 <- g1 + theme_bw()
# # g1
# 
# task <- paste0("visu_effect_noisy_genes_cells")
# ggsave(filename=file.path(DIR, paste0(task, ".eps")), 
#        plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
# 
# 
# g2 <- ggplot(subset(data_to_plot_V, rep==2), aes(x=comp1, y=comp2)) 
# g2 <- g2 + geom_point(aes(col=gene_cluster))
# g2 <- g2 + facet_wrap(~method, scales="free")
# g2 <- g2 + xlab("") + ylab("") + ggtitle("")
# g2 <- g2 + theme_bw()
# # g2
# 
# task <- paste0("visu_effect_noisy_genes_genes")
# ggsave(filename=file.path(DIR, paste0(task, ".eps")), 
#        plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
