### stability of selection

## src
source_dir(RENV) # util functions
library(pCMF)
source(file.path(PROJDIR, "experiment", 
                 "simulated_data", "src_experiment.R")) # code for experiments

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
    
    SEED <- parameters$rep[rep]
    
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
     file=file.path(DIR, "stability_selection.RData"))

################################################################
if(II) {
    ## loading result
    DIR <- file.path(RESDIR, "experiment", "simulated_data")
    tmp <- load(file=file.path(DIR, "stability_selection.RData"))
    
    ## printing result
    library(ggplot2)
    DIR <- file.path(PROJDIR, "figure", "experiment", "simulated_data")
    ensure_dir(DIR)
    
    ## plot TODO
}

# ## printing result
# library(ggplot2)
# 
# DIR <- file.path(RESDIR, "experiment", "simulated_data")
# load(file=file.path(DIR, "visu_effect_dropout.RData"))
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
# task <- paste0("visu_effect_dropout_cells")
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
# task <- paste0("visu_effect_dropout_genes")
# ggsave(filename=file.path(DIR, paste0(task, ".eps")), 
#        plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
