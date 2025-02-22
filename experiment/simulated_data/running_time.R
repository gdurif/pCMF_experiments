### computation time, effect of n, p and ncomp

## src
source_dir(RENV) # util functions
library(pCMF)
source(file.path(PROJDIR, "experiment", 
                 "simulated_data", "src_experiment.R")) # code for experiment

################################################################
## run experiment or directly load result (c.f. below)

## parameters
ncomp <- 2
ncompval <- c(2,10,20,30,40)
drop_out <- 0.7
group_separation <- 0.8
prop_noise_feature <- 0.6
nrep <- 20

n <- 1000
nval <- c(100, 500, 1000, 5000, 10000)
p <- 6000
pval <- c(600, 1200, 2500, 5000)
K <- 40
ngroup <- 3

## run experiment
parameters <- NULL

tmp_parameters <- expand.grid(ncomp=ncomp, drop_out=drop_out, 
                              group_separation=group_separation,
                              prop_noise_feature=prop_noise_feature,
                              rep=1:nrep, n=nval, p=p, K=K, ngroup=ngroup, 
                              KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)
parameters <- rbind(parameters, tmp_parameters)

tmp_parameters <- expand.grid(ncomp=ncomp, drop_out=drop_out, 
                              group_separation=group_separation,
                              prop_noise_feature=prop_noise_feature,
                              rep=1:nrep, n=n, p=pval, K=K, ngroup=ngroup, 
                              KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)
parameters <- rbind(parameters, tmp_parameters)

tmp_parameters <- expand.grid(ncomp=ncompval, drop_out=drop_out, 
                              group_separation=group_separation,
                              prop_noise_feature=prop_noise_feature,
                              rep=1:nrep, n=n, p=p, K=K, ngroup=ngroup, 
                              KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)
parameters <- rbind(parameters, tmp_parameters)

print(paste0("number of experiment = ", nrow(parameters)))

data_to_plot_U <- NULL
data_to_plot_V <- NULL
result <- NULL
for(rep in 1:nrow(parameters)) {
    
    print(paste0("##### replicate ", rep))
    
    tmp_res <- experiment1b(parameters[rep,])
    
    if(!is.null(tmp_res)) {
        result <- rbind(result, tmp_res$result)
        data_to_plot_U <- rbind(data_to_plot_U, tmp_res$data_to_plot_U)
        data_to_plot_V <- rbind(data_to_plot_V, tmp_res$data_to_plot_V)
    }
}

## averaging result
library(plyr)
library(dplyr)

library(doParallel)
registerDoParallel(cores=as.integer(NCORE))

property <- c("exp_dev", "run_time", "adjRI_U", "adjRI_V", "prop0", "ngene", "prop_noise")
feature <- colnames(result)[! colnames(result) %in% c(property, "rep", "replicate")]

average_result <- ddply(result, feature, 
                        function(x) data.frame(myMean(x[property], na.rm=TRUE),
                                               myVar(x[property], na.rm=TRUE)),
                        .parallel=TRUE)

## saving
DIR <- file.path(RESDIR, "experiment", "simulated_data")
ensure_dir(DIR)
save(list=c("data_to_plot_U", "data_to_plot_V", 
            "result", "average_result"), 
     file=file.path(DIR, "running_time.RData"))

################################################################
if(II) {
    ## loading result
    DIR <- file.path(RESDIR, "experiment", "simulated_data")
    tmp <- load(file=file.path(DIR, "running_time.RData"))
    
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
# load(file=file.path(DIR, "running_time.RData"))
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
# task <- paste0("visu_running_time")
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
# task <- paste0("visu_running_time")
# ggsave(filename=file.path(DIR, paste0(task, ".eps")), 
#        plot=g1, width=350, height=100, dpi=1200, scale=1, units="mm")
