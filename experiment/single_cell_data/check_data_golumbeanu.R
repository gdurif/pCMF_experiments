### Baron data set
### checking data
ulimit::memory_limit(6000)

## src
source_dir(RENV)
# library(SingleCellExperiment, quietly=TRUE)
# library(M3Drop, quietly=TRUE)
library(pCMF, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(scatter)

#### data
filename <- file.path(DATADIR, "single_cell", "Golumbeanu2018", "GSE111727_sc_lat_model_normalized_counts_genes.txt")
sce <- t(read.table(filename, stringsAsFactors=FALSE))

# print(colnames(rowData(sce)))
# print(colnames(colData(sce)))
print("groups")
print(table(colData(sce)$Truth))

## X
X <- exprs(sce)
indiv_label <- rep(1, ncol(X))
feature_label <- rep(1, ncol(X))
rm(sce)

print(paste0("mean X = ", mean(X)))
print(paste0("max X = ", max(X)))

## check
print(paste0("dim X = ", paste0(dim(X), collapse=", ")))
print(paste0("prop 0 = ", sum(X!=0)/prod(dim(X))))

## pre-filtering
Xmax <- apply(X,2,max)

print(paste0("nb genes max>0 = ", sum(Xmax>0)))
print(paste0("nb genes max<quantile 95 = ", sum(Xmax<=quantile(Xmax, probs=0.95))))
print(paste0("nb genes expressed in at least 5% of cells = ", sum(apply(X,2,function(x) sum(x!=0))>=0.05*nrow(X))))

keep_genes <- (Xmax>0) & (Xmax<=quantile(Xmax, probs=0.95)) & (apply(X,2,function(x) sum(x!=0))>=0.05*nrow(X))

X <- X[,keep_genes]
print(paste0("dim X = ", paste0(dim(X), collapse=", "), " after pre-filtering"))
print(paste0("prop 0 = ", sum(X!=0)/prod(dim(X)), " after pre-filtering"))

## filtering
prior_S <- 1-exp(-apply(X,2,sd)/mean(X[X!=0]))
# hist(prior_S)

threshold_value <- seq(0.1,0.5,0.1)
tmp <- sapply(threshold_value, function(threshold) {
    filt <- prior_S > threshold
    print(paste0("nb pre-selected genes = ", sum(filt), " for threshold ", threshold))
})

filt <- prior_S > 0.2
X <- X[,filt]
print(paste0("dim X = ", paste0(dim(X), collapse=", "), " after filtering"))
print(paste0("prop 0 = ", sum(X!=0)/prod(dim(X)), " after filtering"))