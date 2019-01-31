### Llorens-Bobadilla 2015 data set
### checking data
ulimit::memory_limit(16000)

## src
source_dir(RENV)
# library(M3Drop, quietly=TRUE)
library(pCMF, quietly=TRUE)
library(ggplot2, quietly=TRUE)

#### data
dataset <- "LlorensBobadilla2015"

tmp <- load(file=file.path(DATADIR, "single_cell", 
                           "LlorensBobadilla2015", 
                           paste0(dataset, "_genes_counts.RData")))

cell_id <- read.table(file.path(DATADIR, "single_cell", 
                                "LlorensBobadilla2015", 
                                paste0(dataset, "_cells_and_types.csv")),
                      header=TRUE, sep=",", stringsAsFactors=FALSE)
cell_id <- cell_id[,-1]

print(paste0("Original dimension of ", dataset, ": ", paste0(dim(neurons_genes), collapse=", ")))

# ctx, str, qNSC1, qNSC2, aNSC1, aNSC2, TAP
groups_to_keep <- c(paste0("qNSC", 1:2), paste0("aNSC", 1:2), "CTX", "TAP")

## X
X <- round(as.matrix(neurons_genes))

## cells to keep
cells_to_keep <- cell_id$cell[cell_id$type %in% groups_to_keep]
keep_cells <- rownames(sce) %in% cells_to_keep
X <- X[keep_cells,]
indiv_labels <- merge(data.frame(cell=rownames(X)), cell_id, by="cell")$type

print("groups")
print(table(indiv_labels))

feature_label <- rep(1, ncol(X))
rm(neurons_genes)

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