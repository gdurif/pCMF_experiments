## ----options, include=FALSE, echo=FALSE------------------------------------
library(BiocStyle)
knitr::opts_chunk$set(warning=FALSE, error=FALSE, message=FALSE)

## ----construct-------------------------------------------------------------
library(SingleCellExperiment)
counts <- matrix(rpois(100, lambda = 10), ncol=10, nrow=10)
sce <- SingleCellExperiment(assays = list(counts = counts))
sce

## ----coerce----------------------------------------------------------------
se <- SummarizedExperiment(list(counts=counts))
as(se, "SingleCellExperiment")

## ----fluidigm--------------------------------------------------------------
library(scRNAseq)
data(allen)
allen
sce <- as(allen, "SingleCellExperiment")
sce

## ----spikes----------------------------------------------------------------
isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
sce

## --------------------------------------------------------------------------
table(isSpike(sce, "ERCC"))
spikeNames(sce)

## ----spikes2---------------------------------------------------------------
isSpike(sce, "Adam") <- grepl("^Adam[0-9]", rownames(sce))
sce
table(isSpike(sce, "Adam"))
spikeNames(sce)

## --------------------------------------------------------------------------
table(isSpike(sce))

## --------------------------------------------------------------------------
temp <- sce
isSpike(temp, "Adam") <- NULL
spikeNames(temp)
temp <- clearSpikes(temp)
spikeNames(temp)

## ----sizeFactors-----------------------------------------------------------
sizeFactors(sce) <- colSums(assay(sce))
head(sizeFactors(sce))

## ----sizeFactors2----------------------------------------------------------
sizeFactors(sce, "ERCC") <- colSums(assay(sce)[isSpike(sce, "ERCC"),])
head(sizeFactors(sce, "ERCC"))
head(sizeFactors(sce)) # same as before

## ----metadata--------------------------------------------------------------
colData(sce)
rowData(sce)

## ----metadata2-------------------------------------------------------------
colData(sce, internal=TRUE)
rowData(sce, internal=TRUE)

## ----subset----------------------------------------------------------------
library(magrittr)
assay(sce) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)

sce_sub <- sce[names(vars[1:100]),]
sce_sub

## ----pca-------------------------------------------------------------------
library(Rtsne)
set.seed(5252)

pca_data <- prcomp(t(log1p(assay(sce_sub))))
tsne_data <- Rtsne(pca_data$x[,1:50], pca = FALSE)

reducedDims(sce_sub) <- SimpleList(PCA=pca_data$x, TSNE=tsne_data$Y)
sce_sub

## --------------------------------------------------------------------------
reducedDims(sce_sub)
reducedDimNames(sce_sub)
head(reducedDim(sce_sub, "PCA")[,1:2])
head(reducedDim(sce_sub, "TSNE")[,1:2])

## --------------------------------------------------------------------------
dim(reducedDim(sce_sub, "PCA"))
dim(reducedDim(sce_sub[,1:10], "PCA"))

## --------------------------------------------------------------------------
counts(sce) <- assay(sce, "tophat_counts")
sce
dim(counts(sce))

## --------------------------------------------------------------------------
sessionInfo()

