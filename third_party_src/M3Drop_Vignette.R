### R code from vignette source 'M3Drop_Vignette.Rnw'

###################################################
### code chunk number 1: M3Drop_Vignette.Rnw:20-22
###################################################
library(M3Drop)
library(M3DExampleData)


###################################################
### code chunk number 2: M3Drop_Vignette.Rnw:28-33
###################################################
Normalized_data <- M3DropCleanData(Mmus_example_list$data, 
		  labels = Mmus_example_list$labels, 
		  is.counts=TRUE, min_detected_genes=2000)
dim(Normalized_data$data)
length(Normalized_data$labels)


###################################################
### code chunk number 3: M3Drop_Vignette.Rnw:40-41
###################################################
fits <- M3DropDropoutModels(Normalized_data$data)


###################################################
### code chunk number 4: M3Drop_Vignette.Rnw:48-54
###################################################
# Sum absolute residuals
data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr,
	   DoubleExpo=fits$ExpoFit$SAr) 
# Sum squared residuals
data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr,
	   DoubleExpo=fits$ExpoFit$SSr)


###################################################
### code chunk number 5: M3Drop_Vignette.Rnw:63-65
###################################################
DE_genes <- M3DropDifferentialExpression(Normalized_data$data, 
			mt_method="fdr", mt_threshold=0.01)


###################################################
### code chunk number 6: M3Drop_Vignette.Rnw:78-80
###################################################
heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, Normalized_data$data, 
			cell_labels = Normalized_data$labels)


###################################################
### code chunk number 7: M3Drop_Vignette.Rnw:88-91
###################################################
cell_populations <- M3DropGetHeatmapCellClusters(heat_out, k=4)
library("ROCR") 
marker_genes <- M3DropGetMarkers(Normalized_data$data, cell_populations)


###################################################
### code chunk number 8: M3Drop_Vignette.Rnw:98-100
###################################################
head(marker_genes[marker_genes$Group==4,],20) 
marker_genes[rownames(marker_genes)=="Cdx2",] 


###################################################
### code chunk number 9: M3Drop_Vignette.Rnw:110-111
###################################################
HVG <- BrenneckeGetVariableGenes(Normalized_data$data)


###################################################
### code chunk number 10: M3Drop_Vignette.Rnw:122-124
###################################################
heat_out <- M3DropExpressionHeatmap(HVG, Normalized_data$data, 
		cell_labels = Normalized_data$labels)


