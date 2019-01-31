install.packages("devtools")
devtools::install_github("rhondabacher/SCnorm")
devtools::install_github("lme4/lme4")
devtools::install_github("BatzoglouLabSU/SIMLR")
devtools::install_github("Vivianstats/scImpute")

install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8", update = FALSE)
BiocManager::install("readat", version = "3.8", update = FALSE)
BiocManager::install("SCnorm", version = "3.8", update = FALSE)
BiocManager::install("sva", version = "3.8", update = FALSE)
BiocManager::install("SIMLR", version = "3.8", update = FALSE)
BiocManager::install("tximport", version = "3.8", update = FALSE)

install.packages("R2admb")
install.packages("glmmADMB", 
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")

devtools::install("../scRNAtools_tmp")