
pkg_list <- c("devtools", 
              "rprojroot",
              "ggplot2",
              "MASS",
              "mclust",
              "RhpcBLASctl",
              "reticulate",
              "Rtsne",
              "plyr",
              "dplyr",
              "doParallel")

print(paste0("Installing the following R packages:", paste0(pkg_list, collapse=", ")))
install.packages(pkg_list)

warning("You should install ZIFA python environment with the './set_python_zifa.sh' script")

pkg_list <- c("BH",
              "fields",
              "Rcpp",
              "RcppEigen",
              "testthat")

install.packages(pkg_list)

devtools::install_git("https://gitlab.inria.fr/gdurif/pCMF", subdir="pkg", branch="prod")

BiocManager::install("SingleCellExperiment", version = "3.8", update = FALSE)
BiocManager::install("M3Drop", version = "3.8", update = FALSE)
BiocManager::install("M3DExampleData", version = "3.8", update = FALSE)
BiocManager::install("scRNAseq", version = "3.8", update = FALSE)

BiocManager::install("impute", version = "3.8", update = FALSE)
install.packages("PMA")

devtools::install_github("krlmlr/ulimit")
