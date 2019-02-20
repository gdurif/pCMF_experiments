# Experiments related to pCMF development and evaluation

In this repository, you can find the codes and datasets that were used to write [1].

The pCMF package is available here: <https://gitlab.inria.fr/gdurif/pCMF>


[1] Durif, Ghislain, Laurent Modolo, Jeff E. Mold, Sophie Lambert-Lacroix, et Franck Picard. « Probabilistic Count Matrix Factorization for Single Cell Expression Data Analysis ». arXiv:1710.11028 [stat], 2019. http://arxiv.org/abs/1710.11028.


## COPYRIGHT

Copyright (c) [2019] [Ghislain DURIF] 

Applies to all files except the ones located in the 'third_party_src' and 'data/single_cell' 
sub-directories, which belong to the respective authors of the respective softwares 
or data.

Licensing: all files except the ones located in the 'third_party_src' and 'data/single_cell' 
sub-directories are distributed under the CECILL v2.1 license. You can check the license 
agreement in the './LICENSE.txt' file.


## CONTACT

If you encounter any issue, you can contact us at 'gd DOT dev AT libertymail DOT net'.


## HOW TO USE

This is some research code, that is not packaged for distribution, however it should be quite easy to use it to replicate the experiment from the paper, at least on Unix OS, tested on GNU/Linux (Ubuntu 16.04 and 18.04) but not on MacOS.

Download the repository:
```bash
git clone https://gitlab.inria.fr/gdurif/pCMF_experiments
cd pCMF_experiments
```

Set the environment (R and Python packages are installed locally in `pCMF_experiments/.R_libs` and `pCMF_experiments/.pyenv` respectively):
```bash
cd pCMF_experiments
source set_env.sh
Rscript set_pkg.R
bash set_python_zifa.sh
```

To benefit from mutli-threading (which increases the performance), you can set the `NCORE` variable in the `pCMF_experiments/.Renviron` file or set the environment variable `NCORE` before sourcing the bash script `pCMF_experiments/set_env.sh`:
```bash
export NCORE=8
source set_env.sh
```

You can then run scripts in the `pCMF_experiments/experiment` sub-folders, either in Rstudio (there is a specific project file `pCMF_experiments/pCMF_experiments.Rproj`), or with R command line tools (source the file `pCMF_experiments/set_env.sh` will set environment variables pointing to the project specific `.Renviron` and `.Rprofile` files.

Source codes defining experiment are located in the corresponding sub-folders. Some utilitary functions are defined in `pCMF_experiments/src/R`.

We mention that, for the moment, our codes use the low level functions `run_XXX` of the pCMF package, and not the wrapper `pCMF`. A new version based on this more simple wrapper will be soon published.
