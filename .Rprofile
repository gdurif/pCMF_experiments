## project
PROJDIR <- tryCatch(rprojroot::find_root(rprojroot::has_file(".Rprofile")), 
                    error=function(e) return(NULL))

if(is.null(PROJDIR)) {
    warning("You must run the `set_env.sh` and `set_pkg.R` scripts to set the environment and to install the required packages. /!\\ Otherwise, the scripts provided will not work.")
}

## additional variables
VAR <- c("DATADIR", "RESDIR", "NCORE")
tmp <- sapply(VAR, function(v) assign(x=v, value=Sys.getenv(v), envir=.GlobalEnv))
rm(tmp)

RENV <- file.path(PROJDIR, "src", "R")
PYENV <- file.path(PROJDIR, ".pyenv")
RESDIR <- file.path(PROJDIR, "result")
NCORE <- as.integer(NCORE)

SEED <-NULL

## preload library
library(parallel)

## interactive mode ?
II <- interactive()

################################################################################
## defines custom functions

# function for paraload use: write in chunkout when using paraload
chunkout_flush = function(FILE, ARGS, message) {
    if(!is.null(FILE)) {
        output = paste0(ARGS, collapse=" / ")
        write.table(data.frame(paste0("NOK: ", output)), file=FILE, 
                    append=FALSE, row.names=FALSE, 
                    col.names=FALSE, quote=FALSE)
    }
    stop(message)
}

# custom error catching for tryCatch
custom_error1 <- function(e) {print(e); return(e);}
custom_error2 <- function(e) {print(e); return(NULL);}
custom_error3 <- function(e) {return(NA);}

check_parallel_error = function(seq) {
    tmp <- sapply(1:length(seq), 
                  function(ind) return(any(is.character(seq[[ind]]))))
    check_error <- sapply(1:length(seq), function(ind) {
        if(tmp[ind]) {
            return(any(grepl(pattern="[Ee]rror", x=seq[[ind]])))
        } else {
            return(FALSE)
        }
    })
    if(any(check_error)) {
        msg=head(seq[check_error][[1]],1)
    } else {
        msg=NULL
    }
    return(list(error=any(check_error), nbError=sum(check_error), msg=msg))
}

# source an entire directory
source_dir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        if(trace) cat(nm,":")           
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
    }
}

# create directory
ensure_dir <- function(folder_name) {
    if(!dir.exists(folder_name)) dir.create(folder_name, recursive = TRUE)
}

# custom column wise mean and variance
myMean <- function(x, ...) {
    tmp <- colwise(mean)(x, ...)
    colnames(tmp) <- paste0(colnames(x), ".m")
    return(tmp)
}

myVar <- function(x, ...) {
    tmp <- colwise(var)(x, ...)
    colnames(tmp) <- paste0(colnames(x), ".sd2")
    return(tmp)
}

## load data
load_data <- function(filename) {
    if(grepl("RData", filename)) {
        tmp <- load(filename)
        return(sce)
    }
    else if(grepl("rds", filename)) {
        sce <- readRDS(filename)
        return(sce)
    } else {
        return(NULL)
    }
}

