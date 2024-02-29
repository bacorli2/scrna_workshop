

# Set current directory to the location of this file
# In RStudio, Menu >> Session >> Set Workspace Directory >> To Source File Location

# Base Packages
base_packages <- c("cli","tidyverse", "BiocManager", "cowplot", "here", "R.utils",
                   "HGNChelper","openxlsx", "devtools","VGAM", "metap")
base_installs <- install.packages(setdiff(base_packages, 
                                          rownames(installed.packages())), 
                                  dependencies = TRUE)
# Load all base packages
base_loads <- lapply(base_packages, library, character.only = TRUE)
# Bioconductor Packages
biocm_packages <-  c("BiocParallel", "Seurat", "celldex", "scAnnotatR", 
                     "scRNAseq", "multtest", "SingleR")
bioc_installs <- setdiff(biocm_packages, rownames(installed.packages()))
if (length(bioc_installs)) {BiocManager::install(bioc_installs) }
# Load all Bioconductor packages
# bioc_loads <- lapply(biocm_packages, library, character.only = TRUE)


# Print out what version of bioconductor manager you have, should by 3.18
BiocManager::version()
# Biocmanager will check if any packages will need to be updates
# If they do, it will print out a command to the console below that you can copy
# and run the console
BiocManager::valid()


# Installing from Github
#-------------------------------------------------------------------------------
# Below are packages installed directly from Github instead of CRAN
# If you request to download these packages too many times, you will get 
# "timed-out" and be forced to wait 30 minutes before you can download from 
# github. 
# If you get timed out while trying to debug the installs, you can get access 
# by setting up a personal github access token (you need a github account to do 
# so). Follow the directions here
# https://emilyriederer.github.io/projmgr/articles/github-pat.html


# Seurat supporting packages
#-------------------------------------------------------------------------------
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat-wrappers")

# Monocle 3 Developer Package and Dependencies
#-------------------------------------------------------------------------------
# Install github released version of monocle3 then supporting packages
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr', "qlcMatrix"))
# Then install github developer release second
devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")
# Presto is used for accelerating cluster marker identification
devtools::install_github('immunogenomics/presto')


# scPred Developer Package
#-------------------------------------------------------------------------------
devtools::install_github("immunogenomics/harmony")
devtools::install_github("powellgenomicslab/scPred")


# Set base directory of for workshop
here::i_am("README.md")

# Create data directory for workshop
dir.create(here::here("_temp_data"))
