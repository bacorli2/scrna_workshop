

# Base Packages
base_packages <- c("tidyverse", "BiocManager", "cowplot", "here", "R.utils", "HGNChelper","openxlsx", "devtools")
base_installs <- install.packages(setdiff(base_packages, rownames(installed.packages())), 
                                  dependencies = TRUE)
# Load all base packages
base_loads <- lapply(base_packages, library, character.only = TRUE)
# Bioconductor Packages
biocm_packages <-  c("BiocParallel", "Seurat", "SeuratWrappers","monocle","celldex")
bioc_installs <- setdiff(biocm_packages, rownames(installed.packages()))
if (length(bioc_installs)) {BiocManager::install(bioc_installs) }
# Load all Bioconductor packages
bioc_loads <- lapply(biocm_packages, library, character.only = TRUE)

# Custom pacakges
devtools::install_github('immunogenomics/presto')

# Set base directory of for workshop
here::i_am("README.md")

# Create data directory for workshop
dir.create(here::here("temp_data"))
