################################################################################
########### Follow Seurat website guidelines for install (below) ###############
################################################################################
# Enter commands in R (or R studio)
install.packages('Seurat')
library(Seurat)

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

# Check that you can load Seurat
library(Seurat)
library(dplyr)
library(tidyverse)
library(Azimuth)

################################################################################
########### Install some additional packages for tutorial ######################
################################################################################
# Enter commands in R (or R studio)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(org.Hs.eg.db)

BiocManager::install("biomaRt")
library(biomaRt)

install.packages("nnls")
library(nnls)



