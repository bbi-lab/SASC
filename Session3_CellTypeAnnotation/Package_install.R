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
library(SeuratData)
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

install.packages("tidyverse")
library(tidyverse)

install.packages("dplyr")
library(dplyr)

install.packages("pheatmap")
library(pheatmap)

################################################################################
######################       Download data for tutorial   ######################
################################################################################

# We recommend create a directory for the following data so that it's easily accessible during the workshop
system("wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", TRUE)  # This will need to be untarred 
system("wget https://krishna.gs.washington.edu/content/members/hGastruloid_website/public/Profiling/RA_hGas_120h.RDS", TRUE)
system("wget https://krishna.gs.washington.edu/content/members/weiy/2024_sasc_workshop/moca_e9.5_subset.rds", TRUE)
system("wget https://krishna.gs.washington.edu/content/members/weiy/2024_sasc_workshop/NNLS_helper.R", TRUE)
