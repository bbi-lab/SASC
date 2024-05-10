################################################################################
########### SASC User Group Meeting #2 - 'Batch Effects' #######################
################################################################################
#Tutorial brought to you by the BAT-Lab Single Cell Genomics Group
#Authors: Anh Vo and Mary O'Neill
#2024-05-10

#Setup
setwd("~/Documents/SASC/") #CHANGE-ME!

#load some packages
library(monocle3)
library(ggplot2)
library(tidyverse)
library(viridis)
library(randomcoloR)
library(kBET)
library(cowplot)
library(harmony)
library(uwot)
library(batchelor)

################################################################################
########### READ DATA ##########################################################
################################################################################
#Read in the cell data set containing a random sub-sampling of 35K barcodes from
sasc <- readRDS("BBI_heart_hs_mix_36601humangenes_35000barcodes.RDS")

#a recent experiment
class(sasc) #cell_data_set
sasc #dim: 36601 35000 

#access the count data
dim(counts(sasc))
counts(sasc)[1:10,1:10] #sparse matrix, "." = 0; single-cell data is zero inflated! 

?detect_genes
sasc <- detect_genes(sasc)
expressed <- data.frame(rowData(sasc)) %>% arrange(desc(num_cells_expressed))
counts(sasc)[head(expressed$id, n=10), 1:10]

#access cell meta data
head(pData(sasc)) #pData, phenotype data
tail(colData(sasc)) #colData preserved

summary(sasc$n.umi) #data has not been cleaned, everything over 100 umis

table(sasc$sample) #sample is a deliberate mix of nuclei from two donor's heart tissue
table(sasc$donor) #we have used SNPs captured in the RNA-seq reads to distinguish (SoupOrCell)
table(sasc$batch) #we have data from two snRNA-seq assays: Parse and Scale

#access gene meta data
head(fData(sasc))
tail(rowData(sasc))

################################################################################
########### QC #################################################################
################################################################################

#calculate mito content
fData(sasc)$MT <- grepl("MT-", rowData(sasc)$gene_short_name)
table(fData(sasc)$MT)

pData(sasc)$MT_reads <- Matrix::colSums(exprs(sasc)[fData(sasc)$MT,])
pData(sasc)$MT_perc <- pData(sasc)$MT_reads/Matrix::colSums(exprs(sasc))*100
summary(sasc$MT_perc)

#plot some basic QC metrics with ggplot
ggplot(data.frame(pData(sasc)), aes(x=n.umi, y=num_genes_expressed)) +
  facet_wrap(~batch, nrow = 1) +
  geom_point(size=0.5, alpha=0.3, aes(color=MT_perc)) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size = 14),
        aspect.ratio = 1) +
  scale_color_viridis() +
  xlab("UMIs") +
  ylab("Number of Genes Captured") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope=1, color="grey") +
  geom_hline(yintercept = 200, linetype="dotted", color="red") + #change to thresholds
  geom_vline(xintercept = 300, linetype="dotted", color="red") #change to thresholds

#mito content differing between assays; however, samples were the same!
#the pipeline or the assay the cause?
ggplot(data.frame(pData(sasc)), aes(x=donor, y=MT_perc)) +
  facet_wrap(~batch, nrow=1, drop=TRUE, scales="free_x") +
  geom_violin(aes(fill=batch)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size = 14)) +
  xlab("") +
  ylab("MT %") +
  theme(legend.position="none")

#'chattiness' observed in lower quality data from Parse
ggplot(data.frame(pData(sasc)), aes(x=n.umi, y=num_genes_expressed)) +
  facet_wrap(~batch, nrow = 1) +
  geom_point(size=0.5, alpha=0.3, aes(color=organism)) +
  theme_light() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size = 14),
        aspect.ratio = 1) +
  scale_color_viridis(discrete=T) +
  xlab("UMIs") +
  ylab("Number of Genes Captured") +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline(slope=1, color="grey") +
  geom_hline(yintercept = 200, linetype="dotted", color="blue") + #change to thresholds
  geom_vline(xintercept = 300, linetype="dotted", color="blue")

table(sasc$batch, sasc$organism)

#below 300 UMIs, hard to genetically differentiate
ggplot(data.frame(pData(sasc)), aes(x=donor, y=n.umi)) +
  facet_wrap(~batch, nrow=1, drop=TRUE, scales="free_x") +
  geom_violin(aes(fill=batch)) +
  geom_boxplot(notch=T, fill="white", width=0.25, alpha=0.3, outlier.shape=NA) +
  theme_light() +
  geom_hline(yintercept = 300, linetype="dashed", color="blue") + #change to thresholds
  theme(axis.text.x=element_text(angle=45, hjust=1, size = 14),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size = 14)) +
  scale_y_log10() +
  xlab("") +
  ylab("UMIs") +
  theme(legend.position="none")

data.frame(pData(sasc)) %>% 
  group_by(batch, donor) %>%
  summarise(min=min(n.umi),
            q1=quantile(n.umi, probs=c(0.25)),
            med=median(n.umi),
            mean=mean(n.umi),
            q3=quantile(n.umi, probs=c(0.75)),
            max=max(n.umi)) %>%
  data.frame()

#Let's remove everything under 300 UMIs, 200 genes, and over 10% mito
sasc$qcflag <- ifelse(sasc$n.umi >= 300 & sasc$num_genes_expressed >= 200 & sasc$MT_perc < 10, "PASS", "FAIL")
table(sasc$qcflag)

sasc <- sasc[,sasc$qcflag == "PASS"] #filter out failing barcodes
sasc #dim: 36601 24820

#Let's remove things mapping to the wrong genome, those nonassignable to a donor, or deemed a doublet
data.frame(pData(sasc)) %>% 
  group_by(batch) %>%
  count(organism) %>%
  data.frame()

data.frame(pData(sasc)) %>% 
  group_by(batch) %>%
  count(donor) %>%
  data.frame()

sasc <- sasc[,sasc$organism == "human" & sasc$donor %in% c(0, 1)]
sasc #dim: 36601 24216 

#our semi-cleaned data!
table(sasc$batch)
table(sasc$batch, sasc$donor)

#if you fell behind, you can read in the cleaned dataset
#saveRDS(sasc, "BBI_heart_hs_mix_36601humangenes_24216barcodes_12meta.RDS")
#sasc <- readRDS("BBI_heart_hs_mix_36601humangenes_24216barcodes_12meta.RDS")

################################################################################
########### PREPROCESS  ########################################################
################################################################################

#remove non-expressed/non-captured genes 
hist(fData(sasc)$num_cells_expressed)
table(fData(sasc)$num_cells_expressed > 25)

sasc <- sasc[fData(sasc)$num_cells_expressed > 25, ] #filter out genes not expressed in at least 25 cells
sasc #dim: 21279 24216 

#standard mini workflow in monocle3
?estimate_size_factors
sasc <- estimate_size_factors(sasc)

?preprocess_cds
set.seed(1000)
sasc <- preprocess_cds(sasc) #this may take a few minutes

?plot_pc_variance_explained
plot_pc_variance_explained(sasc)

?reduce_dimension
set.seed(1000)
sasc <- reduce_dimension(sasc) #this may take a few minutes
sasc #note the reducedDims!

?plot_cells
plot_cells(sasc, color_cells_by= "batch")
plot_cells(sasc, color_cells_by= "donor")

#group cells into clusters
#you can spend a lot of time tweaking clustering (and should!)
?cluster_cells

#for large datasets, I often find the defaults to be too fine-resolution
#as a starting point, I use the following and round to a nice even number
#it is advisable to try several k values (and/or resolutions)
ceiling(sqrt(dim(sasc)[2])*0.25) # 39 -> 40

#run the function
sasc <- cluster_cells(sasc, k=40, cluster_method="leiden", random_seed=1000) #this may take a few minutes
table(clusters(sasc))

colData(sasc)$k40_leiden_clusters = clusters(sasc) # add cluster information to colData for each cell 

plot_cells(sasc) #by default it will plot the clusters

#what about louvain clustering?
sasc <- cluster_cells(sasc, k=40, cluster_method="louvain", random_seed=1000) #this may take a few minutes
#note: random_see does not seem to be working -> bug report to Monocle3?
table(clusters(sasc)) 

plot_cells(sasc)

# top marker genes
?top_markers
top_marker_genes <- top_markers(sasc, group_cells_by="k40_leiden_clusters") #this may take a few minutes

keep <- top_marker_genes %>%
  filter(fraction_expressing >= 0.30) %>%
  group_by(cell_group) %>%
  top_n(3, marker_score) %>% #can change how many genes we look at as well as criteria to select
  pull(gene_short_name) %>%
  unique()

plot_genes_by_group(sasc,
                    c(keep),
                    group_cells_by="k40_leiden_clusters", #"partition", "cluster"
                    ordering_type="maximal_on_diag",
                    max.size=3)

#saving/catch-up poit
#saveRDS(sasc, "BBI_heart_hs_mix_36601humangenes_24216barcodes_14meta.RDS")
#sasc <- readRDS("BBI_heart_hs_mix_36601humangenes_24216barcodes_14meta.RDS")

################################################################################
########### DATA VIZ  ##########################################################
################################################################################
#move out of monocle3 for plotting flexibility
#add UMAP coordinates to the colData for easy plotting
sasc$UMAP1 <- reducedDim(sasc, "UMAP")[,1]
sasc$UMAP2 <- reducedDim(sasc, "UMAP")[,2]

#generate a distinguishable color scheme
set.seed(1000)
colpal <- randomcoloR::distinctColorPalette(k=12)

#Plot! Try changing some variables! (live coding)
ggplot(data.frame(pData(sasc)), aes(x=UMAP1, y=UMAP2, color=cell_type)) + # try coloring by cell_type
  facet_wrap(~batch+donor) + # try faceting on batch, donor, and batch+donor
  geom_point(size=0.5, alpha=0.5) +
  theme_bw() + #try changing themes, e.g. theme_light()
  scale_color_manual(values=colpal) +
  theme(legend.position="bottom", #move the legend
        aspect.ratio = 1,
        panel.grid=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=8, alpha=1)))

#leiden clusters 3, 4, 2, and 9 appear to be related cell types (ventricular cardiomyocytes)

#CLEARLY THERE IS BOTH TECHNICAL (BATCH) AND BIOLOGICAL (DONOR) VARIATION 

################################################################################
########### Quantifying a batch effect #########################################
################################################################################

#kBET - k-nearest neighbour batch effect test
?kBET
data <- reducedDim(sasc) #test is slow so running on PCA
batch <- sasc$batch
subset_size <- 0.1 #subsample to 10% of the data for speed
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
set.seed(1000)
batch.estimate <- kBET(data[subset_id,], batch[subset_id]) #this may take a few minutes

#By default, it plots. Can extract results if you save it.
str(batch.estimate)
batch.estimate$summary

#convince yourself this is real
set.seed(1000)
randombatch <- sample(sasc$batch, dim(sasc)[2])
batch.estimate.fake <- kBET(data[subset_id,], randombatch[subset_id]) #this may take a few minutes
#okay, I believe the 'expected'

#Do we have a batch effect? Yes, clearly!

################################################################################
################## BATCH CORRECTION ############################################
################################################################################

### Batch correction built-in to Monocle3
# Wrapper around "reducedMNN" function in the Batchelor package from Marioni Lab
?align_cds
?reducedMNN
set.seed(1000)
bc_cds <- align_cds(sasc, alignment_group = "batch", k=50) #this may take a minute

bc_cds #note the "Aligned" in reducedDim

reducedDim(bc_cds, "Aligned")[1:10, 1:2]
reducedDim(sasc, "PCA")[1:10, 1:2]

# Note that is did not modify the count matrix
# (The original MNN correction method did.)
identical(counts(sasc), counts(bc_cds))
identical(reducedDim(sasc, "PCA"), reducedDim(bc_cds, "PCA"))
identical(reducedDim(sasc, "UMAP"), reducedDim(bc_cds, "UMAP"))

# We can run reduce_dimensions to generate a UMAP from the 'aligned' PCA
set.seed(1000)
bc_cds <- reduce_dimension(bc_cds, reduction_method = "UMAP", preprocess_method = "Aligned")
identical(reducedDim(sasc, "UMAP"), reducedDim(bc_cds, "UMAP"))
plot_cells(bc_cds, color_cells_by = "batch")

sasc$aligned_UMAP1 <- reducedDim(bc_cds, "UMAP")[,1] #save these in our original cds
sasc$aligned_UMAP2 <- reducedDim(bc_cds, "UMAP")[,2] #save these in our original cds

#Use kBET to quantitatively ask if it removes the batch effect
data <- reducedDim(bc_cds, "UMAP") #note that we are running this on the UMAP
batch <- bc_cds$batch
subset_size <- 0.1 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
set.seed(1000)
batch.estimate <- kBET(data[subset_id,], batch[subset_id]) #this may take a few minutes
batch.estimate$summary

#While the UMAP looks much better, the kBET metric is telling us that there is
#still a batch effect. If you run kBET on the 'Aligned' PCsm the rejection rate
#is close to one. May we have introduced artifacts?

### Batch correction with Harmony
?RunHarmony
set.seed(1000)
harm_cds <- RunHarmony(sasc, 'batch') #this may take a few minutes

harm_cds #note the "HARMONY" in reducedDim
sasc

reducedDim(harm_cds, "HARMONY")[1:10, 1:2]
reducedDim(sasc, "PCA")[1:10, 1:2]

# Note that is did not modify the count matrix
identical(counts(sasc), counts(harm_cds))
identical(reducedDim(sasc, "PCA"), reducedDim(harm_cds, "PCA"))
identical(reducedDim(sasc, "UMAP"), reducedDim(harm_cds, "UMAP"))

# Monocle3 won't play nice with the other reduced dimensions 
reduce_dimension(harm_cds, reduction_method = "UMAP", preprocess_method = "HARMONY")

# So let's do it ourselves. Under the hood, Monocle3 is using the uwot package
# to genereate UMAPs
?umap
harmony_umap <- umap(reducedDim(harm_cds, "HARMONY"), seed=1000) #this may take a minute
sasc$harmony_UMAP1 <- harmony_umap[,1] #save these to our original cds
sasc$harmony_UMAP2 <- harmony_umap[,2] #save these to our original cds

#Plot
ggplot(data.frame(pData(sasc)), aes(x=harmony_UMAP1, y=harmony_UMAP2, color=batch)) + #cycle through the UMAPs
  geom_point(size=0.5, alpha=0.5) +
  theme_bw() + 
  scale_color_viridis(discrete=T, begin=0.1, end=0.9, option="A") +
  theme(legend.position="bottom", 
        aspect.ratio = 1,
        panel.grid=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=8, alpha=1)))

#Use kBET to quantitatively ask
data <- harmony_umap #note, we could alternatively run kBET at the level of the corrected PCs 
batch <- harm_cds$batch
subset_size <- 0.1 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
set.seed(1000)
batch.estimate <- kBET(data[subset_id,], batch[subset_id]) #this may take a few minutes
batch.estimate$summary

#Once again, the UMAP looks better but the kBET metric suggests that it is still
#far from perfect, especially when run on the corrected PCs. 

#So what benefit does batch correction offer??? 

#That is debatable, but certainly one thing it can do is help in identifying
#cell types. We gave you a cds with cell types already annotated. But what if
#you did not know this?

#cluster cells that have been aligned and plot these on our original UMAP
bc_cds <- cluster_cells(bc_cds, k=40, cluster_method="leiden", random_seed=1000) #this may take a few minutes

table(clusters(bc_cds))
table(partitions(bc_cds))

sasc$aligned_clusters <- clusters(bc_cds) #save to original cds object
sasc$aligned_partitions <- partitions(bc_cds) #save to original cds object

ggplot(data.frame(pData(sasc)), aes(x=UMAP1, y=UMAP2, color=k40_leiden_clusters)) + #also color by aligned_clusters & aligned_partitions
  geom_point(size=0.5, alpha=0.5) +
  theme_bw() + 
  scale_color_manual(values=colpal) +
  theme(legend.position="bottom", 
        aspect.ratio = 1,
        panel.grid=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=8, alpha=1)))

data.frame(colData(sasc)) %>%
  group_by(cell_type) %>%
  count(aligned_partitions) %>%
  spread(aligned_partitions,n)

# Even with imperfect batch correction methods, we can start to see that our 
# original clusters 1, 3,and 7 are related (ventricular cardiomyocytes).

#saving/catch-up poit
#saveRDS(sasc, "BBI_heart_hs_mix_36601humangenes_24216barcodes_22meta.RDS")
#sasc <- readRDS("BBI_heart_hs_mix_36601humangenes_24216barcodes_22meta.RDS")

################################################################################
########## Differential Gene Expression ########################################
################################################################################

#We are going to focus specifically on the ventricular cardiomyocytes for the remainder of today

#combine batch and donor as a new column 
sasc$id <- paste(sasc$batch, sasc$donor, sep="_")
head(colData(sasc))

#subset ventricular cardiomyocyte data only 
table(sasc$cell_type)
cds_vent <- sasc[,sasc$cell_type == "Ventricular Cardiomyocyte" & sasc$k40_leiden_clusters %in% c("3", "4", "9", "2")]
cds_vent #dim: 21279 10337 

#sanity Check 
table(cds_vent$cell_type)

#a look at some cherry-picked genes
plot_cells(cds_vent, genes=c("LINC00486", "LINC-PINT", "FN1", "XIST"), scale_to_range=F)

#DE can take a long time, so we will just run this on a subset of pre-defined genes
gene_list <- c("LINC00486", "TTN", "LINC-PINT", "TAS2R14", "MT-CO1", #top markers
               "MT-ND4", "FN1", "LAMA2", "XIST", "PDK4", "ZBTB16", #top markers
               "PPP1R3E", "TMTC1", "NT5DC3", "RBX1", "MRPL45", "ESR2", "TUBGCP4", #selected to be uninteresting
               "MYH7", "MYL2", "MB", "ACTC1", "TPM1", "MYH6") #classical CM genes

cds_subset <- cds_vent[rowData(cds_vent)$gene_short_name %in% gene_list,]
cds_subset #dim: 24 10298 

#plot expression levels of these genes split by donor
plot_genes_violin(cds_subset, group_cells_by="donor", ncol=4) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

?fit_models
donor_model <- fit_models(cds_subset,
                           model_formula_str = "~donor",
                           expression_family="negbinomial")

?coefficient_table
coefficient_table(donor_model) %>% 
  filter(term == "donor1") %>%
  filter(q_value < 0.05) %>%
  select(id, gene_short_name, term, q_value, estimate) %>%
  arrange(gene_short_name)

#controlling for batch effects: the best way in my opinion :) 
donor_batch_model <- fit_models(cds_subset,
                                 model_formula_str = "~donor + batch",
                                 expression_family="negbinomial")

coefficient_table(donor_batch_model) %>% 
  filter(term == "donor1") %>%
  filter(q_value < 0.05) %>%
  select(id, gene_short_name, term, q_value, estimate) %>%
  arrange(gene_short_name)

#comparing models of gene expression
?compare_models
compare_models(donor_batch_model, donor_model) %>% select(gene_short_name, q_value) %>% data.frame()

#plot outside monocle
cntmtx <- normalized_counts(cds_subset)
cds_subset$LINC00486 <- cntmtx["ENSG00000230876",]
cds_subset$TTN <- cntmtx["ENSG00000155657",]
cds_subset$XIST <- cntmtx["ENSG00000229807",]

#the case of LINC00486
a <- ggplot(data.frame(pData(cds_subset)), aes(x=donor, y=LINC00486)) + 
  geom_violin(aes(fill=donor)) +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

b <- ggplot(data.frame(pData(cds_subset)), aes(x=batch, y=LINC00486)) + 
  geom_violin(aes(fill="salmon")) +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

c <- ggplot(data.frame(pData(cds_subset)), aes(x=id, y=LINC00486, fill=donor)) + 
  geom_violin() +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

plot_grid(a + theme(legend.position="none"),
          b + theme(legend.position="none") + ylab(""),
          c + theme(legend.position="none") + ylab(""), 
          labels=c("A", "B", "C"),
          nrow=1)

coefficient_table(donor_model) %>% 
  filter(gene_short_name == "LINC00486" & term =="donor1") %>%
  select(id, gene_short_name, term, q_value, estimate) 

coefficient_table(donor_batch_model) %>% 
  filter(gene_short_name == "LINC00486" & term %in% c("donor1", "batchScale")) %>%
  select(id, gene_short_name, term, q_value, estimate) 

#the case of TNN
x <- ggplot(data.frame(pData(cds_subset)), aes(x=donor, y=TTN)) + 
  geom_violin(aes(fill=donor)) +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

y <- ggplot(data.frame(pData(cds_subset)), aes(x=batch, y=TTN)) + 
  geom_violin(aes(fill="salmon")) +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

z <- ggplot(data.frame(pData(cds_subset)), aes(x=id, y=TTN, fill=donor)) + 
  geom_violin() +
  geom_boxplot(width=0.2, fill="white", alpha=0.3, outlier.shape=NA) +
  theme_bw() 

plot_grid(x + theme(legend.position="none"),
          y + theme(legend.position="none") + ylab(""),
          z + theme(legend.position="none") + ylab(""), 
          labels=c("A", "B", "C"),
          nrow=1)

coefficient_table(donor_model) %>% 
  filter(gene_short_name == "TTN" & term =="donor1") %>%
  select(id, gene_short_name, term, q_value, estimate) 

coefficient_table(donor_batch_model) %>% 
  filter(gene_short_name == "TTN" & term %in% c("donor1", "batchScale")) %>%
  select(id, gene_short_name, term, q_value, estimate) 

################################################################################
########################### KEY TAKEAWAYS ######################################
################################################################################
################################################################################
### DO NOT BLINDLY APPLY BATCH CORRECTION! #####################################
### (YOU RISK INTRODUCING MORE ARTIFACTS THAN YOU REMOVE) ######################
################################################################################
################################################################################
### WOULD RECOMMEND ONLY USING 'BATCH CORRECTION METHODS' FOR VIZUALIZATION ####
### AND (CAREFUL) USE IN DATA EXPLORATION, CELL TYPE ANNOTATION, ETC ###########
################################################################################
################################################################################
### IF YOU HAVE A BATCH EFFECT, MODEL THIS IN YOUR DOWNSTREAM ANALYSES #########
### E.G. INCLUDE AS A COVARIATE IN DE ANALYSIS #################################
################################################################################
################################################################################
############################### THE END ########################################
################################################################################
