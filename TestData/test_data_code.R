## 2025-08-07
## Uploading this to github with sample code for subset object

################
## 2025-01-29 ##
################
## Code used to reprocess Olah et al (2020) microglia scRNA-seq data and generate plots
# This version is processed using Seurat v5.0.2 and subset to remove non-microglial cells 
# (authors annotation) and then reclustered

## Load packages
library(tidyverse)
library(readr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(viridis)
library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)
library(cowplot)
library(scCustomize)
library(gridExtra)
library(SCpubr)
library(hypeR)
library(pheatmap)

packageVersion("Seurat")
set.seed(42)

datobj_olah_mic <- readRDS("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2024-10-31.Human_replication/archive/olah2020_microglia_reprocessed_mic_clusters_only.rds")

setwd("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/")

################
## 2025-08-07 ##
################
## Randomly sample columns from csv to make test dataset
alldat = read.csv("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2025-08-07 Github upload/SupplementaryData14.csv",as.is=T,row.names=1,header=T)
sampled_columns <- sample(colnames(alldat), 3000, replace = FALSE)
alldat = alldat[grep("^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS",rownames(alldat),invert=T), sampled_columns]

## Randomly sample cells from object to make test dataset
# sampled_cells <- sample(colnames(olah2020_microglia_reprocessed_012925), 3000, replace = FALSE)
# olah2020_microglia_reprocessed_012925_subset = subset(olah2020_microglia_reprocessed_012925, cells = sampled_cells)
# DimPlot(olah2020_microglia_reprocessed_012925_subset, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, group.by = "our_cluster_assignment") + theme(legend.position = "right")

################
## 2024-12-17 ##
################
## Object processing

################################################
## Process using similar code from Olah paper ##
################################################
###Step 1: load data, remove genes often associated with technical artifacts, and filter by total UMI count###
# alldat = read.csv("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2025-08-07 Github upload/SupplementaryData14.csv",as.is=T,row.names=1,header=T)
# alldat = alldat[grep("^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS",rownames(alldat),invert=T),]
alldat = alldat[,which(colSums(alldat)>=1000)]
alldat_cpm = sweep(alldat,2,colSums(alldat),"/")*10^6
olah2020_annot <- read_tsv("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2025-08-07 Github upload/cluster_annotation.txt")
batchval = read.csv("/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2024-10-31.Human_replication/41467_2020_19737_MOESM17_ESM.csv",header=T,row.names=1)
batchval_clusters <- batchval %>% 
  rownames_to_column("cell_id") %>%
  left_join(olah2020_annot, by = c("cluster_label" = "cluster"))
rownames(batchval_clusters) <- batchval_clusters$cell_id

head(batchval_clusters)

###Step 2: run first-level clustering, ranging over parameters
datobj_olah=CreateSeuratObject(counts = alldat, project = "olah", meta.data = batchval_clusters)
olah_clusters <- as.factor(batchval_clusters$cluster_label)
names(olah_clusters) <- rownames(batchval_clusters)
# SetIdent(datobj_olah, value = olah_clusters)

# Idents(datobj_olah) = "annotation"
Idents(datobj_olah) = "cluster_label"
head(Idents(datobj_olah))

datobj_olah <- NormalizeData(datobj_olah)
datobj_olah <- ScaleData(datobj_olah,vars.to.regress=c("nCount_RNA", "batch"), features = rownames(datobj_olah)) #nUMI not in object
rownames(datobj_olah@assays$RNA@data) = rownames(alldat)
vargenes=rownames(datobj_olah@assays$RNA@layers$data)[which(apply(datobj_olah@assays$RNA@layers$data,1,var)>apply(datobj_olah@assays$RNA@layers$data,1,mean))]
####
datobj_olah <- FindVariableFeatures(object = datobj_olah, selection.method = "vst", verbose = TRUE)
datobj_olah <- RunPCA(object = datobj_olah, pc.genes = vargenes, do.print = FALSE)
datobj_olah <- RunUMAP(object = datobj_olah, reduction = 'pca', dims = 1:20)
datobj_olah <- FindNeighbors(object = datobj_olah, reduction = "pca", dims = 1:20)
## combinations of principal components (5 to 15) and resolution parameters (0.2, 0.4, 0.6, and 0.8) tested
num_pcs=15
resval=0.6
datobj_olah <- FindClusters(datobj_olah, dims = 1:num_pcs, algorithm=1, resolution = resval)
clustids=as.numeric(datobj_olah@active.ident)

## add clustering labels from paper - seurat_clusters is wrong because this is new clustering so use annotation column and reassign seurat_clusters
DimPlot(datobj_olah, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) + theme(legend.position = "right")

datobj_olah@meta.data$subtype = ""
for(j in unique(batchval_clusters$cluster)){
  cl_type = batchval_clusters[batchval_clusters$cluster==j,]; 
  datobj_olah@meta.data$subtype[datobj_olah@meta.data$cluster_label == j] = as.character(cl_type$annotation[j])
}

Idents(datobj_olah) <- "annotation"

table(datobj_olah$annotation, datobj_olah$seurat_clusters)
DimPlot(datobj_olah, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, group.by = "annotation") + theme(legend.position = "right")

###################################################
## Recluster cells after removing non-microglial ##
###################################################
microglia_cells <- batchval_clusters$cell_id[which(batchval_clusters$cluster_label %notin% c(10,11,12,13,14))]
nonmic_cells <- batchval_clusters$cell_id[which(batchval_clusters$cluster_label %in% c(10,11,12,13,14))]

## Removed 146 non-microglial cells
a <- DimPlot(object = datobj_olah, reduction = "umap", group.by = "cluster_label")
b <- DimPlot(object = datobj_olah, reduction = "umap", group.by = "cluster_label", cells = microglia_cells)

cowplot::plot_grid(a, b)

our_cells <- colnames(olah2020_microglia_reprocessed_012925)

datobj_olah_mic_reclust <- subset(datobj_olah, cells = our_cells)

## Recluster cells
datobj_olah_mic_reclust <- NormalizeData(datobj_olah_mic_reclust)
datobj_olah_mic_reclust <- ScaleData(datobj_olah_mic_reclust)
datobj_olah_mic_reclust <- RunPCA(datobj_olah_mic_reclust, features = vargenes)
datobj_olah_mic_reclust <- FindNeighbors(datobj_olah_mic_reclust, dims = 1:30)
datobj_olah_mic_reclust <- FindClusters(datobj_olah_mic_reclust, resolution = resval)
datobj_olah_mic_reclust <- RunUMAP(object = datobj_olah_mic_reclust, dims = 1:30)

a = DimPlot(object = datobj_olah_mic_reclust, reduction = "umap")
b = DimPlot(object = datobj_olah_mic, reduction = "umap")

cowplot::plot_grid(a, b)

table(datobj_olah_mic_reclust$seurat_clusters, datobj_olah_mic_reclust$annotation)
# table(datobj_olah_mic$seurat_clusters, datobj_olah_mic$annotation)
table(datobj_olah_mic_reclust$seurat_clusters, olah2020_microglia_reprocessed_012925_subset$seurat_clusters)

################
## 2025-01-29 ##
################
datobj_olah_mic <- olah2020_microglia_reprocessed_012925

our_cell_assignments <- datobj_olah_mic@meta.data |>
  dplyr::select(new_cluster_name) %>%
  dplyr::filter(rownames(.) %in% colnames(datobj_olah_mic_reclust))

## Rename clusters ##
olah_metadata <- datobj_olah_mic_reclust@meta.data %>%
  mutate(new_cluster_name_seu = paste0("MG", as.numeric(seurat_clusters)))

unique(olah_metadata$new_cluster_name)

datobj_olah_mic_reclust <- AddMetaData(datobj_olah_mic_reclust, metadata = olah_metadata$new_cluster_name, col.name = "new_cluster_name")
datobj_olah_mic_reclust <- AddMetaData(datobj_olah_mic_reclust, metadata = olah_metadata$our_cluster_assignment, col.name = "our_cluster_assignment")

# Idents(datobj_olah_mic_reclust) <- "new_cluster_name"
Idents(datobj_olah_mic_reclust) <- "our_cluster_assignment"

## Reorder clusters
cluster_order <- c(
  "MG1",
  "MG2",
  "MG3",
  "MG4",
  "MG5",
  "MG6",
  "MG7",
  "MG8",
  "MG9")

datobj_olah_mic_reclust <- SetIdent(datobj_olah_mic_reclust, value = factor(Idents(datobj_olah_mic_reclust), levels = cluster_order))

DimPlot(datobj_olah_mic_reclust, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, group.by = "our_cluster_assignment") + theme(legend.position = "right")

## Dotplot of combined genelist with reprocessed object from 12-08-24
combined_genelist <- c("SPI1", "APOE","CLEC7A","CST7","SPP1",
                       "CSF1R","CX3CR1","TMEM119","P2RY12",
                       "CD28","PDCD1","CD5","DKK2")

DotPlot(datobj_olah_mic_reclust, features = combined_genelist, col.min = -1, col.max = 1,
        scale.max = 20) +
  scale_size(range = c(0.5, 5)) +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) +
  labs(x = "", y = "") +
  coord_fixed(ratio = 0.75) +
  geom_hline(yintercept = c(3.5, 6.5), linetype = "dotted", color = "black", size = 0.5) +
  geom_vline(xintercept = c(1.5, 5.5, 9.5), linetype = "dotted", color = "black", size = 0.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(family = 'Helvetica', size = 10),
    legend.text = element_text(family = 'Helvetica', size = 10),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"))

##################
## Marker genes ##
##################
## microglial subset cell type markers from scrnaseq papers
genesets_from_papers_consensus <- read_xlsx("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2024-10-31.Human_replication/20241003.Myeloid_Genesets_JMC.xlsx", sheet = "genesets")
genesets_from_papers <- read_xlsx("~/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2024-10-31.Human_replication/20240903_Myeloid_Genesets_APD.xlsx", sheet = "genesets")
genesets_from_papers$Consensus_DAM_signature <- genesets_from_papers_consensus$Consensus_DAM_signature

## convert table of genesets to genelists
genelists <- list()
for(i in 1:ncol(genesets_from_papers)) {   
  genelists[[i]] <- genesets_from_papers[ , i]
}; names(genelists) <- colnames(genesets_from_papers) 

genelists_clean <- lapply(genelists, function(col)col[!is.na(col)])

homeostatic_genes <- c("MS4A6A", "CSF1R","CD33","CX3CR1","AIF1","TMEM119","P2RY12")

## Very few cells expressing lymphoid marker genes
# FeaturePlot(datobj_olah_mic, features = c("CD28","PDCD1","CD52","PRKCQ","CD5",
#                                           "DKK2"), label = TRUE, label.size = 5, 
#             order = TRUE, repel = TRUE, pt.size = 0.5, cols = c("lightgrey", "red"))

##################################
## Check enrichment of genesets ##
##################################
datobj_olah_mic <- datobj_olah_mic_reclust

## Consensus homeostatic signature
Idents(datobj_olah_mic) <- "our_cluster_assignment"
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(homeostatic_genes, rownames(datobj_olah_mic))), name = "Homeostatic_signature", ctrl = 100, seed = 123456)

FeaturePlot(datobj_olah_mic, features = "Homeostatic_signature1", label = TRUE, repel = TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## Consensus DAM signature
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(genelists_clean$Consensus_DAM_signature, rownames(datobj_olah_mic))), name = "Consensus_DAM_signature", ctrl = 100, seed = 123456)

FeaturePlot(datobj_olah_mic, features = "Consensus_DAM_signature1", label = TRUE, repel = TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## Lymphoid signature
lymphoid_genes1 <- c("CD28", "PDCD1", "CD52", "CD72", "CD5")
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(lymphoid_genes1, rownames(datobj_olah_mic))), name = "Lymphoid_signature1_", ctrl = 100, seed = 123456)

lymphoid_genes2 <- c("CD28", "PDCD1", "CD52", "CD72", "CD5", "DKK2", "TNFSF13B")
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(lymphoid_genes2, rownames(datobj_olah_mic))), name = "Lymphoid_signature2_", ctrl = 100, seed = 123456)

lymphoid_genes3 <- c("CD28", "PDCD1", "CD52", "CD72", "CD5", "DKK2", "TNFSF13B", "PRKCQ")
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(lymphoid_genes3, rownames(datobj_olah_mic))), name = "Lymphoid_signature3_", ctrl = 100, seed = 123456)

lymphoid_genes4 <- c("CD28", "PDCD1", "CD52", "CD72", "CD5", "PRKCQ")
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(lymphoid_genes4, rownames(datobj_olah_mic))), name = "Lymphoid_signature4_", ctrl = 100, seed = 123456)

## Selected lymphoid signature
lymphoid_genes5 <- c("CD28", "PDCD1", "CD5", "DKK2")
datobj_olah_mic <- AddModuleScore(datobj_olah_mic, features = list(intersect(lymphoid_genes5, rownames(datobj_olah_mic))), name = "Lymphoid_signature_selected", ctrl = 100, seed = 123456)

###################
## Module scores ##
###################
## Check module scores using reclustered data not annotation from original paper
dam_score <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Consensus_DAM_signature1) |> 
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_dam = mean(Consensus_DAM_signature1, na.rm = TRUE))

hom_score <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Homeostatic_signature1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_hom = mean(Homeostatic_signature1, na.rm = TRUE))

lymphoid_score_1 <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature1_1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_1 = mean(Lymphoid_signature1_1, na.rm = TRUE))

lymphoid_score_2 <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature2_1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_2 = mean(Lymphoid_signature2_1, na.rm = TRUE))

lymphoid_score_3 <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature3_1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_3 = mean(Lymphoid_signature3_1, na.rm = TRUE))

lymphoid_score_4 <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature4_1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_4 = mean(Lymphoid_signature4_1, na.rm = TRUE))

lymphoid_score_selected <- datobj_olah_mic@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature_selected1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_selected = mean(Lymphoid_signature_selected1, na.rm = TRUE))

pu1_expression <- FetchData(datobj_olah_mic, vars = "SPI1", layer = "data") |> #lognorm counts
  rownames_to_column("cell_id") |>
  distinct(cell_id, .keep_all = TRUE) |>
  left_join(datobj_olah_mic@meta.data |> rownames_to_column("cell_id")) |> 
  dplyr::select(our_cluster_assignment, SPI1) |> 
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_pu1_expr = mean(SPI1, na.rm = TRUE))

module_score_comparison <- dam_score |>
  left_join(hom_score) |>
  # left_join(lymphoid_score_1) |>
  # left_join(lymphoid_score_2) |>
  # left_join(lymphoid_score_3) |>
  # left_join(lymphoid_score_4) |>
  left_join(lymphoid_score_selected) |>
  left_join(pu1_expression)

module_score_comparison |>
  mutate(dam_hom_ratio = avg_score_dam / avg_score_hom) |>
  data.table()

## DAM score ratio
module_score_comparison <- module_score_comparison |>
  mutate(dam_hom_ratio = avg_score_dam / avg_score_hom)

###########################################################
## Group clusters based on PU.1/lymphoid gene expression ##
###########################################################
olah_metadata_groups <- datobj_olah_mic@meta.data |>
  mutate(lymphoid_group = case_when(
    our_cluster_assignment %in% c("MG1", "MG2", "MG7") ~ "Homeostatic",
    our_cluster_assignment %in% c("MG6", "MG8", "MG9") ~ "Lymphoid_Positive_DAM",
    our_cluster_assignment %in% c("MG3", "MG4", "MG5") ~ "Lymphoid_Negative_DAM",
    TRUE ~ NA))

unique(olah_metadata_groups$lymphoid_group)

datobj_olah_mic <- AddMetaData(datobj_olah_mic, metadata = olah_metadata_groups$lymphoid_group, col.name = "lymphoid_group")

table(datobj_olah_mic$our_cluster_assignment, datobj_olah_mic$lymphoid_group)

Idents(datobj_olah_mic) <- "lymphoid_group"

## Add SPI1 expression per group
spi1_expr <- FetchData(datobj_olah_mic, vars = "SPI1", slot = "scale.data") |>
  rownames_to_column("cell_id")

olah_metadata_groups <- olah_metadata_groups |>
  rownames_to_column("cell_id") |>
  left_join(spi1_expr)

## UMAP
DimPlot(datobj_olah_mic, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, group.by = "lymphoid_group") + theme(legend.position = "bottom")

## Reorder clusters
cluster_order <- c(
  "Homeostatic",
  "Lymphoid_Positive_DAM",
  "Lymphoid_Negative_DAM")

datobj_olah_mic <- SetIdent(datobj_olah_mic, value = factor(Idents(datobj_olah_mic), levels = cluster_order))

Idents(datobj_olah_mic) <- "lymphoid_group"

DotPlot(datobj_olah_mic, features = combined_genelist, col.min = -1, col.max = 1,
        scale.max = 20) +
  scale_size(range = c(0.5, 5)) +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) +
  labs(x = "", y = "") +
  coord_fixed(ratio = 0.75) +
  geom_vline(xintercept = c(1.5, 5.5, 9.5), linetype = "dotted", color = "black", size = 0.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(family = 'Helvetica', size = 10),
    legend.text = element_text(family = 'Helvetica', size = 10),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"))

# saveRDS(datobj_olah_mic, "/Users/tulsi/Library/Mobile Documents/com~apple~CloudDocs/Projects/PU1 Project/2024-12-08.Human_replication_final_object_decision/olah2020_microglia_reprocessed_012925.rds")

## Compare reprocessed clusters to original published clusters using scmap
library(SummarizedExperiment)
library(scmap)

## Convert to SCE
sce <- as.SingleCellExperiment(datobj_olah_mic)
unique(datobj_olah_mic$cluster_label)

## Split into reference and query based on the annotation
colData(sce)$original_clusters <- datobj_olah_mic$annotation
colData(sce)$new_clusters <- datobj_olah_mic$our_cluster_assignment
ref_sce <- sce
query_sce <- sce

## Set original cluster annotation as the reference label
ref_sce$original_clusters <- colData(ref_sce)$original_clusters

## Set new cluster annotation as the query label
query_sce$new_clusters <- colData(query_sce)$new_clusters
rowData(query_sce)$feature_symbol <- rownames(query_sce)

## Select features to use (all genes)
rowData(ref_sce)$feature_symbol <- rownames(ref_sce)
ref_sce <- selectFeatures(ref_sce, suppress_plot = TRUE)
table(rowData(ref_sce)$feature_symbol)

## Index reference clusters
ref_sce <- indexCluster(ref_sce, cluster_col = "original_clusters")

## Run scmapCluster to map query cells to reference clusters
scmap_results <- scmapCluster(
  projection = query_sce,
  index_list = list(ref = metadata(ref_sce)$scmap_cluster_index))

## Add predicted labels to the query SCE
query_sce$scmap_prediction <- scmap_results$combined_labs

## Confusion matrix between original and new clusters
table(Predicted = query_sce$scmap_prediction, Actual = query_sce$new_clusters)
table(datobj_olah_mic$annotation, datobj_olah_mic$our_cluster_assignment)

plot(
  getSankey(
    ref_sce$original_clusters, 
    query_sce$new_clusters,
    plot_height = 400))

############################
## Confusion matrix plots ##
############################
## Create contingency table and calculate percentages by row
conf_mat <- table(datobj_olah_mic$annotation, datobj_olah_mic$our_cluster_assignment)

conf_df <- as.data.frame(conf_mat) %>%
  group_by(Var1) %>%
  mutate(Percent = 100 * Freq / sum(Freq)) %>%
  ungroup()

## Heatmap
ggplot(conf_df, aes(x = Var2, y = Var1, fill = Percent)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Percent") +
  geom_text(aes(label = sprintf("%.1f", Percent)), size = 3) +
  scale_y_discrete(limits = rev(levels(conf_df$Var1))) +  # REVERSE Y-AXIS
  labs(
    title = "",
    x = "New Olah Clusters",
    y = "Original Olah Clusters") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank())

