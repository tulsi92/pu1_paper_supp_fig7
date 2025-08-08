################
## 2025-08-07 ##
################

## Date of code and file generation
## 2024-11-14
## Code used to reprocess Olah et al (2020) microglia scRNA-seq data and generate plots
## This version as shown before with just subsetting to remove non-microglial cells (their classification)
## no reclustering after this

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

packageVersion("Seurat")

############################################
## Process using old code from Olah paper ##
############################################
###Step 1: load data, remove genes often associated with technical artifacts, and filter by total UMI count###
alldat = read.csv("~/Downloads/20241031.PU1_olah/SupplementaryData14.csv",as.is=T,row.names=1,header=T)
alldat = alldat[grep("^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS",rownames(alldat),invert=T),]
alldat = alldat[,which(colSums(alldat)>=1000)]
alldat_cpm = sweep(alldat,2,colSums(alldat),"/")*10^6
olah2020_annot <- read_tsv("~/Downloads/20241031.PU1_olah/cluster_annotation.txt")
batchval = read.csv("~/Downloads/20241031.PU1_olah/41467_2020_19737_MOESM17_ESM.csv",header=T,row.names=1)
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

datobj_olah = NormalizeData(datobj_olah)
# datobj_olah = ScaleData(datobj_olah,vars.to.regress=c("nCount_RNA","batch"))
datobj_olah = ScaleData(datobj_olah,vars.to.regress=c("batch", "nUMI"), features = rownames(datobj_olah)) #nUMI not in object
rownames(datobj_olah@assays$RNA@layers$data) = rownames(alldat)
datobj_olah <- FindVariableFeatures(object = datobj_olah, selection.method = "vst", verbose = TRUE)
vargenes=rownames(datobj_olah@assays$RNA@layers$data)[which(apply(datobj_olah@assays$RNA@layers$data,1,var)>apply(datobj_olah@assays$RNA@layers$data,1,mean))]
datobj_olah <- RunPCA(object = datobj_olah, pc.genes = vargenes, do.print = FALSE)
datobj_olah <- RunUMAP(object = datobj_olah, reduction = 'pca', dims = 1:20)
datobj_olah <- FindNeighbors(object = datobj_olah, reduction = "pca", dims = 1:20)
## combinations of principal components (5 to 15) and resolution parameters (0.2, 0.4, 0.6, and 0.8) tested
num_pcs=15
resval=0.59
datobj_olah <- FindClusters(datobj_olah, dims = 1:num_pcs, algorithm=1, resolution = resval)
clutstids=as.numeric(datobj_olah@active.ident)

## add clustering labels from paper - seurat_clusters is wrong because this is new clustering so use annotation column and reassign seurat_clusters
DimPlot(datobj_olah, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) + theme(legend.position = "right")

datobj_olah@meta.data$subtype = ""
for(j in unique(batchval_clusters$cluster)){
  cl_type = batchval_clusters[batchval_clusters$cluster==j,]; 
  datobj_olah@meta.data$subtype[datobj_olah@meta.data$cluster_label == j] = as.character(cl_type$annotation[j])
}

Idents(datobj_olah) <- "subtype"
# olah_metadata <- datobj_olah@meta.data

table(datobj_olah$annotation, datobj_olah$seurat_clusters)
DimPlot(datobj_olah, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, group.by = "seurat_clusters") + theme(legend.position = "right")

#########################################
## Subset for microglial clusters only ##
#########################################
Idents(datobj_olah) <- "annotation"

datobj_olah_mic <- subset(datobj_olah, cells = microglia_cells)

# datobj_olah_mic <- subset(datobj_olah, idents = c(
#   "1_Homeostatic",
#   "2_Homeostatic",
#   "3_Stress_related",
#   "4_Interferon_response",
#   "5_Anti_inflammatory",
#   "6_Anti_inflammatory",
#   "7_Antigen_presentation",
#   "8_Unknown",
#   "9_Cell_cycle"))

DimPlot(datobj_olah_mic, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, group.by = "annotation") + theme(legend.position = "right")
DimPlot(datobj_olah_mic, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE, group.by = "seurat_clusters") + theme(legend.position = "right")

##################
## Marker genes ##
##################
high_in_dam_genes_2 <- c("APOE","CLEC7A","LPL","CST7","CD9","SPP1","FABP5")
homeostatic_genes_2 <- c("CSF1R","CD33","CX3CR1","AIF1","TMEM119","P2RY12")
lymphoid_genes_2 <- c("CD28","PDCD1","CD72","CD52","CD69","CXCR4","PRKCQ","CD5","CTLA2A","PF4","DKK2","TNFRSF13B","TNFSF13B","LAG3")

## microglial subset cell type markers from scrnaseq papers
genesets_from_papers_consensus <- read_xlsx("~/Downloads/20241031.PU1_olah/20241003.Myeloid_Genesets_JMC.xlsx", sheet = "genesets")
genesets_from_papers <- read_xlsx("~/Downloads/20241031.PU1_olah/20240903_Myeloid_Genesets_APD.xlsx", sheet = "genesets")
genesets_from_papers$Consensus_DAM_signature <- genesets_from_papers_consensus$Consensus_DAM_signature

## convert table of genesets to genelists
genelists <- list()
for(i in 1:ncol(genesets_from_papers)) {   
  genelists[[i]] <- genesets_from_papers[ , i]
}; names(genelists) <- colnames(genesets_from_papers) 

genelists_clean <- lapply(genelists, function(col)col[!is.na(col)])

###################################
## Plot markers for all clusters ##
###################################
Idents(datobj_olah_mic) <- "seurat_clusters"

DotPlot(datobj_olah_mic, features = c("SPI1", high_in_dam_genes_2, homeostatic_genes_2), assay = "RNA", scale.by = "size") +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) +
  labs(x = "", y = "") +
  coord_fixed(ratio = 1.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(family = 'Helvetica', size = 10),
    legend.text = element_text(family = 'Helvetica', size = 10),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm"))

DotPlot(datobj_olah_mic, features = c("SPI1", lymphoid_genes_2), 
        assay = "RNA", scale.by = "size", scale.max = 40) +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) +
  labs(x = "", y = "") +
  coord_fixed(ratio = 1.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(family = 'Helvetica', size = 10),
    legend.text = element_text(family = 'Helvetica', size = 10),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

DotPlot_scCustom(datobj_olah_mic, features = genelists_clean$Consensus_DAM_signature, assay = "RNA", 
                 scale.by = "size") +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) +
  labs(x = "", y = "")

#############################
## Heatmap of marker genes ##
#############################
expression_data <- as.matrix(GetAssayData(datobj_olah_mic, layer = "RNA", slot = "scale.data"))
clusters <- Idents(datobj_olah_mic)

# genes_of_interest <- unique(c("TMEM119","P2RY12", homeostatic_genes_2, genelists_clean$Consensus_DAM_signature))
# genes_of_interest <- intersect(rownames(expression_data), c(homeostatic_genes_2, high_in_dam_genes_2))
# genes_of_interest <- intersect(rownames(expression_data), unique(c(homeostatic_genes_2,
#                                                                    # genelists_clean$Gazestani_cl_HM_3_UP,
#                                                                    genelists_clean$Gazestani_cl_HM_4_UP)))

expression_subset <- expression_data[genes_of_interest, , drop = FALSE]
expression_df <- as.data.frame(t(expression_subset))  # Transpose to make cells rows and genes columns
expression_df$Cluster <- clusters
## Calculate average expression per cluster for each gene
avg_expr_per_cluster <- expression_df %>%
  group_by(Cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))  # Calculate mean expression per gene

## Reshape the data to a matrix (genes as rows, clusters as columns)
collapsed_expr_matrix <- avg_expr_per_cluster %>%
  column_to_rownames(var = "Cluster") %>%
  as.matrix()

Heatmap(collapsed_expr_matrix,
        name = "Average Expression",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        column_names_rot = 45,
        column_title = "Clusters",
        row_title = "Genes")

###################################################
## Recluster cells after removing non-microglial ##
###################################################
microglia_cells <- batchval_clusters$cell_id[which(batchval_clusters$cluster_label %notin% c(10,11,12,13,14))]
nonmic_cells <- batchval_clusters$cell_id[which(batchval_clusters$cluster_label %in% c(10,11,12,13,14))]

## Removed 146 non-microglial cells
a = DimPlot(object = datobj_olah, reduction = "umap", group.by = "cluster_label")
b = DimPlot(object = datobj_olah_mic, reduction = "umap", group.by = "cluster_label", cells = microglia_cells)

cowplot::plot_grid(a, b)

datobj_olah_mic_reclust <- subset(datobj_olah_mic, cells = microglia_cells)

## Recluster cells
datobj_olah_mic_reclust <- RunPCA(datobj_olah_mic_reclust)
datobj_olah_mic_reclust <- FindNeighbors(datobj_olah_mic_reclust, dims = 1:30)
datobj_olah_mic_reclust <- FindClusters(datobj_olah_mic_reclust, resolution = 0.4)
datobj_olah_mic_reclust <- RunUMAP(object = datobj_olah_mic_reclust, dims = 1:30)

DimPlot(object = datobj_olah_mic_reclust, reduction = "umap")

table(datobj_olah_mic_reclust$seurat_clusters, datobj_olah_mic_reclust$annotation)
