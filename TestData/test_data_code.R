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

datobj_olah_mic_subset <- readRDS("olah2020_microglia_reprocessed_mic_clusters_only.rds")

################
## 2025-08-07 ##
################
## Randomly sample columns from csv to make test dataset
alldat = read.csv("SupplementaryData14.csv",as.is=T,row.names=1,header=T)
sampled_columns <- sample(colnames(alldat), 3000, replace = FALSE)
alldat = alldat[grep("^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS",rownames(alldat),invert=T), sampled_columns]

## Randomly sample cells from seurat object to make test dataset
# sampled_cells <- sample(colnames(datobj_olah_mic_subset), 3000, replace = FALSE)
# datobj_olah_mic_subset = subset(datobj_olah_mic_subset, cells = sampled_cells)

datobj_olah_mic_subset <- readRDS("microglia_reprocessed_subset.rds")

our_cluster_assignments <- read_tsv("our_cluster_assignments.txt")

datobj_olah_mic_subset <- AddMetaData(datobj_olah_mic_subset, metadata = our_cluster_assignments$new_cluster_name, col.name = "new_cluster_name")

Idents(datobj_olah_mic_subset) <- "new_cluster_name"

## Reorder clusters
cluster_order <- c(
  "MG1","MG2","MG3",
  "MG4","MG5","MG6",
  "MG7","MG8","MG9")

datobj_olah_mic_subset <- SetIdent(datobj_olah_mic_subset, value = factor(Idents(datobj_olah_mic_subset), levels = cluster_order))

DimPlot(datobj_olah_mic_subset, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE, group.by = "our_cluster_assignment") + theme(legend.position = "right")

## Dotplot of combined genelist with reprocessed object from 12-08-24
combined_genelist <- c("SPI1", "APOE","CLEC7A","CST7","SPP1",
                       "CSF1R","CX3CR1","TMEM119","P2RY12",
                       "CD28","PDCD1","CD5","DKK2")

DotPlot(datobj_olah_mic_subset, features = combined_genelist, col.min = -1, col.max = 1,
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

##################################
## Check enrichment of genesets ##
##################################
## microglial subset cell type markers from scrnaseq papers
homeostatic_genes <- c("MS4A6A","CSF1R","CD33","CX3CR1","AIF1","TMEM119","P2RY12")

datobj_olah_mic_subset <- datobj_olah_mic_subset

## Consensus homeostatic signature
Idents(datobj_olah_mic_subset) <- "our_cluster_assignment"
datobj_olah_mic_subset <- AddModuleScore(datobj_olah_mic_subset, features = list(intersect(homeostatic_genes, rownames(datobj_olah_mic_subset))), name = "Homeostatic_signature", ctrl = 100, seed = 123456)

FeaturePlot(datobj_olah_mic_subset, features = "Homeostatic_signature1", label = TRUE, repel = TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## Consensus DAM signature
Consensus_DAM_signature <- c("APOE","B2M","CCL3","CTSB","CTSA","CTSD","CTSL","CTSZ","LGALS3","MITF",
                             "TIMP2","TREM2","SPP1","CD9","ITGAX","LPL","LILRB4","CD68","CST7","BHLHE40",
                             "CXCR4","IGF1","HIF1A","ANXA5","CTSS")

datobj_olah_mic_subset <- AddModuleScore(datobj_olah_mic_subset, features = list(intersect(Consensus_DAM_signature, rownames(datobj_olah_mic_subset))), name = "Consensus_DAM_signature", ctrl = 100, seed = 123456)

FeaturePlot(datobj_olah_mic_subset, features = "Consensus_DAM_signature1", label = TRUE, repel = TRUE, order = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## Lymphoid signature
lymphoid_genes <- c("CD28", "PDCD1", "CD5", "DKK2")
datobj_olah_mic_subset <- AddModuleScore(datobj_olah_mic_subset, features = list(intersect(lymphoid_genes5, rownames(datobj_olah_mic_subset))), name = "Lymphoid_signature", ctrl = 100, seed = 123456)

###################
## Module scores ##
###################
## Check module scores using reclustered data not annotation from original paper
dam_score <- datobj_olah_mic_subset@meta.data |> 
  dplyr::select(our_cluster_assignment, Consensus_DAM_signature1) |> 
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_dam = mean(Consensus_DAM_signature1, na.rm = TRUE))

hom_score <- datobj_olah_mic_subset@meta.data |> 
  dplyr::select(our_cluster_assignment, Homeostatic_signature1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_hom = mean(Homeostatic_signature1, na.rm = TRUE))

lymphoid_score_selected <- datobj_olah_mic_subset@meta.data |> 
  dplyr::select(our_cluster_assignment, Lymphoid_signature_selected1) |>
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_score_lymphoid_selected = mean(Lymphoid_signature_selected1, na.rm = TRUE))

pu1_expression <- FetchData(datobj_olah_mic_subset, vars = "SPI1", layer = "data") |> #lognorm counts
  distinct(cell_id, .keep_all = TRUE) |>
  left_join(datobj_olah_mic_subset@meta.data |> rownames_to_column("cell_id")) |> 
  dplyr::select(our_cluster_assignment, SPI1) |> 
  group_by(our_cluster_assignment) |>
  dplyr::summarize(avg_pu1_expr = mean(SPI1, na.rm = TRUE))

module_score_comparison <- dam_score |>
  left_join(hom_score) |>
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
olah_metadata_groups <- datobj_olah_mic_subset@meta.data |>
  mutate(lymphoid_group = case_when(
    our_cluster_assignment %in% c("MG1", "MG2", "MG7") ~ "Homeostatic",
    our_cluster_assignment %in% c("MG6", "MG8", "MG9") ~ "Lymphoid_Positive_DAM",
    our_cluster_assignment %in% c("MG3", "MG4", "MG5") ~ "Lymphoid_Negative_DAM",
    TRUE ~ NA))

unique(olah_metadata_groups$lymphoid_group)

datobj_olah_mic_subset <- AddMetaData(datobj_olah_mic_subset, metadata = olah_metadata_groups$lymphoid_group, col.name = "lymphoid_group")

table(datobj_olah_mic_subset$our_cluster_assignment, datobj_olah_mic_subset$lymphoid_group)

Idents(datobj_olah_mic_subset) <- "lymphoid_group"

## Add SPI1 expression per group
spi1_expr <- FetchData(datobj_olah_mic_subset, vars = "SPI1", slot = "scale.data") |>
  rownames_to_column("cell_id")

olah_metadata_groups <- olah_metadata_groups |>
  rownames_to_column("cell_id") |>
  left_join(spi1_expr)

## UMAP
DimPlot(datobj_olah_mic_subset, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE, group.by = "lymphoid_group") + theme(legend.position = "bottom")

## Reorder clusters
cluster_order <- c(
  "Homeostatic",
  "Lymphoid_Positive_DAM",
  "Lymphoid_Negative_DAM")

datobj_olah_mic_subset <- SetIdent(datobj_olah_mic_subset, value = factor(Idents(datobj_olah_mic_subset), levels = cluster_order))

Idents(datobj_olah_mic_subset) <- "lymphoid_group"

DotPlot(datobj_olah_mic_subset, features = combined_genelist, col.min = -1, col.max = 1,
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
