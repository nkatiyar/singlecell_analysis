#Script to analyze single cell RNA-Seq data using Seurat.
library(dplyr)
library(Seurat)
library(cowplot)
library(harmony)

scRNA_seq_expt <- readRDS("../RDS_files/seurat_black6_scRNA_May26_2021.rds")

new.cluster.ids = c(`0`="Tcells_naive",`1`="Tcells_mem", `2`="Bcells", `3`="Bcells", `4`="Tcells_naive", `5`="Fibroblasts", `6`="Tcells_naive", `7`="Myoepithelial", `8`="Luminal", `9`="Dendritic/Macrophages",`10`="Luminal", `11`="Dendritic/Macrophages", `12`="Vascular", `13` = "Bcells", `14`="Doublet", `15`="Tcells_mem", `16`="Pericytes", `17`="Plasma")

#names(new.cluster.ids) <- levels(scRNA_seq_expt)
scRNA_seq_expt <- RenameIdents(scRNA_seq_expt, new.cluster.ids)
pdf("Annotated_clusters_lumAV_comb_UMAP.pdf")
DimPlot(scRNA_seq_expt, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#8AB6E1", "#ff99cc", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#DimPlot(scRNA_seq_expt, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#ff99cc", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

scRNA_seq_expt_new <- subset(scRNA_seq_expt, idents = c("Doublet"), invert = TRUE)
pdf("Annotated_clusters_doublet_rmvd_UMAP_new.pdf")
DimPlot(scRNA_seq_expt_new, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#8AB6E1", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

################------------------------##############################
# Save a single object to a file
saveRDS(scRNA_seq_expt, "seuratObject_rename_lum_comb_Mar14_2022.rds")

#################-------------------------############################
#######---------------------------------------------------------################

