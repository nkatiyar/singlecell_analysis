#Script to analyze single cell RNA-Seq data using Seurat.
library(dplyr)
library(Seurat)
library(cowplot)
library(harmony)

dir.create("../Plots")
dir.create("../Tables")

#mm_3_18_comb_before_filter_batch1.rds
scRNA_seq_expt1 <- readRDS("../RDS_files/mm_3_18_comb_before_filter_batch1.rds")
scRNA_seq_expt2 <- readRDS("../RDS_files/mm_3_18_comb_before_filter_batch2.rds")

print(str(scRNA_seq_expt1@meta.data))
print(str(scRNA_seq_expt2@meta.data))

scRNA_seq_expt1$expt <- sample(c("expt1"), size = ncol(scRNA_seq_expt1), replace = T)
scRNA_seq_expt2$expt <- sample(c("expt2"), size = ncol(scRNA_seq_expt2), replace = T)

###################################################################
scRNA_seq_expt_comb <- merge(scRNA_seq_expt1, y = c(scRNA_seq_expt2), add.cell.ids=c("expt1", "expt2"), project = "Mouse aging cancer")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
scRNA_seq_expt_comb[["percent.mt"]] <- PercentageFeatureSet(scRNA_seq_expt_comb, pattern = "^mt-")

# Visualize QC metrics as a violin plot
pdf("../Plots/Violin_plot_scRNA_seq_RNA_feature.pdf")
VlnPlot(scRNA_seq_expt_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
saveRDS(scRNA_seq_expt_comb, "scRNA_seq_expt_comb_before_filter_May26_2021.rds")

scRNA_seq_expt_comb <- readRDS("scRNA_seq_expt_comb_before_filter_May26_2021.rds")
pdf("../Plots/Plots_qc_orig.pdf", height=10, width=10)
plot1 <- FeatureScatter(scRNA_seq_expt_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA_seq_expt_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#scRNA_seq_expt <- subset(scRNA_seq_expt_comb, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & nCount_RNA < 100000 & percent.mt < 40)
scRNA_seq_expt_comb <- subset(scRNA_seq_expt_comb, subset = nFeature_RNA > 500 & percent.mt < 10)

#Normalizing the data
scRNA_seq_expt <- NormalizeData(scRNA_seq_expt_comb, normalization.method = "LogNormalize", scale.factor = 10000)
#scRNA_seq_expt <- NormalizeData(scRNA_seq_expt)

#Identification of highly variable features (feature selection)
scRNA_seq_expt <- FindVariableFeatures(scRNA_seq_expt, selection.method = "vst", nfeatures = length(rownames(scRNA_seq_expt)))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA_seq_expt), 10)
top5 <- head(VariableFeatures(scRNA_seq_expt), 5)

# plot variable features with and without labels
pdf("../Plots/Plot_variable_feature.pdf", height=10, width=10)
plot1 <- VariableFeaturePlot(scRNA_seq_expt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

saveRDS(scRNA_seq_expt, file="../RDS_files/scRNA_seq_expt_orig_May26_2021.rds")

scRNA_seq_expt <- readRDS("../RDS_files/scRNA_seq_expt_orig_May26_2021.rds")
all.genes <- rownames(scRNA_seq_expt)
scRNA_seq_expt <- ScaleData(scRNA_seq_expt, features = all.genes)

#Perform linear dimensional reduction
scRNA_seq_expt <- RunPCA(scRNA_seq_expt, features = VariableFeatures(object = scRNA_seq_expt))
print(scRNA_seq_expt[["pca"]], dims = 1:5, nfeatures = 5)

pdf("./Plots/Plot_Viz_loading.pdf")
VizDimLoadings(scRNA_seq_expt, dims = 1:2, reduction = "pca")
dev.off()

pdf("./Plots/Plot_DimPlot.pdf")
DimPlot(scRNA_seq_expt, reduction = "pca")
dev.off()

pdf("./Plots/Heatmap.pdf")
DimHeatmap(scRNA_seq_expt, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf("./Plots/Heatmap_dim18.pdf")
DimHeatmap(scRNA_seq_expt, dims = 1:18, cells = 500, balanced = TRUE)
dev.off()

pdf("./Plots/Plot_ElbowPlot.pdf")
ElbowPlot(scRNA_seq_expt)
dev.off()

saveRDS(scRNA_seq_expt, file = "../RDS_files/scRNA_seq_expt_before_clustering_May26_2021.rds")

#Cluster the cells
scRNA_seq_expt <- FindNeighbors(scRNA_seq_expt, dims = 1:10)
scRNA_seq_expt <- FindClusters(scRNA_seq_expt, resolution = 0.5)

#Run non-linear dimensional reduction (UMAP/tSNE)
scRNA_seq_expt <- RunUMAP(scRNA_seq_expt, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("../Plots/Plot_UMAP.pdf")
DimPlot(scRNA_seq_expt, reduction = "umap", label=TRUE)
dev.off()

saveRDS(scRNA_seq_expt, file = "../RDS_files/scRNA_seq_expt_out_May26_2021.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
scRNA_seq_expt <- readRDS("../RDS_files/scRNA_seq_expt_out_May26_2021.rds")
#save(scRNA_seq_expt, file="scRNA_seq_expt_new.RData")
# Visualization
p1 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE)

pdf("../Plots/Dimplot_by_time_and_orig_cluster.pdf", height=10, width=20)
plot_grid(p1,p2)
dev.off()

# Visualization
p3 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "expt", label = TRUE)
p4 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE)

pdf("../Plots/Dimplot_by_expt_and_orig_cluster.pdf", height=10, width=20)
plot_grid(p3,p4)
dev.off()

scRNA_seq_expt <- scRNA_seq_expt %>% RunHarmony("expt", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scRNA_seq_expt, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
pdf("../Plots/DimPlot_after_harmony.pdf")
p1 <- DimPlot(object = scRNA_seq_expt, reduction = "harmony", pt.size = .1, group.by = "expt")
p2 <- VlnPlot(object = scRNA_seq_expt, features = "harmony_1", group.by = "expt", pt.size = .1)
plot_grid(p1,p2)
dev.off()

scRNA_seq_expt <- scRNA_seq_expt %>% 
    RunUMAP(reduction = "harmony", dims = 1:10) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

options(repr.plot.height = 4, repr.plot.width = 6)
pdf("../Plots/Dimplot_umap_after_harmony.pdf")
DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, pt.size = .1, label.size=5)
dev.off()

# Visualization
p3 <- DimPlot(scRNA_seq_expt, reduction = "umap", group.by = "expt", label = TRUE, label.size=5)
p4 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("../Plots/Dimplot_by_expt_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p3,p4)
dev.off()

p5 <- DimPlot(scRNA_seq_expt, reduction = "umap", split.by = "expt", label = TRUE, label.size=5)
p6 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("../Plots/Dimplot_split_expt_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p5,p6)
dev.off()

p7 <- DimPlot(scRNA_seq_expt, reduction = "umap", split.by = "orig.ident", label = TRUE, label.size=5)
p8 <- DimPlot(scRNA_seq_expt, reduction = "umap", label = TRUE, label.size=5)

pdf("../Plots/Dimplot_split_age_and_orig_cluster_after_harmony.pdf", height=10, width=20)
plot_grid(p7,p8)
dev.off()

saveRDS(scRNA_seq_expt, file = "../RDS_files/seurat_black6_scRNA_May26_2021.rds")

################Find all markers #####################
scRNA_seq_expt.markers <- FindAllMarkers(scRNA_seq_expt, only.pos = FALSE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(scRNA_seq_expt.markers, "../Tables/Table_markers_all_prelim_pos_neg.txt", sep="\t", quote=FALSE)

top10 <- scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10, "../Tables/Table_markers_top10.txt", sep="\t", quote=FALSE)

top5 <- scRNA_seq_expt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(top5, "../Tables/Table_markers_top5.txt", sep="\t", quote=FALSE)

pdf("../Plots/Heatmap_top_expressed.pdf", height=20, width=30)
DoHeatmap(scRNA_seq_expt, features = top10$gene) + NoLegend()
dev.off()

table_expt_rep <- table(scRNA_seq_expt$expt, scRNA_seq_expt$replicate)
write.table(table_expt_rep, "../Tables/Table_expt_rep.txt", sep="\t", quote=F)

table_expt_age <- table(scRNA_seq_expt$expt, scRNA_seq_expt$orig.ident)
write.table(table_expt_age, "../Tables/Table_expt_age.txt", sep="\t", quote=F)

saveRDS(scRNA_seq_expt, file = "../RDS_files/seurat_black6_scRNA_all_markers.rds")

######################################
#cluster.averages <- AverageExpression(my_data, return.seurat = TRUE)
#marker_list <- c("Pecam1", "Cdh5", "Eng", "Krt8", "Krt18", "Krt19", "Krt17", "Krt5", "Krt14", "Acta2", "Myl9", "Mylk", "Myh11", "Prlr", "Cited1", "Esr1", "Pgr", "Prom1","Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf", "Wap", "Glycam1", "Olah", "Mfge8", "Trf", "Wfdc18", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Notch3", "Kit", "Aldh1a3", "Cd14", "Aif1", "Itgax", "Fcgr2b", "Itgam", "Blnk", "Cd79a", "Cd3d", "Gzma", "Ncr", "Cd2441", "Cd8a", "Cd8b1", "Cd4", "Cd44", "Sell", "Nkg7", "Gzmh", "Cst7", "Ccl5", "Hla-b","Gnly", "Gzmb")
marker_list <- c("Pecam1", "Cdh5", "Eng", "Krt8", "Krt18", "Krt19", "Krt17", "Krt5", "Krt14", "Acta2", "Myl9", "Mylk", "Myh11", "Prlr", "Cited1", "Esr1", "Pgr", "Prom1","Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf", "Wap", "Glycam1", "Olah", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Notch3", "Rgs5", "Des", "Kit", "Aldh1a3", "Cd14", "Aif1", "Itgax", "Fcgr1", "Adgre1","C1qa", "C1qb", "C1qc", "Fcgr2b", "Cd209a", "Itgam", "Cd24a", "H2-Ab1", "Blnk", "Cd79a", "Cd79b", "Cd3d", "Cd3g", "Cd3e", "Gzma", "Gzmb", "Ncr", "Cd2441", "Cd8a", "Cd8b1", "Il7r", "Cd4", "Cd44", "Sell", "Nkg7", "Gzmh", "Cst7", "Ccl5", "Hla-b","Gnly")

#-------------------------------#
cluster.averages_tp_new <- AverageExpression(scRNA_seq_expt, return.seurat = TRUE)

#cluster.averages <- AverageExpression(my_data, return.seurat = TRUEn)
Clust_Celltype_all <- subset(cluster.averages_tp_new$RNA@scale.data, rownames(cluster.averages_tp_new$RNA@scale.data) %in% marker_list)

library(pheatmap)
pdf("../Plots/Heatmap_all_markers_cluster.pdf", height=15)
pheatmap(Clust_Celltype_all, annotation_legend = TRUE, cluster_cols=TRUE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 8, fontsize_col = 8)
dev.off()

################################################

djamal_marker_list <- c("Ppbp","Pf4", "S100a8","Cd14","Lyz","Fcgr3a","Cd1c","Cst3","Fcer1g","Il7R","Cd3e","Cd3d","Cd8a","Cd8b","Gzma","Gzmb","Gzmh","Gzmk","Nkg7","Xcl1","Ms4a1","Cd79a","Ighm","Ighd","Tnfrsf17","Mzb1","Lilra4","Irf7","Tcf4","Il3ra", "Lef1", "Ccr7")

cluster.averages_tp_new <- AverageExpression(scRNA_seq_expt, return.seurat = TRUE)
Clust_djamal_Celltype_all <- subset(cluster.averages_tp_new$RNA@scale.data, rownames(cluster.averages_tp_new$RNA@scale.data) %in% djamal_marker_list)

Clust_djamal_Celltype_all_immune <- Clust_djamal_Celltype_all[,c(1,2,3,4,9,11,13,14,16)]

pdf("../Plots/Heatmap_djamal1_markers_row_col_cluster.pdf")
pheatmap(Clust_djamal_Celltype_all_immune, annotation_legend = TRUE, cluster_rows=TRUE, cluster_cols=TRUE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 8, fontsize_col = 8)
dev.off()

###############################################

marker_list_master2 <- c("Pecam1", "Cdh5", "Eng", "Krt8", "Krt18", "Krt19", "Krt17", "Krt5", "Krt14", "Acta2", "Myl9", "Mylk", "Myh11", "Prlr", "Cited1", "Esr1", "Pgr", "Prom1","Mfge8", "Trf", "Csn3", "Wfdc18", "Elf5", "Ltf", "Wap", "Glycam1", "Olah", "Col1a1", "Col1a2", "Col3a1", "Fn1", "Notch3", "Rgs5", "Des", "Kit", "Aldh1a3", "Cd14", "Aif1", "Itgax", "Fcgr1", "Adgre1","C1qa", "C1qb", "C1qc", "Fcgr2b", "Cd209a", "Itgam", "Cd24a", "H2-Ab1", "Blnk", "Cd79a", "Cd79b", "Cd3d", "Cd3g", "Cd3e", "Gzma", "Gzmb", "Ncr", "Cd2441", "Cd8a", "Cd8b1", "Il7r", "Cd4", "Cd44", "Sell", "Nkg7", "Gzmh", "Cst7", "Ccl5", "Hla-b","Gnly", "Xcl1","Il2rb","Cd3g","Cd28","Cd3e","Cd3d","Cd7","Lck","Cd27","Cd2","Icos","Cxcr3","Eomes","Ptprc","Fyn","Il7r","Cd247","Id2","Zap70","Ifngr1","Cd8b1","Fasl","Cd8a","Tcf7","Jak1","Il2rg","Thy1","S1pr1","H2-Q7","Sell","Ccr7","Cd40","Cd69","Icosl","Cxcr6","Il18r1","Ccr2","Il2ra","Cd40lg","Tnfrsf25","Ccr6","Cd44","Il2rb","Il17re","Rorc","Ly6e","Il23r","Tnfsf14","Il27ra","Cd84","Hmgb2","Tnfrsf4","Cdca3","Pcna","Ctla4","Xcl1","Cd4","Ighm","Ighd","Cd79a","Igkc","Cd79b","Iglc2","Iglc3","H2-Ob","Foxp3","Pax5","H2-Oa","Cd74","Fcrla","Cd83","H2-Aa","Cd19","Cd5","Btla","H2-Eb1","H2-Ab1","Iglc1","H2-Q6","Tnfrsf13b","Tnfrsf1b","Irf8","Il21r","Trbc1","Smad7","Iglc2","Cd74","Ighd","Iglc3","H2-Aa","H2-Eb1","H2-Ab1","Igkc","Cd83","Ighm","Cd19","Tnfrsf13c","Irf8","Ccr7","Cxcr4","Irf4","Cd38","Sell","Cd69","Icosl","Cd86","Il21r","Cd86","Lyz2","Fcer1g","Cd74","Ighd","Iglc3","H2-Aa","H2-Eb1","H2-Ab1","Igkc,Cd83,Ighm","Cd19","Tnfrsf13c","Irf8","Ccr7","Cxcr4","Irf4","Cd38","Sell","Cd69","Icosl","Cd86","Il21r","Lyz2","Fcer1g","Plbd1","Csf1r","Cd68","Csf2ra","Sirpa","Clec4a3","Ccl4","Cd14","Cd209a","Clec4a1","Clec12a","Fcgr2b","Fcgr3","Clec10a","Csf2rb","Cx3cr1","Tlr2","Flt3","Mrc1","Batf3","Cd40","Il15","H2-Q6","Cd274","Ncf1","H2-Q7","Cd209f","Cd36","Cd209d","Cd163","Cd209g","Fcgrt","Lamp2","Prkcd","Clec4a2","Qk","Mcl1","Clec4b1","Xbp1","Tbx2","Il34","H2-Ob","Cxcr5","S100a9","S100a8","Csf3r","Ccl6","Ccl3","Cd33","H2-Q10","Csf1","Ly6g","Fcgr4","Cd80")

cluster.averages_tp_new <- AverageExpression(scRNA_seq_expt, return.seurat = TRUE)
Clust_master2_Celltype_all <- subset(cluster.averages_tp_new$RNA@scale.data, rownames(cluster.averages_tp_new$RNA@scale.data) %in% marker_list_master2)

pdf("../Plots/Heatmap_master2_markers_row_col_cluster.pdf", width=10, height=15)
pheatmap(Clust_master2_Celltype_all, annotation_legend = TRUE, cluster_rows=TRUE, cluster_cols=TRUE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 6, fontsize_col = 8)
dev.off()

pdf("../Plots/Heatmap_master2_markers_col_cluster.pdf", width=10, height=12)
pheatmap(Clust_master2_Celltype_all, annotation_legend = TRUE, cluster_rows=FALSE, cluster_cols=TRUE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 8, fontsize_col = 8)
dev.off()

pdf("../Plots/Heatmap_master2_markers_row_cluster.pdf", width=10, height=12)
pheatmap(Clust_master2_Celltype_all, annotation_legend = TRUE, cluster_rows=TRUE, cluster_cols=FALSE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 8, fontsize_col = 8)
dev.off()

#################################################
scRNA_seq_expt <- readRDS("../RDS_files/seurat_black6_scRNA_May26_2021.rds")

pdf("../Plots/FeaturePlot_Immune_Ptprc.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ptprc"))
dev.off()

pdf("../Plots/FeaturePlot_Epithelial_Epcam.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Epcam"))
dev.off()

pdf("../Plots/FeaturePlot_Dendritic_Itgax.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Itgax"))
dev.off()

pdf("../Plots/FeaturePlot_Bcells_Blnk.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Blnk"))
dev.off()

pdf("../Plots/FeaturePlot_Bcells_Cd79a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd79a"))
dev.off()

pdf("../Plots/FeaturePlots_Bcells_Cd79b.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd79b"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_naive_S100a4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("S100a4"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_naive_Cd44.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd44"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_naive_Sell.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sell"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_naive_Il7r.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Il7r"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_Cd3d.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd3d"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_Cd3e.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd3e"))
dev.off()

pdf("../Plots/FeaturePlot_Tcells_Cd3g.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd3g"))
dev.off()

pdf("../Plots/FeaturePlot_CD4_CD8_pos_Cd4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd4"))
dev.off()

pdf("../Plots/FeaturePlot_CD4_CD8_pos_Cd8a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd8a"))
dev.off()

pdf("../Plots/FeaturePlot_CD4_CD8_pos_Cd8b1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd8b1"))
dev.off()

pdf("../Plots/FeaturePlot_naive_Bcells_marker.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ighd"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Igha.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Igha"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Ighg1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ighg1"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Ighg3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ighg3"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Cd27.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd27"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Ccr7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ccr7"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Sell.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sell"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Ms4a1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ms4a1"))
dev.off()

pdf("../Plots/FeaturePlot_mem_Bcells_marker_Cd19.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd19"))
dev.off()

pdf("../Plots/FeaturePlot_plasma_Bcells_marker_Mzb1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mzb1"))
dev.off()

pdf("../Plots/FeaturePlot_plasma_Bcells_marker_Jchain.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Jchain"))
dev.off()

pdf("../Plots/FeaturePlot_plasma_Bcells_marker_Gzmb.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzmb"))
dev.off()

pdf("../Plots/FeaturePlot_plasma_Bcells_marker_Cd20.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd20"))
dev.off()

pdf("../Plots/FeaturePlot_plasma_Bcells_marker_Ms4a1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ms4a1"))
dev.off()

pdf("../Plots/FeaturePlot_S100A4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("S100a4"))
dev.off()

pdf("../Plots/FeaturePlot_Cd24.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd24"))
dev.off()

pdf("../Plots/FeaturePlot_Cd38.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd38"))
dev.off()

pdf("../Plots/FeaturePlot_Cxcr5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cxcr5"))
dev.off()

pdf("../Plots/FeaturePlot_Pax5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pax5"))
dev.off()

pdf("../Plots/FeaturePlot_Spi-B.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Spib"))
dev.off()

pdf("../Plots/FeaturePlot_Pou2af1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pou2af1"))
dev.off()

pdf("../Plots/FeaturePlot_Cd40.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd40"))
dev.off()

pdf("../Plots/FeaturePlot_Cd95.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd95"))
dev.off()

pdf("../Plots/FeaturePlot_Fas.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fas"))
dev.off()

pdf("../Plots/FeaturePlot_Prdm1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Prdm1"))
dev.off()

pdf("../Plots/FeaturePlot_Irf4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Irf4"))
dev.off()

pdf("../Plots/FeaturePlot_Xbp1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Xbp1"))
dev.off()

pdf("../Plots/FeaturePlot_E2a_or_Tcf3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Tcf3"))
dev.off()

######-----------------------------#################

pdf("../Plots/FeaturePlot_MacrophageM1_M2_Cd32.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd32"))
dev.off()

pdf("../Plots/FeaturePlot_MacrophageM1_M2_Cd86.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd86"))
dev.off()

pdf("../Plots/FeaturePlot_MacrophageM1_M2_Cd115.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd115"))
dev.off()

pdf("../Plots/FeaturePlot_MacrophageM1_M2_Cd163.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd163"))
dev.off()

pdf("../Plots/FeaturePlot_Bcell_maybe_naive.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd27"))
dev.off()

pdf("../Plots/FeaturePlot_naive_Sell.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sell"))
dev.off()

pdf("../Plots/FeaturePlot_naive_Lef1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Lef1"))
dev.off()

pdf("../Plots/FeaturePlot_naive_Ccr7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ccr7"))
dev.off()

pdf("../Plots/FeaturePlot_some_immune_cells_Pdgfra.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pdgfra"))
dev.off()

pdf("../Plots/FeaturePlot_some_immune_cells_Pdgfrb.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pdgfrb"))
dev.off()

pdf("../Plots/FeaturePlot_some_immune_cells_Fap.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fap"))
dev.off()

pdf("../Plots/FeaturePlot_luminal_Krt19.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt19"))
dev.off()

pdf("../Plots/FeaturePlot_luminal_Krt18.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt18"))
dev.off()

pdf("../Plots/FeaturePlot_luminal_Krt8.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt8"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-HS_Prlr.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Prlr"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-HS_Cited1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cited1"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-HS_Pgr.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pgr"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-HS_Prom1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Prom1"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-HS_Esr1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Esr1"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Mfge8.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mfge8"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Trf.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Trf"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Csn3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Csn3"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Wfdc18.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Wfdc18"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Elf5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Elf5"))
dev.off()

pdf("../Plots/FeaturePlot_luminal-AV_Ltf.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ltf"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Krt17.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt17"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Krt14.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt14"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Krt5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Krt5"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Acta2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Acta2"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Myl9.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Myl9"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Mylk.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mylk"))
dev.off()

pdf("../Plots/FeaturePlot_Myoepithelial_Myh11.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Myh11"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-ECM_producing_Fn1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fn1"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-ECM_producing_Col1a1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Col1a1"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-ECM_producing_Col1a2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Col1a2"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-ECM_producing_Col3a1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Col3a1"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-Contractile_Acta2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Acta2"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-Contractile_Myl9.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Myl9"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-Contractile_Mylk.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mylk"))
dev.off()

pdf("../Plots/FeaturePlot_Fibroblasts-Contractile_Myh11.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Myh11"))
dev.off()

pdf("../Plots/FeaturePlot_Vascular_Pecam1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pecam1"))
dev.off()

pdf("../Plots/FeaturePlot_Vascular_Cdh5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cdh5"))
dev.off()

pdf("../Plots/FeaturePlot_Vascular_Eng.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Eng"))
dev.off()

pdf("../Plots/FeaturePlot_Vascular_endothelial_Sox17.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sox17"))
dev.off()

pdf("../Plots/FeaturePlot_Vascular_endothelial_Sele.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Sele"))
dev.off()

pdf("../Plots/FeaturePlot_endothelial_Mmrn1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mmrn1"))
dev.off()

pdf("../Plots/FeaturePlot_endothelial_Prox1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Prox1"))
dev.off()

pdf("../Plots/FeaturePlot_endothelial_Flt4.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Flt4"))
dev.off()

pdf("../Plots/FeaturePlot_endothelial_Ccl21a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ccl21a"))
dev.off()

pdf("../Plots/FeaturePlot_Pericytes_Rgs5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Rgs5"))
dev.off()

pdf("../Plots/FeaturePlot_Pericytes_Des.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Des"))
dev.off()

pdf("../Plots/FeaturePlot_Pericytes_Notch3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Notch3"))
dev.off()

pdf("../Plots/FeaturePlot_Immune.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ptprc"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Cd14.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd14"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Aif1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Aif1"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Itgax.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Itgax"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Csf1r.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Csf1r"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Fcgr3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fcgr3"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Adgre1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Adgre1"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Ms4a7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ms4a7"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Ma_Mrc1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mrc1"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Ma_Cd209f.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd209f"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Ma_Cd163.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd163"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Mb_Mmp12.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mmp12"))
dev.off()

pdf("../Plots/FeaturePlot_Macrophages_Mb_Mmp13.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mmp13"))
dev.off()

pdf("../Plots/FeaturePlot_Macropahges1_C1qa.pdf")
FeaturePlot(scRNA_seq_expt, features = c("C1qa"))
dev.off()

pdf("../Plots/FeaturePlot_Macropahges1_C1qb.pdf")
FeaturePlot(scRNA_seq_expt, features = c("C1qb"))
dev.off()

pdf("../Plots/FeaturePlot_Macropahges1_C1qc.pdf")
FeaturePlot(scRNA_seq_expt, features = c("C1qc"))
dev.off()

pdf("../Plots/FeaturePlot_Dendritic_Fcgr2b.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fcgr2b"))
dev.off()

pdf("../Plots/FeaturePlot_Dendritic_Cd209a.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd209a"))
dev.off()

pdf("../Plots/FeaturePlot_Dendritic_Itgam.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Itgam"))
dev.off()

pdf("../Plots/FeaturePlot_Dendritic1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ccr7"))
dev.off()

pdf("../Plots/FeaturePlot_NaturalKiller_Ncr.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ncr"))
dev.off()

pdf("../Plots/FeaturePlot_NaturalKiller_Cd2441.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd2441"))
dev.off()

pdf("../Plots/FeaturePlot_NaturalKiller_Cd56.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd56"))
dev.off()

pdf("../Plots/FeaturePlot_NaturalKiller_Nkg7.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Nkg7"))
dev.off()

pdf("../Plots/FeaturePlot_NaturalKiller_Xcl1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Xcl1"))
dev.off()

pdf("../Plots/FeaturePlot_Cytotoxic_Gzma.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzma"))
dev.off()

pdf("../Plots/FeaturePlot_Cytotoxic_Gzmb.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzmb"))
dev.off()

pdf("../Plots/FeaturePlot_Cytotoxic_Gzmh.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzmh"))
dev.off()

pdf("../Plots/FeaturePlot_Cytotoxic_Gzmk.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Gzmk"))
dev.off()

pdf("../Plots/FeaturePlot_HS-AV_Kit.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Kit"))
dev.off()

pdf("../Plots/FeaturePlot_HS-AV_Aldh1a3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Aldh1a3"))
dev.off()

pdf("../Plots/FeaturePlot_HS-AV_Cd14.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd14"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_genes_Wap.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Wap"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_genes_Glycam1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Glycam1"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_genes_Olah.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Olah"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_related_Mfge8.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mfge8"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_related_Trf.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Trf"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_related_Csn3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Csn3"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_related_Wfdc18.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Wfdc18"))
dev.off()

pdf("../Plots/FeaturePlot_Milk_related_Ltf.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ltf"))
dev.off()

pdf("../Plots/FeaturePlot_Cd138.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd138"))
dev.off()

pdf("../Plots/FeaturePlot_Cd43.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd43"))
dev.off()

pdf("../Plots/FeaturePlot_Cd5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd5"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Ighm.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ighm"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_R3g1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("R3g1"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Mic3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Mic3"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Cr2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cr2"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Cd22.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cd22"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Cr1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Cr1"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Fcer2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Fcer2"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Pax5.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Pax5"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Ebf1.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Ebf1"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Tcf3.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Tcf3"))
dev.off()

pdf("../Plots/FeaturePlot_marginal_zoneB_Slc22a2.pdf")
FeaturePlot(scRNA_seq_expt, features = c("Slc22a2"))
dev.off()

####################################################################
#scRNA_seq_expt = readRDS("seurat_black6_scRNA.rds")

new.cluster.ids = c(`0`="Tcells_naive",`1`="Tcells_mem", `2`="Bcells", `3`="Bcells", `4`="Tcells_naive", `5`="Fibroblasts", `6`="Tcells_naive", `7`="Myoepithelial", `8`="Luminal-AV", `9`="Dendritic/Macrophages",`10`="Luminal-HS", `11`="Dendritic/Macrophages", `12`="Vascular", `13` = "Bcells", `14`="Doublet", `15`="Tcells_mem", `16`="Pericytes", `17`="Plasma")

#names(new.cluster.ids) <- levels(scRNA_seq_expt)
scRNA_seq_expt <- RenameIdents(scRNA_seq_expt, new.cluster.ids)
pdf("Annotated_clusters_UMAP.pdf")
DimPlot(scRNA_seq_expt, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#ff99cc", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

scRNA_seq_expt_new <- subset(scRNA_seq_expt, idents = c("Doublet"), invert = TRUE)
pdf("Annotated_clusters_doublet_rmvd_UMAP_new.pdf")
DimPlot(scRNA_seq_expt_new, cols = c("#FFCC33", "#FF9933", "#746CB1", "#1BBDC1", "#09713A", "#3DB54A", "#c00000", "#A0D082", "#8AB6E1", "#262262", "#BD77B2"), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

#################---------------------------########################

#scRNA_seq_expt = readRDS("seurat_black6_scRNA.rds")
Idents(scRNA_seq_expt_new) <- scRNA_seq_expt@meta.data$orig.ident

scRNA_seq_expt_new$orig.ident <- factor(x = scRNA_seq_expt_new$orig.ident, levels = c("mm10_3mths", "mm10_18mths"))
p9 <- DimPlot(scRNA_seq_expt_new, cols = c("blue", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE)
p10 <- DimPlot(scRNA_seq_expt_new, cols = c("blue", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size=2)

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_new2.pdf", height=10, width=20)
plot_grid(p9)
dev.off()

pdf("Dimplot_split_age_and_orig_cluster_after_harmony_ptsize2.pdf", height=10, width=20)
plot_grid(p10)
dev.off()

################------------------------##############################
# Save a single object to a file
saveRDS(scRNA_seq_expt, "seuratObject_rename_ident_May26_2021.rds")

#################-------------------------############################
#######---------------------------------------------------------################

mydata_new = readRDS("seuratObject_rename_ident.rds")
mydata_new_rmv_doublet <- subset(mydata_new, idents = c("Doublet"), invert = TRUE)
#mydata_new_rmv_doublet <- subset(mydata_new_rmv_doublet, idents = c("Dend/Macro/Vascular"), invert = TRUE)
scRNA_seq_expt.markers_new <- FindAllMarkers(mydata_new_rmv_doublet, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
scRNA_seq_expt.markers_new %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

write.table(scRNA_seq_expt.markers_new, "Table_markers_all_pos_neg_rename_clusters.txt", sep="\t", quote=FALSE)
top10_new <- scRNA_seq_expt.markers_new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

pdf("Heatmap_top_expressed_rename_clusters_new.pdf", height=15, width=20)
DoHeatmap(mydata_new, features = top10_new$gene, size=3)
dev.off()

#.....................#
######################################
Clust_Celltype_clusters_orig <- subset(cluster.averages_tp_new_orig$RNA@scale.data, rownames(cluster.averages_tp_new_orig$RNA@scale.data) %in% marker_list)

cluster.averages_tp_new <- AverageExpression(mydata_new, return.seurat = TRUE)

Clust_Celltype_all <- subset(cluster.averages_tp_new$RNA@scale.data, rownames(cluster.averages_tp_new$RNA@scale.data) %in% top10_new$gene)

library(pheatmap)
pdf("Heatmap_all_markers_cluster.pdf")
pheatmap(Clust_Celltype_all, annotation_legend = TRUE, cluster_cols=FALSE, cellwidth = 20,treeheight_row = 0, treeheight_col = 0, fontsize_row = 4, fontsize_col = 8)
dev.off()

########################################################
############--------------------#####################

