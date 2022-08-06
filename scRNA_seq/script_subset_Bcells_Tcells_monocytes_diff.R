
library(dplyr)
library(Seurat)
library(cowplot)
library(harmony)

scRNA_seq_expt <- readRDS("../RDS_files/seurat_black6_scRNA.rds")
scRNA_seq_expt_Bcells <- subset(scRNA_seq_expt, idents = c("2", "3", "13", "17"))
scRNA_seq_expt_Bcells_pos_markers <- FindAllMarkers(scRNA_seq_expt_Bcells, only.pos = TRUE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
scRNA_seq_expt_Bcells_pos_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(scRNA_seq_expt_Bcells_pos_markers, "Table_markers_Bcells_prelim_pos_only.txt", sep="\t", quote=FALSE)

top100 <- scRNA_seq_expt_Bcells_pos_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100 <- top100 %>% filter(p_val_adj <= 0.05)
top100 <- arrange(top100, desc(avg_log2FC), group_by = cluster)
#top100 <- top100[order(-top100$avg_log2FC),]
write.table(top100, "Table_markers_Bcells_pos_only_top100.txt", sep="\t", quote=FALSE)

top50 <- scRNA_seq_expt_Bcells_pos_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50 <- top50 %>% filter(p_val_adj <= 0.05)
top50 <- arrange(top50, desc(avg_log2FC), group_by = cluster)
#top50 <- top50[order(-top50$avg_log2FC),]
write.table(top50, "Table_markers_Bcells_pos_only_top50.txt", sep="\t", quote=FALSE)

top10 <- scRNA_seq_expt_Bcells_pos_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- top10 %>% filter(p_val_adj <= 0.05)
top10 <- arrange(top10, desc(avg_log2FC), group_by = cluster)
#top10 <- top10[order(-top10$avg_log2FC),]
write.table(top10, "Table_markers_Bcells_pos_only_top10.txt", sep="\t", quote=FALSE)

top5 <- scRNA_seq_expt_Bcells_pos_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- top5 %>% filter(p_val_adj <= 0.05)
top5 <- arrange(top5, desc(avg_log2FC), group_by = cluster)
#top5 <- top5[order(-top5$avg_log2FC),]
write.table(top5, "Table_markers_Bcells_pos_only_top5.txt", sep="\t", quote=FALSE)

#############3--------------------------------------######################

scRNA_seq_expt_Tcells <- subset(scRNA_seq_expt, idents = c("0", "1", "4", "6", "15"))
scRNA_seq_expt_Tcells_pos_markers <- FindAllMarkers(scRNA_seq_expt_Tcells, only.pos = TRUE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
scRNA_seq_expt_Tcells_pos_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(scRNA_seq_expt_Tcells_pos_markers, "Table_markers_Tcells_prelim_pos_only.txt", sep="\t", quote=FALSE)

top100 <- scRNA_seq_expt_Tcells_pos_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100 <- top100 %>% filter(p_val_adj <= 0.05)
top100 <- arrange(top100, desc(avg_log2FC), group_by = cluster)
#top100 <- top100[order(-top100$avg_log2FC),]
write.table(top100, "Table_markers_Tcells_pos_only_top100.txt", sep="\t", quote=FALSE)

top50 <- scRNA_seq_expt_Tcells_pos_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50 <- top50 %>% filter(p_val_adj <= 0.05)
top50 <- arrange(top50, desc(avg_log2FC), group_by = cluster)
#top50 <- top50[order(-top50$avg_log2FC),]
write.table(top50, "Table_markers_Tcells_pos_only_top50.txt", sep="\t", quote=FALSE)

top10 <- scRNA_seq_expt_Tcells_pos_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- top10 %>% filter(p_val_adj <= 0.05)
top10 <- arrange(top10, desc(avg_log2FC), group_by = cluster)
#top10 <- top10[order(-top10$avg_log2FC),]
write.table(top10, "Table_markers_Tcells_pos_only_top10.txt", sep="\t", quote=FALSE)

top5 <- scRNA_seq_expt_Tcells_pos_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- top5 %>% filter(p_val_adj <= 0.05)
top5 <- arrange(top5, desc(avg_log2FC), group_by = cluster)
#top5 <- top5[order(-top5$avg_log2FC),]
write.table(top5, "Table_markers_Tcells_pos_only_top5.txt", sep="\t", quote=FALSE)

#############------------------------------------#########################

scRNA_seq_expt_monocytes <- subset(scRNA_seq_expt, idents = c("9", "11"))
scRNA_seq_expt_monocytes_pos_markers <- FindAllMarkers(scRNA_seq_expt_monocytes, only.pos = TRUE, min.pct = 0.25, test.use = 'LR', latent.vars = 'expt', logfc.threshold = 0.25)
scRNA_seq_expt_monocytes_pos_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.table(scRNA_seq_expt_monocytes_pos_markers, "Table_markers_monocytes_prelim_pos_only.txt", sep="\t", quote=FALSE)

top100 <- scRNA_seq_expt_monocytes_pos_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top100 <- top100 %>% filter(p_val_adj <= 0.05)
top100 <- arrange(top100, desc(avg_log2FC), group_by = cluster)
#top100 <- top100[order(-top100$avg_log2FC),]
write.table(top100, "Table_markers_monocytes_pos_only_top100.txt", sep="\t", quote=FALSE)

top50 <- scRNA_seq_expt_monocytes_pos_markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50 <- top50 %>% filter(p_val_adj <= 0.05)
top50 <- arrange(top50, desc(avg_log2FC), group_by = cluster)
#top50 <- top50[order(-top50$avg_log2FC),]
write.table(top50, "Table_markers_monocytes_pos_only_top50.txt", sep="\t", quote=FALSE)

top10 <- scRNA_seq_expt_monocytes_pos_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- top10 %>% filter(p_val_adj <= 0.05)
top10 <- arrange(top10, desc(avg_log2FC), group_by = cluster)
#top10 <- top10[order(-top10$avg_log2FC),]
write.table(top10, "Table_markers_monocytes_pos_only_top10.txt", sep="\t", quote=FALSE)

top5 <- scRNA_seq_expt_monocytes_pos_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5 <- top5 %>% filter(p_val_adj <= 0.05)
top5 <- arrange(top5, desc(avg_log2FC), group_by = cluster)
#top5 <- top5[order(-top5$avg_log2FC),]
write.table(top5, "Table_markers_monocytes_pos_only_top5.txt", sep="\t", quote=FALSE)


