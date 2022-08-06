
library(dplyr)
library(Seurat)
library(cowplot)

mm_3_18 = readRDS("../RDS_files/seurat_black6_scRNA.rds")
cells_per_ident_per_cluster =table(Idents(mm_3_18), mm_3_18@meta.data$orig.ident)
write.table(cells_per_ident_per_cluster, "Table_cells_per_ident_per_cluster.txt", sep="\t", quote=FALSE)

#df_cells_cluster = as.data.frame.matrix(cells_per_cluster_timepoint)
df_cells_cluster = as.data.frame.matrix(cells_per_ident_per_cluster)

Cell_type = rownames(df_cells_cluster)
df_cells_cluster <- cbind(Cell_type, df_cells_cluster)

#Convert ito long format.
library(reshape2)
df_cells_long_format2 = melt(df_cells_cluster, id.vars=c("Cell_type"))
names(df_cells_long_format2) = c("Cell_Type", "Time_point", "prop")

library(ggplot2)
library(stringr)
write.table(df_cells_long_format2,"Table_before_plot.txt", sep="\t")

data_n = read.table("Table_before_plot.txt", header=T, sep="\t")
data_n$Time_point = factor(data_n$Time_point, levels = c("mm10_3mths","mm10_18mths"))

#data_new <- data_n[data_n$Cell_Type != "B/Tcells", ]
data_new <- data_n
#data_new$Time_point = factor(data_new$Time_point, levels = c("M18_rep1","M3_rep1"))
chk_rep = aggregate(data_new$prop, by=list(Time_point=data_new$Time_point), FUN=sum)
chk_rep$Cell_Type = c("ZTotal", "ZTotal")
names(chk_rep) = c("Time_point", "prop", "Cell_Type")
data_new = rbind(data_new, chk_rep)

#data_new$Cell_Type <- factor(data_n_rep1_new$Cell_Type, levels = c("CD4/CD8naive","CD4naive","CD8naive","B cells","Effector_mem", "CD8effector", "NK cells", "Dendritic/Macrophages","Macrophages", "Luminal-AV", "Luminal-HS", "Fibroblasts","Myoepithelial","Vascular/Lymphatic", "ZTotal"))

data_new$Cell_Type <- factor(data_new$Cell_Type, levels = c("2", "3", "13", "14", "15", "1", "6", "0", "4", "17", "9", "11", "7", "8", "10", "16", "5", "12", "ZTotal"))

pdf("Barplot_cell_cluster_per_timepoint_color.pdf", width=25)
ggplot(data_new, aes(fill=Time_point, y=prop, x=Cell_Type)) + geom_bar(position="fill", stat="identity", color="black", size=0.01)+ggtitle("Proportion of cells in cluster by time point \n") + theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16), plot.title=element_text(size=20,face="bold")) + scale_fill_brewer(palette = "Paired") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.text = element_text(size=12)) + geom_hline(yintercept=0.50, linetype="dashed", color = "black")
dev.off()

######################################################################################
cells_per_ident_per_cluster_reps =table(Idents(mm_3_18), mm_3_18@meta.data$replicate)
write.table(cells_per_ident_per_cluster_reps, "Table_cells_per_replicate_per_cluster.txt", sep="\t", quote=FALSE)

#df_cells_cluster = as.data.frame.matrix(cells_per_cluster_timepoint)
df_cells_cluster_reps = as.data.frame.matrix(cells_per_ident_per_cluster_reps)

Cell_type = rownames(df_cells_cluster_reps)
df_cells_cluster_reps <- cbind(Cell_type, df_cells_cluster_reps)
df_cells_cluster_rep1 = df_cells_cluster_reps[,c(1,2,5)]
df_cells_cluster_rep2 = df_cells_cluster_reps[,c(1,3,6)]
df_cells_cluster_rep3 = df_cells_cluster_reps[,c(1,4,7)]

#Convert to long format.
library(reshape2)
df_cells_long_format2_rep1 = melt(df_cells_cluster_rep1, id.vars=c("Cell_type"))
names(df_cells_long_format2_rep1) = c("Cell_Type", "Time_point", "prop")

df_cells_long_format2_rep2 = melt(df_cells_cluster_rep2, id.vars=c("Cell_type"))
names(df_cells_long_format2_rep2) = c("Cell_Type", "Time_point", "prop")

df_cells_long_format2_rep3 = melt(df_cells_cluster_rep3, id.vars=c("Cell_type"))
names(df_cells_long_format2_rep3) = c("Cell_Type", "Time_point", "prop")

library(ggplot2)
library(stringr)
write.table(df_cells_long_format2_rep1,"Table_before_plot_rep1.txt", sep="\t")
write.table(df_cells_long_format2_rep2,"Table_before_plot_rep2.txt", sep="\t")
write.table(df_cells_long_format2_rep3,"Table_before_plot_rep3.txt", sep="\t")

data_n_rep1 = read.table("Table_before_plot_rep1.txt", header=T, sep="\t")
data_n_rep1_new <- data_n_rep1[data_n_rep1$Cell_Type != "B/Tcells",]
data_n_rep1_new$Time_point = factor(data_n_rep1_new$Time_point, levels = c("M18_rep1","M3_rep1"))
chk_rep1 = aggregate(data_n_rep1_new$prop, by=list(Time_point=data_n_rep1_new$Time_point), FUN=sum)
chk_rep1$Cell_Type = c("ZTotal", "ZTotal")
names(chk_rep1) = c("Time_point", "prop", "Cell_Type")
data_n_rep1_new = rbind(data_n_rep1_new, chk_rep1)

data_n_rep1_new$Cell_Type <- factor(data_n_rep1_new$Cell_Type, levels = c("CD4/CD8naive","CD4naive","CD8naive","B cells","Effector_mem", "CD8effector", "NK cells", "Dendritic/Macrophages","Macrophages", "Luminal-AV", "Luminal-HS", "Fibroblasts","Myoepithelial","Vascular/Lymphatic", "ZTotal"))

pdf("Barplot_cell_cluster_per_timepoint_color_rep1.pdf", width=25)
ggplot(data_n_rep1_new, aes(fill=Time_point, y=prop, x=Cell_Type)) + geom_bar(position="fill", stat="identity", color="black", size=0.01)+ggtitle("Proportion of cells in cluster by time point \n") + theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16), plot.title=element_text(size=20,face="bold")) + scale_fill_brewer(palette = "Paired") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.text = element_text(size=12)) + geom_hline(yintercept=0.544, linetype="dashed", color = "black")
dev.off()


