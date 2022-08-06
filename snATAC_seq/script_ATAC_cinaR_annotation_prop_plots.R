
###############################################
#Script to check for overlap between pseudpbulk and single cell DE genes (up and downregulated respectively)
library(Seurat)
library(purrr)
library(dplyr)
library(ggplot2)

dir.create("ATAC_barplots")

get_DE_overlap_results <- function(num){

mydata_all <- names(table(Idents(mydata)))
celltype <- mydata_all[num]
#celltype <- gsub("_", "-", celltype)
celltype <- gsub("-", "_", celltype)
celltype <- gsub("/", "_or_", celltype)

#Read the ATACbulk DE genes data...
path_ATAC <- "/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scATAC_new/pseudobulk_DA_analysis/DA_analysis/Results/DESeq2_files/"
celltype_ATAC <- paste0(celltype, "_DA_DESeq2.txt")

celltype_ATAC_path <- paste0(path_ATAC,celltype_ATAC)

celltype_ATACbulk_all_log2FC0 <- read.delim(celltype_ATAC_path, sep="\t", header=T)

celltype_df <- celltype_ATACbulk_all_log2FC0[c(1,7)]

celltype_df$M18_M3.annotation <- gsub(".*Exon.*", "Exon", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*Intron.*", "Intron", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*Promoter.*", "Promoter", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*3' UTR.*", "3' UTR", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*5' UTR.*", "5' UTR", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*Distal Intergenic.*", "Distal Intergenic", celltype_df$M18_M3.annotation)
celltype_df$M18_M3.annotation <- gsub(".*Downstream.*", "Downstream", celltype_df$M18_M3.annotation)

celltype_df$celltype <- celltype
res <- celltype_df %>% group_by(celltype,M18_M3.annotation) %>% summarise(Freq=n())

df_comb <- rbind(df_comb, res)

#For DA peaks - opening, closing and total.

celltype_df2 <- celltype_ATACbulk_all_log2FC0[c(1,7,17)]
celltype_df2$M18_M3.annotation <- gsub(".*Exon.*", "Exon", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*Intron.*", "Intron", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*Promoter.*", "Promoter", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*3' UTR.*", "3' UTR", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*5' UTR.*", "5' UTR", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*Distal Intergenic.*", "Distal Intergenic", celltype_df2$M18_M3.annotation)
celltype_df2$M18_M3.annotation <- gsub(".*Downstream.*", "Downstream", celltype_df2$M18_M3.annotation)

celltype_df2_filter <- celltype_df2 %>% filter(abs(M18_M3.logFC) > 1) # Filter DA peaks total

celltype_ATACbulk_open <- celltype_df2_filter %>% filter(M18_M3.logFC > 1) # Filter opening peaks only
celltype_ATACbulk_close <- celltype_df2_filter %>% filter(M18_M3.logFC < 1) # filter closing peaks only

celltype_df_DA <- celltype_df2_filter [c(1,2)]
celltype_df_open <- celltype_ATACbulk_open [c(1,2)]
celltype_df_close <- celltype_ATACbulk_close [c(1,2)]

celltype_df_DA$celltype <- celltype
res_DA <- celltype_df_DA %>% group_by(celltype,M18_M3.annotation) %>% summarise(Freq=n())
res_DA$type <- "All_DA_peaks"

celltype_df_open$celltype <- celltype
res_open <- celltype_df_open %>% group_by(celltype,M18_M3.annotation) %>% summarise(Freq=n())
res_open$type <- "Opening"

celltype_df_close$celltype <- celltype
res_close <- celltype_df_close %>% group_by(celltype,M18_M3.annotation) %>% summarise(Freq=n())
res_close$type <- "Closing"

res_comb1 <- rbind(res_DA, res_open)
res_comb2 <- rbind(res_comb1, res_close)

res_comb2$type <- factor(res_comb2$type, levels = c("Opening", "Closing", "All_DA_peaks"))
filename_comb <- paste0("./ATAC_barplots/", celltype, "_DA_annotation_combined.pdf")
pdf(filename_comb)
print(ggplot(res_comb2, aes(fill = M18_M3.annotation, y = Freq, x = type)) + geom_bar(position="stack", stat="identity", width=0.5) + scale_fill_manual(values = c("#63C29A", "#C15E44", "#F0C144", "#4CB8D8", "#AC67E9", "#F089E7", "#F1E6A1")))
dev.off()

#Barplots for immune cells with labels (order by age).
filename_plot <- paste0("./ATAC_barplots/", celltype, "_barplot_celltype_annotation_prop.pdf")
pdf(filename_plot)
print(ggplot(res_comb2, aes(fill=M18_M3.annotation, y=Freq, x=type)) + geom_bar(position="fill", stat="identity", colour="black", width=0.5) + scale_fill_manual(values=c("#63C29A", "#C15E44", "#F0C144", "#4CBBD8", "#AC67E9", "#F089E7", "#F1E6A1")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size=6)) + scale_y_continuous(expand = c(0,0)))
dev.off()

filename_plot <- paste0("./ATAC_barplots/", celltype, "_barplot_celltype_annotation_prop_label.pdf")
pdf(filename_plot)
print(ggplot(res_comb2, aes(fill=M18_M3.annotation, y=Freq, x=type)) + geom_bar(position="fill", stat="identity", colour="black", width=0.5) + scale_fill_manual(values=c("#63C29A", "#C15E44", "#F0C144", "#4CBBD8", "#AC67E9", "#F089E7", "#F1E6A1")) + geom_text(aes(label=Freq), position = position_fill(vjust=0.5)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size=6)) + scale_y_continuous(expand = c(0,0)))
dev.off()


################-----------------######################################
return(df_comb)

}

###########################################
#mydata <- readRDS("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/RDS_files/seuratObject_rename_ident_May26_2021.rds")

#############--------------------------------########################
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
{
print("Running overlap analysis...")
}

df_comb <- data.frame()

mydata = readRDS(args[1])
mydata = subset(mydata, idents = c("Doublet", "Plasma"), invert=TRUE)

mydata_all <- names(table(Idents(mydata)))

k <- map(1:length(mydata_all), get_DE_overlap_results)

df <- do.call(Map, c(f=c, k))
k_df <- as.data.frame(df)

k_df$celltype <- factor(k_df$celltype, levels = c("Tcells_naive", "Tcells_mem", "Bcells", "Dendritic_or_Macrophages", "Luminal_AV", "Luminal_HS", "Myoepithelial", "Fibroblasts", "Vascular", "Pericytes"))

pdf("./ATAC_barplots/All_celltypes_annotation_summary.pdf", width=20, height=10)
print(ggplot(k_df, aes(fill = M18_M3.annotation, y = Freq, x = celltype)) + geom_bar(position="stack", stat="identity", width=0.5) + scale_fill_manual(values = c("#63C29A", "#C15E44", "#F0C144", "#4CB8D8", "#AC67E9", "#F089E7", "#F1E6A1")) + theme(text = element_text(size=15), axis.text.x = element_text(angle=45, hjust=1)))
dev.off()

