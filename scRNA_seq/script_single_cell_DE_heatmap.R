
########################################################
library(Seurat)
library(purrr)
library(dplyr)
dir.create("DE_analysis")

get_scDE_results <- function(num, A, B, x, y){
contrast = "18_vs_3"

celltype = paste0(mydata_all[num], "_", contrast)
celltype <- gsub("/", "_", celltype)

ident1 = paste0(mydata_all[num], "_", A)
ident2 = paste0(mydata_all[num], "_", B)

print(paste0("DE analysis for ", mydata_all[num], "...."))
var1 <- paste0(celltype, "_log1.markers")
mydata = y
#var1.markers <- FindMarkers(mydata, ident.1 = ident1, ident.2 = ident2, test.use = 'LR', latent.vars = 'expt', verbose = FALSE, logfc.threshold = 0)

#var1_file <- paste0("./DE_analysis/Diff_",var1, ".txt")
#write.table(var1.markers, var1_file, sep="\t", quote=FALSE)

#celltype_log1.markers <- read.table("./DE_analysis/Diff_T_or_NKcells_18_vs_3mths_log1_LR.txt", sep="\t", header=T)
#var1.markers$gene = rownames(var1.markers)
#celltype_DE_log0_pos = var1.markers %>% filter(avg_log2FC > 0, p_val_adj <= 0.05)
#celltype_DE_log0_neg = var1.markers %>% filter(avg_log2FC < 0, p_val_adj <= 0.05)
#celltype_DE_log0_neg <- celltype_DE_log0_neg[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
#celltype_DE_log0_pos <- celltype_DE_log0_pos[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

file1_log0_neg <- paste0("../DE_analysis/", celltype, "_DE_neg_log0_LR.txt")
file1_log0_pos <- paste0("../DE_analysis/", celltype, "_DE_pos_log0_LR.txt")

celltype_DE_log0_neg <- read.table(file1_log0_neg, sep="\t", header=T)
celltype_DE_log0_pos <- read.table(file1_log0_pos, sep="\t", header=T)

write.table(celltype_DE_log0_neg, file1_log0_neg, sep="\t", quote=FALSE)
write.table(celltype_DE_log0_pos, file1_log0_pos, sep="\t", quote=FALSE)

DE_genes_heatmap <- c(celltype_DE_log0_neg$gene, celltype_DE_log0_pos$gene)
file_heatmap <- paste0("Heatmap_single_cell_", celltype, ".pdf")

mydata_subset <- subset(mydata, idents = c(ident1, ident2))

pdf(file_heatmap, height=80)
DoHeatmap(mydata_subset, features = DE_genes_heatmap, size=2) + NoLegend()
dev.off()

}

#############--------------------------------########################
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
{
print("Running differential analysis...")
}

mydata = readRDS(args[1])
#mydata = readRDS("/projects/anczukow-lab/neerja_projects/Integrate_expt1_expt2_only/Analysis_PC_10_mt10/RDS_files/seuratObject_rename_ident.rds")
mydata$celltype.stim <- paste(Idents(mydata), mydata$orig.ident, sep = "_")
mydata$celltype <- Idents(mydata)
Idents(mydata) <- "celltype.stim"
saveRDS(mydata, "seuratObject_rename_ident_DE_tp.rds")

mydata_all <- names(table(mydata$celltype))

map(1:length(mydata_all), get_scDE_results, A= "mm10_18mths", B="mm10_3mths", x = mydata_all, y = mydata)


