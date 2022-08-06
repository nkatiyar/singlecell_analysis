
########################################################
library(Seurat)
library(purrr)
library(dplyr)
library(cinaR)
library(qusage)

dir.create("../DE_enrichment")

get_scDE_results <- function(num, x, y){

#celltype = paste0(mydata_all[num], "_", contrast)
celltype <- gsub("/", "_or_", celltype)

print(paste0("DE enrichment analysis for ", mydata_all[num], "...."))
mydata = y

celltype_path <- "/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/"
geneset_path <- "/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/"

var1 <- paste0(celltype_path, "/", celltype, "_mm10_18mths_vs_mm10_3mths_all_genes.csv")

file1_c2 <- paste0(geneset_path, "c2.cp.kegg.v7.4.symbols.gmt")
gmtfile <- read.gmt(file1_c2)

celltype_DE_log0 <- read.table(var1, sep=",", header=T)

#celltype_DE_log0_neg <- read.table(file1_log0_neg, sep="\t", header=T)
#celltype_DE_log0_pos <- read.table(file1_log0_pos, sep="\t", header=T)

celltype_genes_human <- mouse2human(celltype_DE_log0$gene)
celltype_DE_log0$gene <- celltype_genes_human$HGNC.symbol

df_inp <- celltype_DE_log0[,c(1,3)] 
df_inp_sorted <- df_inp[order(df_inp[,2]),]

res_gsea <- GSEA(df_inp_sorted, gmtfile)

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

map(1:length(mydata_all), get_scDE_results, x = mydata_all, y = mydata)


