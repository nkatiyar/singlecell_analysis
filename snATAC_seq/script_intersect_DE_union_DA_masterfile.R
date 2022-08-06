
###############################################
#Script to check for overlap between pseudpbulk and single cell DE genes (up and downregulated respectively)
library(Seurat)
library(purrr)
library(dplyr)

get_DE_overlap_results <- function(num){

mydata_all <- names(table(Idents(mydata)))
celltype <- mydata_all[num]
#celltype <- gsub("_", "-", celltype)
celltype <- gsub("-", "_", celltype)
celltype <- gsub("/", "_or_", celltype)

dir_name <- paste0("./", celltype, "_dir_DE_DA")
print(dir_name)

dir.create(dir_name)

celltype_ATAC <- paste0(celltype, "_DA_DESeq2.txt")

#Read the ATACbulk DE genes data...
path_ATAC <- "/scATAC_new/pseudobulk_DA_analysis/DA_analysis/Results/DESeq2_files/"
celltype_ATAC_path <- paste0(path_ATAC,celltype_ATAC)

celltype_ATACbulk_all_log2FC0 <- read.delim(celltype_ATAC_path, sep="\t", header=T)

celltype_ATACbulk_all_log2FC0_filtered <- celltype_ATACbulk_all_log2FC0 %>% filter(M18_M3.logFC > 1)

#Read the single cell DE genes data...
path_DE_union <- "/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/"

Celltype_DE_union_all <- paste0(celltype, "_dir_no_filter/", celltype, "_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt")
celltype_DE_union_file_all <- paste0(path_DE_union, Celltype_DE_union_all)

celltype_DE_union_all_log2FC0 <- read.table(celltype_DE_union_file_all, sep="\t", header=T)
names(celltype_DE_union_all_log2FC0) <- "gene"

############-------------------#########################################

celltype_merged_union_all <- merge(celltype_DE_union_all_log2FC0, celltype_ATACbulk_all_log2FC0_filtered, by = "gene", all=TRUE) 

union_file_all <- paste0(dir_name, "/", celltype, "_DE_union_vs_DA_intersect_pos_log2FC0.txt")

write.table(union_file_all, celltype_merged_union_all, sep="\t", quote=F)

################-----------------######################################

}

###########################################

#file_obj <- readRDS("/scRNAseq/RDS_files/seuratObject_rename_ident_May26_2021.rds")

#############--------------------------------########################
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
{
print("Running overlap analysis...")
}

mydata = readRDS(args[1])
mydata = subset(mydata, idents = c("Doublet", "Plasma"), invert=TRUE)

mydata_all <- names(table(Idents(mydata)))

map(1:length(mydata_all), get_DE_overlap_results)


