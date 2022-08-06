
library(Seurat)
rna_obj <- readRDS("TCells_RNA_subclusters.rds")

pdf("Test_Sid_UMAP_TCells.pdf")
DimPlot(rna_obj, reduction = "humap3")
dev.off()



