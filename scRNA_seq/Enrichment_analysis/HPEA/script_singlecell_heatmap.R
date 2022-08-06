
library(Seurat)

mydata_orig <- readRDS("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/RDS_files/seuratObject_rename_ident_May26_2021.rds")
mydata <- subset(mydata_orig, idents = c("Luminal-AV"))

mydata$celltype.stim <- paste(Idents(mydata), mydata$orig.ident, sep = "_")
mydata$celltype <- Idents(mydata)
Idents(mydata) <- "celltype.stim"

genes <- c("Ahr", "Anxa1", "Bmpr2", "Ccnd1", "Cdk4", "Cdk6", "Ctnnb1", "E2f1", "Egfr", "Fzd1", "Gsk3b", "Hipk2", "Odc1", "Pigr", "Plk3", "Ppp4r3a", "Smarca4", "Raf1", "Skp1", "Sp1", "Stat1", "Tp53", "Tcf7l1", "Wnt4")
pdf("Heatmap_breast_cancer_pathway_enriched_genes.pdf", height=15)
DoHeatmap(mydata, features = genes)
dev.off()


