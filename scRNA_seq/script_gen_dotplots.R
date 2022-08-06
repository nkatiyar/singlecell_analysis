
library(Seurat)
library(ggplot2)

#Script to generate DotPlots for marker genes.
#mydata_new = readRDS("seuratObject_rename_ident.rds")

dir.create("../DotPlots")
gen_dotplots <- function(mydata_new){
	mydata_new_rmv_doublet <- subset(mydata_new, idents = c("Doublet"), invert = TRUE)
	marker_list_orig <- c("Pecam1", "Eng", "Cdh5", "Rgs5", "Notch3", "Des", "Col1a1", "Fn1", "Ighg1", "Ighg3", "Cd27", "Jchain", "Mzb1", "Gzmb", "Ms4a2", "Ighd", "Blnk","Cd79a", "Cd8a", "Cd8b1", "Cd3d", "Cd4", "Sell", "Cd44", "Itgam", "Csf1r", "Cd86", "Cd163", "Ly75", "Fcgr2b", "Itgax", "Aif1","Myl9", "Acta2", "Krt17", "Prlr", "Cited1", "Esr1",  "Mfge8", "Trf", "Csn3", "Pi16", "Clec4a", "Cd38", "IgG", "Tnfrsf13b", "Tnfrsf17", "Il7r", "Ccr7", "Cd19", "Ms4a1", "Cd34", "B3gat1", "Klrg1", "Eomes", "Prdm1", "Ptprc", "Id2", "Stat4", "T-bet", "Tcf1", "Bcl6", "Id3", "Stat3", "Cx3cr1")
	
	marker_list_subset <- c("Pecam1", "Eng", "Cdh5", "Rgs5", "Notch3", "Des", "Col1a1", "Fn1", "Blnk", "Jchain", "Mzb1","Cd79a", "Cd8a", "Cd8b1", "Cd3d", "Cd4", "Sell", "Cd44", "Il7r", "Itgam", "Csf1r", "Cd86", "Cd163", "Fcgr2b", "Itgax", "Aif1","Myl9", "Acta2", "Krt17", "Prlr", "Cited1", "Esr1",  "Mfge8", "Trf", "Csn3")

	#mydata_new_rmv_doublet@active.ident <- factor(mydata_new_rmv_doublet@active.ident,
        #	                    levels=c("Vascular", "Pericytes", "Fibroblasts", "Plasma", "Bcells", "Tcells_mem", "Tcells_naive", "Dendritic/Macrophages","Myoepithelial", "Luminal-HS", "Luminal-AV"))
	mydata_new_rmv_doublet@active.ident <- factor(mydata_new_rmv_doublet@active.ident,
        	                    levels=c("Luminal-AV", "Luminal-HS", "Myoepithelial", "Dendritic/Macrophages", "Tcells_naive", "Tcells_mem", "Bcells", "Plasma", "Fibroblasts", "Pericytes", "Vascular"))
	#mydata_new_rmv_doublet@active.ident <- factor(mydata_new_rmv_doublet@active.ident,
         #                           levels=c("2", "3", "5", "1", "13", "14", "0", "4", "15", "17", "18", "9", "11", "12", "16", "6", "7", "10", "8"))

        #marker_list_final <- rev(marker_list_subset)

	pdf("./DotPlots/DotPlot_markers_all_clusters.pdf", width=25, height=12)
	print(DotPlot(object = mydata_new_rmv_doublet, features = marker_list_final, cols = c("lightgrey", "red"), dot.scale=10))
	dev.off()
	
	pdf("./DotPlots/DotPlot_markers_all_clusters_flip.pdf", width=15, height=12)
	print(DotPlot(object = mydata_new_rmv_doublet, features = marker_list, cols = c("lightgrey", "red"), dot.scale=10) + coord_flip() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)))  
	dev.off()
	
	pdf("../DotPlots/DotPlot_markers_all_celltypes_flip_new.pdf", width=12, height=18)
	print(DotPlot(object = mydata_new_rmv_doublet, features = marker_list_subset, cols = c("lightgrey", "red"), dot.scale=10) + coord_flip() + theme(legend.position = "bottom", legend.box = "vertical", legend.key.width = unit(1, "cm")) + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) )
	dev.off()

	#######################################
	canonical_epithelial_markers <- c("Krt17", "Acta2", "Myl9", "Prlr", "Cited1", "Esr1", "Mfge8", "Trf", "Csn3", "Wap")
	
	mydata_new_epithelial <- subset(mydata_new, idents = c("Luminal-HS", "Luminal-AV", "Myoepithelial"))
	mydata_new_epithelial@active.ident <- factor(mydata_new_epithelial@active.ident,
        	                    levels=c("Myoepithelial", "Luminal-HS", "Luminal-AV"))
	
	pdf("./DotPlots/DotPlot_epithelial_subset_markers.pdf", width=8, height=5)
	print(DotPlot(object = mydata_new_epithelial, features = rev(canonical_epithelial_markers), cols = c("lightgrey", "red")), dot.scale=10)
	dev.off()
	
	canonical_epithelial_markers_v2 <- c("Krt17", "Krt14", "Krt5", "Acta2", "Myl9", "Mylk", "Myh11", "Krt19", "Krt18", "Krt8", "Prlr", "Cited1", "Pgr", "Esr1", "Mfge8", "Trf", "Csn3", "Elf5", "Ltf", "Kit", "Aldh1a3","Wap", "Glycam1", "Olah")
	
	mydata_new_epithelial_v2 <- subset(mydata_new, idents = c("Luminal-HS", "Luminal-AV", "Myoepithelial"))
	mydata_new_epithelial_v2@active.ident <- factor(mydata_new_epithelial_v2@active.ident,
	                            levels=c("Luminal-AV", "Luminal-HS", "Myoepithelial"))
	
	
	pdf("./DotPlots/DotPlot_epithelial_markers.pdf", width=20, height=5)
	print(DotPlot(object = mydata_new_epithelial_v2, features = canonical_epithelial_markers_v2, cols = c("lightgrey", "red")), dot.scale=10)
	dev.off()
	
	######################################3
	canonical_stromal_markers <- c("Pecam1", "Cdh5", "Eng", "Col1a1", "Fn1", "Acta2", "Myl9", "Notch3", "Cd14", "Aif1", "Itgax", "Fcgr2b", "Cd209a", "Itgam", "Cd3d", "Gzma", "Ncr", "Cd2441", "Cd8a", "Cd8b1", "Blnk", "Cd79a", "Cd79b", "Cd4", "Cd44", "Sell", "Nkg7", "Gzmh", "Cst7", "Ccl5", "Hla-b","Gnly", "Gzmb")
	
	canonical_stromal_markers_final <- rev(canonical_stromal_markers)
	mydata_new_stromal <- subset(mydata_new, idents = c("T/NKcells", "Bcells", "Dendritic/Macrophages", "Fibroblasts", "Vascular"))
	mydata_new_stromal@active.ident <- factor(mydata_new_stromal@active.ident,
        	                    levels=c("T/NKcells", "Bcells", "Dendritic/Macrophages", "Fibroblasts","Vascular"))
	pdf("./DotPlots/DotPlot_stromal_markers.pdf", width=40)
	print(DotPlot(object = mydata_new_stromal, features = canonical_stromal_markers_final, cols = c("lightgrey", "red")), dot.scale=10)
	dev.off()
	
	#####################################3
	canonical_stromal_subset_markers <- c("Pecam1", "Eng", "Cdh5", "Acta2", "Myl9", "Notch3", "Rgs5", "Des", "Pdgfra", "Pdgfrb", "Fap", "Col1a1", "Fn1")
	canonical_stromal_markers_final <- rev(canonical_stromal_subset_markers)
	mydata_new_stromal <- subset(mydata_new, idents = c("Pericytes", "Vascular", "Fibroblasts"))
	mydata_new_stromal@active.ident <- factor(mydata_new_stromal@active.ident,
        	                    levels=c("Vascular", "Pericytes", "Fibroblasts"))
	pdf("./DotPlots/DotPlot_stromal_subset_markers.pdf", height=5, width=15)
	print(DotPlot(object = mydata_new_stromal, features = canonical_stromal_markers_final, cols = c("lightgrey", "red")), dot.scale=10)
	dev.off()
	
	######################################
	canonical_immune_subset_markers <- c("Fcgr2b", "Itgam", "Cd14", "Aif1", "Itgax", "Cd3d", "Cd4", "Cd8a", "Blnk", "Cd79a", "Cd79b", "Jchain")
	canonical_immune_markers_final <- rev(canonical_immune_subset_markers)
	mydata_new_immune <- subset(mydata_new, idents = c("Dendritic/Macrophages", "T/NKcells", "Bcells"))
	mydata_new_immune@active.ident <- factor(mydata_new_immune@active.ident,
	                            levels=c("Bcells", "T/NKcells", "Dendritic/Macrophages"))
	pdf("./DotPlots/DotPlot_immune_subset_markers.pdf", height=5, width=15)
	print(DotPlot(object = mydata_new_immune, features = rev(canonical_immune_markers_final), cols = c("lightgrey", "red")), dot.scale=10)
	dev.off()
	
}

####################################
dir.create("DotPlots")
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else
#paste("seuratObject_rename_ident.rds")
{
print("Generating dotplots using marker genes...")
}

mydata_new = readRDS(args[1])

gen_dotplots(mydata_new)

