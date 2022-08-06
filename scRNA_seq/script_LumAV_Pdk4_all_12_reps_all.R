
#Script to generate expression plots using seurat object.
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(gridBase)
library(gridExtra)
library(tidyr)
library(cowplot)


gene_Vln_Pie <- function(mydata, gene_name, Celltype){
	subset_celltype <- subset(mydata, idents = Celltype)
	subset_celltype$replicate.stim <- paste(Idents(subset_celltype), subset_celltype$replicate, sep = "_")
	subset_celltype$celltype <- Idents(subset_celltype)
	Idents(subset_celltype) <- "replicate.stim"

	df = as.data.frame(subset_celltype[["RNA"]]@data)
	new_df <- df[, df[gene_name,] > 0]
	subset_celltype_gene <- CreateSeuratObject(counts = new_df, assay = "RNA")
	subset_celltype_gene <- AddMetaData(subset_celltype_gene, subset_celltype@meta.data[colnames(subset_celltype_gene), , drop=F])

	exprn_val <- subset(subset_celltype_gene[["RNA"]]@data, rownames(subset_celltype_gene[["RNA"]]@data) %in% gene_name)
	exprn_val_df <- as.data.frame(exprn_val)
	rownames(exprn_val_df) <- colnames(subset_celltype_gene[["RNA"]]@data)
	celltype_age <- subset_celltype_gene@meta.data$replicate.stim
	celltype <- subset_celltype_gene@meta.data$celltype
	df_final_AV <- cbind(exprn_val_df, celltype_age, celltype)
	names(df_final_AV) <- c("exprn", "celltype_age", "celltype")

	#######---------------------############

	a <- as.data.frame(table(subset_celltype_gene$replicate.stim))
	b <- as.data.frame(table(subset_celltype$replicate.stim))
	total = merge(b,a, by="Var1", all=TRUE)
	total[is.na(total)] <- 0
	names(total) <- c("Sample", "Total", "Expressed")
	total$Not_expressed <- total$Total - total$Expressed

	cells_pie <- total[,c("Expressed", "Not_expressed")]
	cells_pie$Celltype <- total$Sample
	cells_pie <- cells_pie[c("Celltype", "Expressed", "Not_expressed")]
	cells_pie_long <- gather(cells_pie, condition, cells, Expressed:Not_expressed, factor_key=TRUE)

	level_order <- c(paste0(Celltype,"_M3_rep1"), paste0(Celltype,"_M3_rep2"), paste0(Celltype,"_M3_rep3"), paste0(Celltype,"_M3_rep4"), paste0(Celltype,"_M3_rep5"), paste0(Celltype,"_M3_rep6"), paste0(Celltype,"_M18_rep1"), paste0(Celltype,"_M18_rep2"), paste0(Celltype,"_M18_rep3"), paste0(Celltype,"_M18_rep4"), paste0(Celltype,"_M18_rep5"), paste0(Celltype,"_M18_rep6"))

	#########Violin Plot ######################
	vlnplot_n <- df_final_AV %>%
  	ggplot(aes(y=exprn, x = factor(celltype_age, level = level_order), fill = factor(celltype_age, level = level_order))) +
  	geom_violin(position = position_dodge(0.8), alpha=0.5) +
    	scale_fill_manual(values=c("springgreen", "springgreen", "springgreen", "springgreen", "springgreen", "springgreen", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue")) + geom_jitter(shape=16, position=position_jitter(0.2),cex=5) +
    	xlab("Sample")+ylab("Expression")+
    	theme(plot.margin = unit(c(4,4,20,10), "cm"))+ggtitle(gene_name)+theme(plot.title = element_text(hjust = 0.5, size=80))+ theme(axis.text.x=element_text(size=60, angle=90), axis.text.y = element_text(size = 60), axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60), legend.title=element_text(size=60), legend.text=element_text(size=60))+theme(legend.key.size = unit(2, "cm"))

	###### Pie chart #########################
	cells_pie_long$Celltype = factor(cells_pie_long$Celltype, levels=level_order)
	pie_all <- ggplot(data = cells_pie_long, aes(x="", y=cells, fill=condition)) + geom_bar(stat = "identity", position = "fill") + theme_void() +
        coord_polar(theta = "y", start=0) + facet_wrap(~Celltype, nrow=1) + scale_fill_grey() + theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'))

	pie_only <- ggplot(data = cells_pie_long, aes(x="", y=cells, fill=condition)) + geom_bar(stat = "identity", position = "fill") + theme_void() + 
	coord_polar(theta = "y", start=0) + facet_wrap(~Celltype, nrow=1) + scale_fill_grey() + theme(legend.position = "right", legend.title = element_text(size = 10), legend.text = element_text(size = 10), legend.key.size = unit(1, 'cm'))

	Celltype = gsub("/", "_or_", Celltype)
	outfile = paste0("./Violin_plots/","VlnPlot_Pie_chart_",Celltype, "_", gene_name, ".pdf")
	outfile1 = paste0("./Violin_plots/","VlnPlot_",Celltype, "_", gene_name, ".pdf")
	outfile2 = paste0("./Violin_plots/","Pie_chart_",Celltype, "_", gene_name, ".pdf")

#pdf("VlnPlot_Pie_chart_comb.pdf", width=120, height=40)
	pdf(file = outfile, width=120, height=40)
	print(ggdraw()+draw_plot(vlnplot_n, x=0, y=0, width=1, height=1)+draw_plot(pie_all, x=0.05, y=-0.4, width=0.7))
	dev.off()

	pdf(file = outfile1, width=120, height=40)
	print(ggdraw()+draw_plot(vlnplot_n, x=0, y=0, width=1, height=1))
	dev.off()

	pdf(file = outfile2, width=30, height=10)
	print(ggdraw()+draw_plot(pie_only))
	dev.off()
}
############################################

dir.create("Violin_plots")
args <- commandArgs(TRUE)
if (length(args) < 3) {
  stop("Seurat object, gene name and Celltype must be supplied (input file).n", call.=FALSE)
} else
{
print("Generating Violinplots and Piecharts for gene expression across replicates...")
}

#args[1] = "/projects/anczukow-lab/neerja_projects/Integrate_expt1_expt2_only/Analysis_PC_10_mt10/RDS_files/seuratObject_rename_ident.rds"
#args[2] = "Pdk4"
#args[3] = "Luminal-AV"

mydata = readRDS(args[1])
gene_name = args[2]
Celltype = args[3]

#Call function to generate plots.
gene_Vln_Pie(mydata, gene_name, Celltype)

