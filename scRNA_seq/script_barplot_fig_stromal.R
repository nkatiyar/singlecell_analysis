

library(ggplot2)
library(dplyr)
dat_DE = read.table("stromal_DE_genes_summary.txt", sep="\t", header=T)
dat_DE$celltype <- factor(dat_DE$celltype, levels=c("Fibroblast", "Vascular", "Pericyte"))

pdf("Plot_DE_stromal_up_down.pdf")
ggplot(dat_DE, aes(x=celltype, y=value, fill=regulation)) + scale_fill_manual("regulation", values = c("Up" = "indianred1", "Down" = "steelblue1")) + geom_col(width=0.4, colour="black") + ylab("Differential Genes") 
dev.off()

###########---------------------------------------####################

dat_DA = read.table("stromal_DA_genes_summary.txt", sep="\t", header=T)
dat_DA$celltype <- factor(dat_DA$celltype, levels=c("Fibroblast", "Vascular", "Pericyte"))

pdf("Plot_DA_stromal_up_down.pdf")
ggplot(dat_DA, aes(x=celltype, y=value, fill=regulation)) + scale_fill_manual("regulation", values = c("Open" = "indianred1", "Close" = "steelblue1")) + geom_col(width=0.4, colour="black") + ylab("Genes for differential Peaks")
dev.off()

