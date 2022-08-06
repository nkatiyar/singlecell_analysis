

library(ggplot2)
library(dplyr)
data_enrich <- read.table("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Scripts/Enrichment_analysis/HPEA/All_samples_lab_meeting_reanalysis/Fibroblasts_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", header=T)
rownames(data_enrich) <- NULL

data_enrich$score = -log(data_enrich$adj.p)
#data_enrich <- data_enrich[order(data_enrich$status,data_enrich$score),]
data_enrich_open <- data_enrich %>% filter(status == "Opening")
data_enrich_open <- data_enrich_open[order(data_enrich_open$status,-data_enrich_open$score),]

data_enrich_close <- data_enrich %>% filter(status == "Closing")
data_enrich_close <- data_enrich_close[order(data_enrich_close$status,data_enrich_close$score),]

data_enrich_all <- rbind(data_enrich_open, data_enrich_close)
print(data_enrich_all)

pdf("Enriched_pathways_fibro_DA_genes.pdf")
ggplot(data_enrich_all) +
  # set overall appearance of the plot
  theme_bw() +
  # Define the dependent and independent variables
  aes(x = module.name, y = score) +
  # From the defined variables, create a vertical bar chart
  geom_col(position = "stack", aes(fill = status)) +
  # Set the legend, main titles, and axis titles
  ggtitle("Cancer Hallmarks") +
  xlab("Pathways") +
  ylab("-log(FDR)") +
  # Flip the x and y axes
  coord_flip()

dev.off()

