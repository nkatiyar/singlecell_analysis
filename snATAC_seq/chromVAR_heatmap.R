library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(universalmotif)
library(dendextend)
library(grid)
library(gridExtra)
library(ggside)
library(patchwork)

# read in the differential motifs dataframe from Signac Vignette
chromVAR.signature.motifs <- read.csv("FinalMotifs1.txt", sep = '\t', header = T)
colnames(chromVAR.signature.motifs) <- c("ID", "p_val", "log_FC_Change", "pct.1", "pct.2", "adj_p_val", "CellType")
chromVAR.signature.motifs$CellType <- gsub(".txt","",chromVAR.signature.motifs$CellType)
chromVAR.signature.motifs$log_FC_Change <- as.numeric(chromVAR.signature.motifs$log_FC_Change)
chromVAR.signature.motifs <- chromVAR.signature.motifs[(chromVAR.signature.motifs$log_FC_Change > 0),]

# create a motif database in order to cluster motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# convert motif ids to motif names using database
chromVAR.signature.motifs <- chromVAR.signature.motifs %>%
  rowwise() %>%
  mutate(motif = pfm[[ID]]@name) %>%
  as.data.frame()

# select the top motifs you want to plot for each group
top_motifs <- chromVAR.signature.motifs %>%
  dplyr::group_by(CellType) %>%
  slice_min(order_by = adj_p_val, n = 10, with_ties = F) %>%
  as.data.frame() %>%
  dplyr::select(motif,ID)

top_motifs <- top_motifs %>%
  distinct()

# get the top motifs as Universal Motif objects so that we can filter and cluster them
motifs <- universalmotif::convert_motifs(pfm[top_motifs$ID])

# just keep human motifs
#motifs <- filter_motifs(motifs, organism = "Mus musculus")

# plots a hierarchical clustering of the motifs to give an idea
#motif_tree(motifs, labels = "name", layout = "rectangular", linecol = "organism")

# find euclidian distance of motifs
comp.matrix <- universalmotif::compare_motifs(motifs, method = "EUCL")

# convert to suitable format for dhc
comp.dist <- as.dist(comp.matrix)

# hierarchical clustering of the motifs for pheatmap
dhc <- hclust(comp.dist)

# just get the top motifs from results and also convert adjusted p value to -log for plotting
df <- chromVAR.signature.motifs %>%
  mutate(CellType = factor(CellType, levels = c("NaiveCD4_LEF1", "CD4_CCR7_JUN", "MemoryCD4_CD44", "CD8_CCR7", "CD8_CCR9", "CD8_GZMK", "CD8_GZMM", "IL7R_ICOS", "NK"))) %>%
  mutate(Log_Adj_P = -log10(adj_p_val)) %>%
  mutate(Log_Adj_P = ifelse(Log_Adj_P > 400, 400, Log_Adj_P)) %>%
  filter(motif %in% dhc$labels) %>%
  as.data.frame()

# convert data frame to matrix for pheatmap
a <- df[,c("CellType","motif","Log_Adj_P")] %>%
  distinct(CellType, motif, .keep_all = T) %>%
  spread(key = "CellType", value = "Log_Adj_P")

rownames(a) <- a$motif
b <- a[,-1]
c <- as.data.frame(t(b))
# pheatmap. a[dhc$labels,-1] is really important as the clustering will not be correct if motifs are not given in the same order as dhc object. Adds stars to boxes based on thresholds
p <- pheatmap::pheatmap(c[,dhc$labels], 
                        cluster_rows = F,
                        cluster_cols = as.hclust(dhc),
                        cutree_cols = 24, scale = "none", na_col = "grey",
                        color = colorRampPalette(c("lightblue1","#56B1F7", "#132B43"))(50),
                        #angle_col = 45, 
                        #display_numbers = matrix(ifelse(a[dhc$labels,-1] > 30, "***", ifelse(a[dhc$labels,-1] > 10, "**",ifelse(a[dhc$labels,-1] > 1.3, "*","" ) )), nrow(a[dhc$labels,-1])),
                        fontsize_number = 12,
                        show_colnames = T, legend = T)
pdf("./TCells_MotifClustering.pdf", height = 4, width = 15)
p
dev.off()


