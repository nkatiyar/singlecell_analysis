
## Using cinaR to do HPEA and hypergeometric tests
## Includes notes from Onur Karakaslar - creator of package 

#  Onur - I updated the pipeline a bit, this line is for updating yours too.
#install.packages("devtools")
#devtools::install_github("eonurk/cinaR", force = TRUE)

# Onur - Restart session to make it work properly
# Onur - If the run_enrichment line does not output any pathways, please close R and restart.
#.rs.restartR()

##Loading Packages
library(cinaR)
library(cinaRgenesets)
library(BiocParallel)
library(ggplot2)
register(DoparParam())

##Using cinaRgenesets to download available datasets
geneset <- cinaRgenesets::vp2008
wiki <- cinaRgenesets::wiki
dice <- cinaRgenesets::dice.major
wiki.inf <- cinaRgenesets::wiki.inf

##Uploading my own gmt lists from the MSIGDB website
#Hallmarks = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/h.all.v7.4.symbols.gmt')
Ontology_MSIG = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c5.all.v7.4.ontology_symbols.gmt')
Oncogenic_MSIG = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c6.all.v7.4.oncogenic_symbols.gmt')
Immunologic_MSIG = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c7.all.v7.4.immunologic_symbols.gmt')
Hallmarks = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/h.all.v7.4.symbols.gmt')
Curated = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c2.all.v7.4.symbols.gmt')
Wiki_MSIG = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c2.cp.wikipathways.v7.4.symbols.gmt')
KEGG_MSIG = cogena::gmt2list('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Genesets/c2.cp.kegg.v7.4.symbols.gmt')

##Subsetting the Hallmarks for Cancer related Pathways
Cancer_Hallmarks <- list(Hallmarks$HALLMARK_APICAL_SURFACE, Hallmarks$HALLMARK_APOPTOSIS, Hallmarks$HALLMARK_DNA_REPAIR, Hallmarks$HALLMARK_E2F_TARGETS,
                         Hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, Hallmarks$HALLMARK_ESTROGEN_RESPONSE_EARLY, Hallmarks$HALLMARK_ESTROGEN_RESPONSE_LATE,
                         Hallmarks$HALLMARK_FATTY_ACID_METABOLISM, Hallmarks$HALLMARK_G2M_CHECKPOINT, Hallmarks$HALLMARK_GLYCOLYSIS, Hallmarks$HALLMARK_HEDGEHOG_SIGNALING,
                         Hallmarks$HALLMARK_HYPOXIA, Hallmarks$HALLMARK_KRAS_SIGNALING_DN, Hallmarks$HALLMARK_KRAS_SIGNALING_UP, Hallmarks$HALLMARK_MITOTIC_SPINDLE,
                         Hallmarks$HALLMARK_MYC_TARGETS_V1, Hallmarks$HALLMARK_MTORC1_SIGNALING, Hallmarks$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                         Hallmarks$HALLMARK_P53_PATHWAY, Hallmarks$HALLMARK_INFLAMMATORY_RESPONSE, Hallmarks$HALLMARK_PI3K_AKT_MTOR_SIGNALING, 
                         Hallmarks$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, Hallmarks$HALLMARK_TGF_BETA_SIGNALING, Hallmarks$HALLMARK_TNFA_SIGNALING_VIA_NFKB, 
			 Hallmarks$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
names(Cancer_Hallmarks) <- c("HALLMARK_APICAL_SURFACE", "HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_E2F_TARGETS",
                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
                             "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_GLYCOLYSIS", "HALLMARK_HEDGEHOG_SIGNALING",
                             "HALLMARK_HYPOXIA", "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MITOTIC_SPINDLE",
                             "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                             "HALLMARK_P53_PATHWAY", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", 
                             "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_WNT_BETA_CATENIN_SIGNALING")

#### You can automate this part using a for loop or an apply function #####
# load your DE genes
Bcells_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Bcells_dir_no_filter/Bcells_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
DendMac_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Dendritic_Macrophages_dir_no_filter/Dendritic_Macrophages_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Tcells_naive_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Tcells_naive_dir_no_filter/Tcells_naive_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Tcells_mem_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Tcells_mem_dir_no_filter/Tcells_mem_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Plasma_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Plasma_dir_no_filter/Plasma_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
LumAV_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Luminal-AV_dir_no_filter/Luminal-AV_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
LumHS_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Luminal-HS_dir_no_filter/Luminal-HS_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Myo_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Myoepithelial_dir_no_filter/Myoepithelial_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Vascular_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Vascular_dir_no_filter/Vascular_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Fibro_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Fibroblasts_dir_no_filter/Fibroblasts_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)
Pericytes_Single_cell_Orig = read.delim("/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Overlap_analysis/Intersect_and_union_DE_singlecell_pseudobulk_DESeq2/No_filter/Pericytes_dir_no_filter/Pericytes_singlecell_vs_pseudobulk_merged_union_all_log2FC0.txt", sep="\t", header=T)

Bcells_Single_cell_select <- Bcells_Single_cell_Orig[Bcells_Single_cell_Orig$DE!="None",]
DendMac_Single_cell_select <- DendMac_Single_cell_Orig[DendMac_Single_cell_Orig$DE!="None",]
Tcells_naive_Single_cell_select <- Tcells_naive_Single_cell_Orig[Tcells_naive_Single_cell_Orig$DE!="None",]
Tcells_mem_Single_cell_select <- Tcells_mem_Single_cell_Orig[Tcells_mem_Single_cell_Orig$DE!="None",]
Plasma_Single_cell_select <- Plasma_Single_cell_Orig[Plasma_Single_cell_Orig$DE!="None",]
LumAV_Single_cell_select <- LumAV_Single_cell_Orig[LumAV_Single_cell_Orig$DE!="None",]
LumHS_Single_cell_select <- LumHS_Single_cell_Orig[LumHS_Single_cell_Orig$DE!="None",]
Myo_Single_cell_select <- Myo_Single_cell_Orig[Myo_Single_cell_Orig$DE!="None",]
Vascular_Single_cell_select <- Vascular_Single_cell_Orig[Vascular_Single_cell_Orig$DE!="None",]
Fibro_Single_cell_select <- Fibro_Single_cell_Orig[Fibro_Single_cell_Orig$DE!="None",]
Pericytes_Single_cell_select <- Pericytes_Single_cell_Orig[Pericytes_Single_cell_Orig$DE!="None",]

#Create a list of lists
Combined_Single_cell_Genes = list(B_Cells = Bcells_Single_cell_select,
				  Plasma = Plasma_Single_cell_select,
                                 DCs_Macrophages = DendMac_Single_cell_select, 
                                 Tcells_naive = Tcells_naive_Single_cell_select,
                                 Tcells_mem = Tcells_mem_Single_cell_select,
                                 Luminal_AV = LumAV_Single_cell_select,
                                 Fibroblasts = Fibro_Single_cell_select,
                                 Myoepithelial = Myo_Single_cell_select,
                                 Pericytes = Pericytes_Single_cell_select,
                                 Luminal_HS = LumHS_Single_cell_select,
                                 Vascular = Vascular_Single_cell_select)

#Change the naming of the columns (and filter for significant FDR)
#Sig_genes is for hypergeometric tests
#All_genes is for HPEA
Sig_genes <- function(x){
  x <- subset(x, (p_val_adj < 0.05 | padj < 0.05))
  colnames(x)[1] <- "gene_name"
  colnames(x)[14] <- "logFC"
  x
}

All_genes <- function(x){
  colnames(x)[1] <- "gene_name"
  colnames(x)[14] <- "logFC"
  x
}
Sig_Single_cell_Genes <- lapply(Combined_Single_cell_Genes, function(x){Sig_genes(x)})
All_Single_cell_Genes <- lapply(Combined_Single_cell_Genes, function(x){All_genes(x)})

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Single_cell_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Single_cell_Genes

##########3------------------####################3
###Altering dot_plot fn in cinaR to change labels
###Changing Hypergeometric test labels will not be as easy
dot_plot_HPEA <- function (results, fdr.cutoff = 0.1, filter.pathways = FALSE) 
{
  if (is.null(results[["Enrichment.Results"]])) {
    stop("Did you run the enrichment pipeline in cinaR? For more info ?cinaR")
  }
  adj.p <- column_label <- log2FoldChange <- module.name <- padj <- prcomp <- status <- NULL
  exp.atac <- FALSE
  atac.check <- unlist(lapply(results[["DA.results"]][["DA.peaks"]], 
                              function(x) {
                                ncol(x) == 3
                              }), use.names = FALSE)
  if (any(atac.check)) {
    message(">> cinaR was run for RNA-Seq experiments!")
    exp.atac <- FALSE
  }
  results <- results[["Enrichment.Results"]]
  df.plot <- dplyr::bind_rows(results, .id = "column_label")
  if ("padj" %in% colnames(df.plot)) {
    df.plot$status <- ifelse(df.plot$NES > 0, "Opening", "Closing")
    colnames(df.plot)[c(2, 4)] <- c("module.name", "adj.p")
  }

  if (!exp.atac) {
    df.plot[, "status"][df.plot[, "status"] == "Opening"] <- "Up"
    df.plot[, "status"][df.plot[, "status"] == "Closing"] <- "Down"
  }
  if (filter.pathways) {
    if (sum(df.plot$adj.p < fdr.cutoff) == 0) {
      stop("You can't filter because there are no pathways to be displayed!")
    }
    df.plot <- subset(df.plot, adj.p < fdr.cutoff)
  }
  plot.dot <- ggplot2::ggplot(df.plot, ggplot2::aes(x = column_label, 
                                                    y = module.name, size = ifelse(adj.p < fdr.cutoff, -log(adj.p), 
                                                                                   NA), color = status))
  plot.dot <- plot.dot + ggplot2::geom_point()
  plot.dot <- plot.dot + ggplot2::labs(x = "Contrast", y = "Pathways", 
                                       color = "Sign", size = "-log10(adj.p)", caption = paste0("FDR < ", 
                                                                                                fdr.cutoff))
  plot.dot <- plot.dot + ggplot2::scale_x_discrete(limits = c("Luminal_AV", "Luminal_HS", "Myoepithelial",
                                                              "DCs_Macrophages", "Tcells_naive", "Tcells_mem", "B_Cells", "Plasma", 
                                                              "Fibroblasts", "Pericytes", "Vascular"))
  plot.dot <- plot.dot + ggplot2::scale_color_manual(values = color_values)
  plot.dot <- plot.dot + ggplot2::theme(text = element_text(size=12),
                                        axis.text.x = element_text(angle=45, hjust=1)) 
  return(plot.dot)
}

###Runing dot_plot fn
##Using cinaR to do enrichment analysis
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = wiki,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "HPEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

#results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
#                                       geneset = wiki,
#                                       experiment.type = "RNA-Seq",
#                                       enrichment.method = "HPEA",
#                                       reference.genome = "mm10",
#                                       enrichment.FDR.cutoff = 0.1)

pdf("Dotplot_wiki_cinaR_HPEA_singlecell_sig_only.pdf", width=15)
dot_plot_HPEA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_sig$Enrichment.Results$Luminal_AV, "LumAV_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Luminal_HS, "LumHS_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Myoepithelial, "Myo_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$B_Cells, "Bcells_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Plasma, "Plasma_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Pericytes, "Pericytes_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Vascular, "Vascular_enrich_HPEA_wiki_sig_only.txt", sep="\t", quote=FALSE)
 
####################-------------------################################

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Single_cell_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Single_cell_Genes

##Using cinaR to do enrichment analysis
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = Cancer_Hallmarks,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "HPEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

#results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
#                                       geneset = Cancer_Hallmarks,
#                                       experiment.type = "RNA-Seq",
#                                       enrichment.method = "HPEA",
#                                       reference.genome = "mm10",
#                                       enrichment.FDR.cutoff = 0.1)

pdf("Dotplot_cancer_hallmarks_HPEA_singlecell_sig_only.pdf", width=15)
dot_plot_HPEA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

#pdf("Dotplot_cancer_hallmarks_HPEA_singlecell_all.pdf", width=15)
#dot_plot_HPEA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
#dev.off()

write.table(results_sig$Enrichment.Results$Luminal_AV, "LumAV_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Luminal_HS, "LumHS_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Myoepithelial, "Myo_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$B_Cells, "Bcells_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Plasma, "Plasma_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Pericytes, "Pericytes_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Vascular, "Vascular_enrich_HPEA_cancer_hallmarks_sig_only.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Single_cell_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Single_cell_Genes

Wiki_MSIG_subset <- list(Wiki_MSIG$WP_IL7_SIGNALING_PATHWAY, Wiki_MSIG$WP_IL5_SIGNALING_PATHWAY, Wiki_MSIG$WP_IL4_SIGNALING_PATHWAY, Wiki_MSIG$WP_IL3_SIGNALING_PATHWAY, Wiki_MSIG$WP_IL2_SIGNALING_PATHWAY, Wiki_MSIG$WP_IL18_SIGNALING_PATHWAY, Wiki_MSIG$WP_TYPE_II_INTERFERON_SIGNALING_IFNG, Wiki_MSIG$WP_TGFBETA_SIGNALING_PATHWAY, Wiki_MSIG$WP_TCELL_RECEPTOR_AND_COSTIMULATORY_SIGNALING, Wiki_MSIG$WP_TCELL_ANTIGEN_RECEPTOR_TCR_SIGNALING_PATHWAY)

names(Wiki_MSIG_subset) <- c("WP_IL7_SIGNALING_PATHWAY", "WP_IL5_SIGNALING_PATHWAY", "WP_IL4_SIGNALING_PATHWAY", "WP_IL3_SIGNALING_PATHWAY", "WP_IL2_SIGNALING_PATHWAY", "WP_IL18_SIGNALING_PATHWAY", "WP_TYPE_II_INTERFERON_SIGNALING_IFNG", "WP_TGFBETA_SIGNALING_PATHWAY", "WP_TCELL_RECEPTOR_AND_COSTIMULATORY_SIGNALING", "WP_TCELL_ANTIGEN_RECEPTOR_TCR_SIGNALING_PATHWAY")

##Using cinaR to do enrichment analysis
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = Wiki_MSIG_subset,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "HPEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

pdf("Dotplot_wiki_msigdb_HPEA_singlecell_sig_only.pdf", width=15)
dot_plot_HPEA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_sig$Enrichment.Results$Luminal_AV, "LumAV_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Luminal_HS, "LumHS_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Myoepithelial, "Myo_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$B_Cells, "Bcells_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Plasma, "Plasma_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Pericytes, "Pericytes_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Vascular, "Vascular_enrich_HPEA_wiki_msigdb_sig_only.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Single_cell_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Single_cell_Genes

##Using cinaR to do enrichment analysis
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = KEGG_MSIG,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "HPEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

#results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
#                                       geneset = KEGG_MSIG,
#                                       experiment.type = "RNA-Seq",
#                                       enrichment.method = "HPEA",
#                                       reference.genome = "mm10",
#                                       enrichment.FDR.cutoff = 0.1)


pdf("Dotplot_kegg_msigdb_HPEA_singlecell_sig_only.pdf", width=15, height=20)
dot_plot_HPEA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

#pdf("Dotplot_kegg_msigdb_HPEA_singlecell_all.pdf", width=15, height=20)
#dot_plot_HPEA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
#dev.off()

write.table(results_sig$Enrichment.Results$Luminal_AV, "LumAV_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Luminal_HS, "LumHS_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Myoepithelial, "Myo_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$B_Cells, "Bcells_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Plasma, "Plasma_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Pericytes, "Pericytes_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Vascular, "Vascular_enrich_HPEA_kegg_msigdb_sig_only.txt", sep="\t", quote=FALSE)

#######################3--------------------###########################3
# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Single_cell_Genes

#results_all <- list()
#results_all[["DA.peaks"]] <- All_Single_cell_Genes

##Using cinaR to do enrichment analysis using VP2008 geneset.
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = geneset,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "HPEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

#results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
#                                       geneset = geneset,
#                                       experiment.type = "RNA-Seq",
#                                       enrichment.method = "HPEA",
#                                       reference.genome = "mm10",
#                                       enrichment.FDR.cutoff = 0.1)


pdf("Dotplot_VP2008_HPEA_singlecell_sig_only.pdf", width=15, height=20)
dot_plot_HPEA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

#pdf("Dotplot_VP2008_msigdb_HPEA_singlecell_all.pdf", width=15, height=20)
#dot_plot_HPEA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
#dev.off()

write.table(results_sig$Enrichment.Results$Luminal_AV, "LumAV_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Luminal_HS, "LumHS_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Myoepithelial, "Myo_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$B_Cells, "Bcells_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Plasma, "Plasma_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Pericytes, "Pericytes_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)
write.table(results_sig$Enrichment.Results$Vascular, "Vascular_enrich_HPEA_VP2008_sig_only.txt", sep="\t", quote=FALSE)


