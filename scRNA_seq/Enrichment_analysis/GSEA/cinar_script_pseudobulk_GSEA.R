
## Using cinaR to do GSEA and hypergeometric tests
## Includes notes from Onur Karakaslar - creator of package 

#  Onur - I updated the pipeline a bit, this line is for updating yours too.
install.packages("devtools")
devtools::install_github("eonurk/cinaR", force = TRUE)

# Onur - Restart session to make it work properly
# Onur - If the run_enrichment line does not output any pathways, please close R and restart.
.rs.restartR()

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
                         Hallmarks$HALLMARK_MYC_TARGETS_V1, Hallmarks$HALLMARK_MYC_TARGETS_V2, Hallmarks$HALLMARK_MTORC1_SIGNALING, Hallmarks$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                         Hallmarks$HALLMARK_P53_PATHWAY, Hallmarks$HALLMARK_INFLAMMATORY_RESPONSE, Hallmarks$HALLMARK_PI3K_AKT_MTOR_SIGNALING, 
                         Hallmarks$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, Hallmarks$HALLMARK_TGF_BETA_SIGNALING, Hallmarks$HALLMARK_TNFA_SIGNALING_VIA_NFKB, 
                         Hallmarks$HALLMARK_UNFOLDED_PROTEIN_RESPONSE, Hallmarks$HALLMARK_UV_RESPONSE_UP, Hallmarks$HALLMARK_UV_RESPONSE_DN, Hallmarks$HALLMARK_WNT_BETA_CATENIN_SIGNALING)
names(Cancer_Hallmarks) <- c("HALLMARK_APICAL_SURFACE", "HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR", "HALLMARK_E2F_TARGETS",
                             "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
                             "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_GLYCOLYSIS", "HALLMARK_HEDGEHOG_SIGNALING",
                             "HALLMARK_HYPOXIA", "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MITOTIC_SPINDLE",
                             "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                             "HALLMARK_P53_PATHWAY", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", 
                             "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
                             "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UV_RESPONSE_UP", "HALLMARK_UV_RESPONSE_DN", 
                             "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

#### You can automate this part using a for loop or an apply function #####
# load your DE genes
Bcells_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Bcells_mm10_18mths_vs_mm10_3mths_all_genes.csv')
DendMac_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Dendritic_or_Macrophages_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Tcells_naive_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Tcells_naive_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Tcells_mem_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Tcells_mem_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Plasma_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Plasma_mm10_18mths_vs_mm10_3mths_all_genes.csv')
LumAV_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Luminal-AV_mm10_18mths_vs_mm10_3mths_all_genes.csv')
LumHS_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Luminal-HS_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Myo_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Myoepithelial_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Vascular_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Vascular_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Fibro_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Fibroblasts_mm10_18mths_vs_mm10_3mths_all_genes.csv')
Pericytes_Pseudobulk_Orig = read.csv('/projects/anczukow-lab/neerja_projects/New_scRNA_ATAC_seq_analysis_PC10/scRNAseq/Pseudobulk_DE_analysis/DESeq2/pairwise/Pericytes_mm10_18mths_vs_mm10_3mths_all_genes.csv')

#Create a list of lists
Combined_Pseudobulk_Genes = list(B_Cells = Bcells_Pseudobulk_Orig,
                                 DCs_Macrophages = DendMac_Pseudobulk_Orig, 
                                 Tcells_naive = Tcells_naive_Pseudobulk_Orig,
                                 Tcells_mem = Tcells_mem_Pseudobulk_Orig,
				 Plasma = Plasma_Pseudobulk_Orig,
                                 Luminal_AV = LumAV_Pseudobulk_Orig,
                                 Fibroblasts = Fibro_Pseudobulk_Orig,
                                 Myoepithelial = Myo_Pseudobulk_Orig,
                                 Pericytes = Pericytes_Pseudobulk_Orig,
                                 Luminal_HS = LumHS_Pseudobulk_Orig,
                                 Vascular = Vascular_Pseudobulk_Orig)

#Change the naming of the columns (and filter for significant FDR)
#Sig_genes is for hypergeometric tests
#All_genes is for GSEA
Sig_genes <- function(x){
  x <- subset(x, padj < 0.05)
  colnames(x)[1] <- "gene_name"
  colnames(x)[3] <- "logFC"
  x
}

All_genes <- function(x){
  colnames(x)[1] <- "gene_name"
  colnames(x)[3] <- "logFC"
  x
}
Sig_Pseudobulk_Genes <- lapply(Combined_Pseudobulk_Genes, function(x){Sig_genes(x)})
All_Pseudobulk_Genes <- lapply(Combined_Pseudobulk_Genes, function(x){All_genes(x)})

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Pseudobulk_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes

results_all_curated <- results_all
##Using cinaR to do enrichment analysis
results_all_curated[["Enrichment.Results"]] <- run_enrichment(results_all_curated,
                                       geneset = Curated,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)

results_all_hallmarks <- results_all
##Using cinaR to do enrichment analysis
results_all_hallmarks[["Enrichment.Results"]] <- run_enrichment(results_all_hallmarks,
                                       geneset = Hallmarks,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)


###Altering dot_plot fn in cinaR to change labels
###Changing Hypergeometric test labels will not be as easy
dot_plot_BA <- function (results, fdr.cutoff = 0.1, filter.pathways = FALSE) 
{
  if (is.null(results[["Enrichment.Results"]])) {
    stop("Did you run the enrichment pipeline in cinaR? For more info ?cinaR")
  }
  adj.p <- column_label <- log2FoldChange <- module.name <- padj <- prcomp <- status <- NULL
  exp.atac <- TRUE
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
    df.plot$status <- ifelse(df.plot$NES > 0, "Up", 
                             "Down")
    colnames(df.plot)[c(2, 4)] <- c("module.name", "adj.p")
  }
  if (!exp.atac) {
    df.plot[, "status"][df.plot[, "status"] == "Up"] <- "Up"
    df.plot[, "status"][df.plot[, "status"] == "Down"] <- "Down"
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
# create your differential analyses object
results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes

results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
                                       geneset = Cancer_Hallmarks,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)

pdf("Dotplot_cancer_hallmarks_GSEA_pseudobulk_all.pdf", width=15)
dot_plot_BA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_all$Enrichment.Results$Luminal_AV, "LumAV_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Luminal_HS, "LumHS_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Myoepithelial, "Myo_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$B_Cells, "Bcells_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Plasma, "Plasma_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Pericytes, "Pericytes_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Vascular, "Vascular_enrich_GSEA_cancer_hallmarks.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes

##Using cinaR to do enrichment analysis
results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
                                       geneset = Wiki_MSIG,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)

pdf("Dotplot_wiki_msigdb_hallmarks_GSEA_pseudobulk_all.pdf", width=15, height=30)
dot_plot_BA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_all$Enrichment.Results$Luminal_AV, "LumAV_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Luminal_HS, "LumHS_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Myoepithelial, "Myo_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$B_Cells, "Bcells_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Plasma, "Plasma_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Pericytes, "Pericytes_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Vascular, "Vascular_enrich_GSEA_wiki_msigdb.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes
##Using cinaR to do enrichment analysis
results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
                                       geneset = KEGG_MSIG,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)

pdf("Dotplot_kegg_msigdb_hallmarks_GSEA_pseudobulk_all.pdf", width=15, height=30)
dot_plot_BA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_all$Enrichment.Results$Luminal_AV, "LumAV_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Luminal_HS, "LumHS_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Myoepithelial, "Myo_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$B_Cells, "Bcells_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Plasma, "Plasma_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Pericytes, "Pericytes_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Vascular, "Vascular_enrich_GSEA_kegg_msigdb.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes
##Using cinaR to do enrichment analysis
results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
                                       geneset = wiki,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 1)

pdf("Dotplot_wiki_Onur_hallmarks_GSEA_pseudobulk_all.pdf", width=15, height=30)
dot_plot_BA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_all$Enrichment.Results$Luminal_AV, "LumAV_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Luminal_HS, "LumHS_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Myoepithelial, "Myo_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$B_Cells, "Bcells_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Plasma, "Plasma_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Pericytes, "Pericytes_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Vascular, "Vascular_enrich_GSEA_Onur_wiki.txt", sep="\t", quote=FALSE)

# create your differential analyses object
results_sig <- list()
results_sig[["DA.peaks"]] <- Sig_Pseudobulk_Genes

results_all <- list()
results_all[["DA.peaks"]] <- All_Pseudobulk_Genes

##Using cinaR to do enrichment analysis using VP2008 geneset.
results_sig[["Enrichment.Results"]] <- run_enrichment(results_sig,
                                       geneset = geneset,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)

results_all[["Enrichment.Results"]] <- run_enrichment(results_all,
                                       geneset = geneset,
                                       experiment.type = "RNA-Seq",
                                       enrichment.method = "GSEA",
                                       reference.genome = "mm10",
                                       enrichment.FDR.cutoff = 0.1)


pdf("Dotplot_VP2008_GSEA_pseudobulk_sig_only.pdf", width=15, height=20)
dot_plot_BA(results_sig, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

pdf("Dotplot_VP2008_msigdb_GSEA_pseudobulk_all.pdf", width=15, height=20)
dot_plot_BA(results_all, fdr.cutoff = 0.1, filter.pathways = T)
dev.off()

write.table(results_all$Enrichment.Results$Luminal_AV, "LumAV_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Luminal_HS, "LumHS_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Myoepithelial, "Myo_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$DCs_Macrophages, "Dendritic_Macrophages_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_naive, "Tcells_naive_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Tcells_mem, "Tcells_mem_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$B_Cells, "Bcells_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Plasma, "Plasma_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Fibroblasts, "Fibroblasts_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Pericytes, "Pericytes_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)
write.table(results_all$Enrichment.Results$Vascular, "Vascular_enrich_GSEA_VP2008.txt", sep="\t", quote=FALSE)

