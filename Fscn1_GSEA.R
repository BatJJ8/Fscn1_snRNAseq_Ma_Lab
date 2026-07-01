# Load libraries
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(tidyverse)
library(writexl)
library(ggrepel)
library(tools)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(stringr)
library(scales)
library(SeuratWrappers)
library(glmGamPoi)
library(presto)
library(hdf5r)
library(readxl)
library(openxlsx)
library(future)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(DO.db)
library(biomaRt)
library(AnnotationDbi)
library(fgsea)
library(msigdbr)
library(msigdbdf)
library(clusterProfiler)
library(gprofiler2)
library(DOSE)
library(dittoSeq)

theme_set(theme_minimal())

options(future.globals.maxSize = 2e10)  

set.seed(888888)

# Load the DEGs
wb_path <- "KO_vs_WT_DGE_MAST.xlsx"   

mast_res  <- read.xlsx(wb_path, sheet = "MAST_All_DEGs", rowNames = TRUE)

# GSEA using fgsea

# Retrieve mouse Hallmark gene sets from msigdbr
msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
# Convert the msigdbr data frame into a list where each element is a gene set vector
pathways <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)

run_fgsea_preranked <- function(markers_df, pathways, ranking_method = "signed_p", top_n = 10) {
  
  required_cols <- c("p_val", "avg_log2FC")
  missing_cols <- setdiff(required_cols, colnames(markers_df))
  if (length(missing_cols) > 0) {
    stop("markers_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Drop rows with missing essentials
  markers_df <- markers_df[!is.na(markers_df$p_val) & !is.na(markers_df$avg_log2FC), , drop = FALSE]
  
  # Compute ranking vector using all markers
  if (ranking_method == "signed_p") {
    ranks <- sign(markers_df$avg_log2FC) * -log10(markers_df$p_val + 1e-10) # add pseudocount
  } else if (ranking_method == "fold_change") {
    ranks <- markers_df$avg_log2FC
  } else {
    stop("Invalid ranking method specified. Choose 'signed_p' or 'fold_change'.")
  }
  
  names(ranks) <- rownames(markers_df)
  
  # Clean
  keep <- !is.na(names(ranks)) & names(ranks) != "" & !is.na(ranks)
  ranks <- ranks[keep]
  ranks <- ranks[!duplicated(names(ranks))]
  
  ranks <- sort(ranks, decreasing = TRUE)
  
  set.seed(123)
  fgseaRes <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks)
  fgseaResTidy <- fgseaRes[order(fgseaRes$pval), ]
  
  # Top pathways + plots
  top_pathways <- head(fgseaResTidy, top_n)
  
  enrichment_plots <- lapply(top_pathways$pathway, function(pw) {
    fgsea::plotEnrichment(pathways[[pw]], ranks) + ggplot2::labs(title = pw)
  })
  
  combined_plot <- patchwork::wrap_plots(enrichment_plots, ncol = 2)
  
  list(fgseaRes = fgseaResTidy, plot = combined_plot)
}


# Helper function to compute non-leading edge genes
get_non_leading_edge <- function(pathway, leadingEdge, filtered_genes, pathways_list) {
  # Retrieve full gene set for pathway
  full_gene_set <- pathways_list[[as.character(pathway)]]
  if (is.null(full_gene_set) || length(full_gene_set) == 0) {
    return("")
  }
  # Compute overlap between filtered marker genes and full gene set
  overlap <- intersect(filtered_genes, full_gene_set)
  if (length(overlap) == 0) {
    return("")
  }
  # Leading edge genes
  le <- if (is.null(leadingEdge)) character(0) else unlist(leadingEdge)
  # Non-leading edge genes
  non_le <- setdiff(overlap, le)
  paste(non_le, collapse = ", ")
}


# Preranked GSEA on all DEGs
gsea_overall <- run_fgsea_preranked(mast_res, pathways, ranking_method = "signed_p", top_n = 5)
print(gsea_overall$plot)


# Get top 10 significant pathways (sorted by p-value)
top10 <- gsea_overall$fgseaRes %>% 
  arrange(pval) %>% 
  head(10) %>% 
  mutate(sig = -log10(padj))  # transform padj for size mapping

# Reformat pathway names
top10 <- top10 %>%
  mutate(pathway = gsub("^HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = toTitleCase(tolower(pathway))
  )

# Categorize pathways as positively enriched or negatively enriched
top10 <- top10 %>%
  mutate(direction = ifelse(NES > 0, "Enriched", "Negatively Enriched"))

gsea_df <- gsea_overall$fgseaRes %>%
  mutate(
    # Convert each leadingEdge entry to a comma-separated string
    leading_edge_str = map_chr(leadingEdge, ~ paste(unlist(.x), collapse = ", ")),
    # Compute non-leading edge genes for each pathway using helper function
    non_leading_edge_str = map2_chr(pathway, leadingEdge, 
                                    ~ get_non_leading_edge(.x, .y, rownames(mast_res), pathways))
  ) %>%
  select(pathway, pval, padj, NES, leading_edge_str, non_leading_edge_str)


write_xlsx(gsea_df, "fgsea_combined.xlsx")