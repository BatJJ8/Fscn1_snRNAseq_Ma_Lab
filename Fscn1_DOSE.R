# Load libraries
library(tidyverse)
library(scales)
library(openxlsx)
library(writexl)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(Orthology.eg.db)
library(clusterProfiler)
library(DOSE)
library(RColorBrewer)
library(viridis)

theme_set(theme_minimal())

set.seed(888888)

# Load the DEGs
wb_path <- "KO_vs_WT_DGE_MAST.xlsx"   

mast_res  <- read.xlsx(wb_path, sheet = "MAST_All_DEGs", rowNames = TRUE)


# Create custom mapping function
mapIt <- function(mouseids, horg, morg, orth){
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- AnnotationDbi::select(orth, mouseg, "Homo_sapiens","Mus_musculus")
  names(mapped) <- c("Mus_egid","Homo_egid")
  husymb <- AnnotationDbi::select(horg, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  return(data.frame(Mus_symbol = mouseids,
                    mapped,
                    Homo_symbol = husymb[,2]))
}

sig_markers <- mast_res[mast_res$p_val < 0.05, ]
mouse_genes <- unique(rownames(sig_markers))
cat("Number of significant mouse genes:", length(mouse_genes), "\n")


# Map mouse gene symbols to human gene symbols
mapped_genes <- mapIt(mouse_genes, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
# Extract unique human gene symbols
human_genes <- unique(mapped_genes$Homo_symbol)
# Remove NA or blank entries
human_genes <- human_genes[!is.na(human_genes) & human_genes != ""]
cat("Number of mapped human genes from DEGs:", length(human_genes), "\n")


# Convert the human gene symbols to human Entrez IDs
gene_ids <- bitr(human_genes, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)
if(is.null(gene_ids) || nrow(gene_ids) == 0){
  stop("No human Entrez IDs could be mapped from the DEGs.")
}
cat("Number of human Entrez IDs from DEGs:", nrow(gene_ids), "\n")


# Disease Ontology Enrichment Analysis using DOSE

# Run DOSE using enrichDO for all significant DEGs
do_results <- enrichDO(gene         = gene_ids$ENTREZID,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05,
                       minGSSize    = 5,
                       maxGSSize    = 500)


# Separate upregulated and downregulated markers
up_markers   <- sig_markers[sig_markers$avg_log2FC > 0, ]
down_markers <- sig_markers[sig_markers$avg_log2FC < 0, ]

# Get unique mouse gene symbols for each group
up_mouse_genes   <- unique(rownames(up_markers))
down_mouse_genes <- unique(rownames(down_markers))

cat("Number of significant upregulated mouse genes:", length(up_mouse_genes), "\n")
cat("Number of significant downregulated mouse genes:", length(down_mouse_genes), "\n")


# Map mouse gene symbols to human gene symbols
# Upregulated genes mapping
mapped_genes_up <- mapIt(up_mouse_genes, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
up_human_genes  <- unique(mapped_genes_up$Homo_symbol)
up_human_genes  <- up_human_genes[!is.na(up_human_genes) & up_human_genes != ""]
cat("Number of mapped human genes from upregulated DEGs:", length(up_human_genes), "\n")

# Downregulated genes mapping
mapped_genes_down <- mapIt(down_mouse_genes, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
down_human_genes  <- unique(mapped_genes_down$Homo_symbol)
down_human_genes  <- down_human_genes[!is.na(down_human_genes) & down_human_genes != ""]
cat("Number of mapped human genes from downregulated DEGs:", length(down_human_genes), "\n")


# 3. Convert the human gene symbols to human Entrez IDs
# Upregulated genes Entrez IDs
gene_ids_up <- bitr(up_human_genes, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)
if (is.null(gene_ids_up) || nrow(gene_ids_up) == 0) {
  stop("No human Entrez IDs could be mapped from the upregulated DEGs.")
}
cat("Number of human Entrez IDs from upregulated DEGs:", nrow(gene_ids_up), "\n")

# Downregulated genes Entrez IDs
gene_ids_down <- bitr(down_human_genes, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)
if (is.null(gene_ids_down) || nrow(gene_ids_down) == 0) {
  stop("No human Entrez IDs could be mapped from the downregulated DEGs.")
}
cat("Number of human Entrez IDs from downregulated DEGs:", nrow(gene_ids_down), "\n")


# Run DO Enrichment Analysis using enrichDO
# Enrichment for upregulated DEGs
do_results_up <- enrichDO(gene         = gene_ids_up$ENTREZID,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05,
                          minGSSize    = 5,
                          maxGSSize    = 500)

if (nrow(as.data.frame(do_results_up)) > 0) {
  cat("Top DO terms:\n")
  print(head(do_results_up))
  
  p <- dotplot(do_results_up, showCategory = 10) +
    ggtitle("DO Enrichment for All Upregulated DEGs") +
    theme_classic(base_size = 16) +
    guides(
      fill = guide_colorbar(order = 1, barwidth = unit(8, "cm"), barheight = unit(0.5, "cm")),
      size = guide_legend(order = 2)
    ) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p)
} else {
  cat("No enriched DO terms found for the significant DEGs.\n")
}


# Convert the enrichResult to a data frame
do_results_up_df <- as.data.frame(do_results_up)

# Convert GeneRatio and BgRatio from "x/y" to numeric values
do_results_up_df$GeneRatioNumeric <- sapply(do_results_up_df$GeneRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

do_results_up_df$BgRatioNumeric <- sapply(do_results_up_df$BgRatio, function(x) {
  parts <- unlist(strsplit(x, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
})

# Calculate fold enrichment
do_results_up_df$FoldEnrichment <- do_results_up_df$GeneRatioNumeric / do_results_up_df$BgRatioNumeric

# Dot plot using fold enrichment as the x-axis
p <- ggplot(do_results_up_df, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
  geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = "black", stroke = 1.5) +
  scale_size_continuous(range = c(4, 10)) +
  scale_fill_gradient(
    low = "red", 
    high = "blue",
    guide = guide_colorbar(
      barwidth = unit(6, "cm"),  
      barheight = unit(0.5, "cm")  
    )
  ) +
  labs(title = "DO Enrichment for Upregulated DEGs",
       x = "Fold Enrichment",
       y = NULL,
       color = "p.adjust") +
  guides(
    fill = guide_colorbar(order = 1, barwidth = unit(6, "cm"), barheight = unit(0.5, "cm")),
    size = guide_legend(order = 2)
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line      = element_line(linewidth = 1),
    axis.ticks     = element_line(linewidth = 1)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 5))

print(p)

# Enrichment for downregulated DEGs
do_results_down <- enrichDO(gene         = gene_ids_down$ENTREZID,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 0.05,
                            minGSSize    = 5,
                            maxGSSize    = 500)

# Check the results
if (nrow(as.data.frame(do_results_down)) > 0) {
  cat("Top DO terms:\n")
  print(head(do_results_down))
  print(dotplot(do_results_down, showCategory = 10) + ggtitle("DO Enrichment for All Downregulated DEGs"))
} else {
  cat("No enriched DO terms found for the significant DEGs.\n")
}

do_results_df   <- as.data.frame(do_results)
do_results_up_df    <- as.data.frame(do_results_up)
do_results_down_df  <- as.data.frame(do_results_down)



# Map Entrez IDs to gene symbols
all_entrez_ids <- unique(unlist(strsplit(do_results@result$geneID, "/")))
id2symbol <- bitr(all_entrez_ids,
                  fromType = "ENTREZID",
                  toType   = "SYMBOL",
                  OrgDb    = org.Hs.eg.db) %>%
  distinct(ENTREZID, SYMBOL)

# Helper to collapse genes and match order of significance
collapse_genes <- function(enrich_result, mapping_df) {
  ordered_terms <- enrich_result@result %>%
    arrange(p.adjust) %>%
    select(ID, Description, p.adjust)
  
  # collapse genes
  collapsed <- enrich_result@result %>%
    select(ID, Description, geneID) %>%
    separate_rows(geneID, sep = "/") %>%
    left_join(mapping_df, by = c("geneID" = "ENTREZID")) %>%
    mutate(GeneSymbol = ifelse(!is.na(SYMBOL), SYMBOL, geneID)) %>%
    group_by(ID, Description) %>%
    summarise(GeneSymbols = str_c(unique(GeneSymbol), collapse = ", ")) %>%
    ungroup()
  
  # restore the order of significance
  collapsed <- left_join(ordered_terms, collapsed, by = c("ID", "Description")) %>%
    arrange(p.adjust) %>%
    select(ID, Description, GeneSymbols)
  
  return(collapsed)
}

genes_all  <- collapse_genes(do_results, id2symbol)
genes_up   <- collapse_genes(do_results_up, id2symbol)
genes_down <- collapse_genes(do_results_down, id2symbol)

sheets <- list(
  "All_DEGs_DO"     = as.data.frame(do_results),
  "Up_DEGs_DO"      = as.data.frame(do_results_up),
  "Down_DEGs_DO"    = as.data.frame(do_results_down),
  "All_DEGs_Genes"  = genes_all,
  "Up_DEGs_Genes"   = genes_up,
  "Down_DEGs_Genes" = genes_down
)

write_xlsx(sheets, "DO_enrichment.xlsx")