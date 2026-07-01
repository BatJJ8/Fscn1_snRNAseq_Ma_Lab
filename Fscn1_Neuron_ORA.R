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

# Load the integrated Seurat object
combined <- readRDS("combined_fascin_seurat_object.rds")

Idents(combined) <- combined$seurat_clusters

# Load the DEGs
wb_path <- "celltype_DGE_MAST.xlsx"
sheet_names <- getSheetNames(wb_path)

# Naming
abbrs <- c("Immature Oligodendrocytes"="IO",
           "Mature Oligodendrocytes"  ="MO",
           "Ependymal Cells"          ="Ependymal")

abbrs_rev <- setNames(names(abbrs), abbrs)

suffix_map <- c(
  " KO DOWN" = "markers_down",
  " KO UP"   = "markers_up",
  " DE"      = "markers"
)
suffixes_sorted <- names(suffix_map)[order(-nchar(names(suffix_map)))]

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

# Parse data
markers_list <- list()

for (sheet in sheet_names) {
  
  matched_suffix <- NULL
  for (suf in suffixes_sorted) {
    if (endsWith(sheet, suf)) {
      matched_suffix <- suf
      break
    }
  }
  
  if (is.null(matched_suffix)) {
    message("Skipping unrecognized sheet: ", sheet)
    next
  }
  
  short_ct <- substr(sheet, 1, nchar(sheet) - nchar(matched_suffix))
  
  ct <- if (short_ct %in% names(abbrs_rev)) abbrs_rev[[short_ct]] else short_ct
  
  slot_name <- suffix_map[[matched_suffix]]
  
  df <- read.xlsx(wb_path, sheet = sheet, rowNames = TRUE)
  
  if (is.null(markers_list[[ct]])) {
    markers_list[[ct]] <- list()
  }
  
  markers_list[[ct]][[slot_name]] <- df
}



# Filter neuron DEGs
neuron_up_genes <- rownames(markers_list[["Neurons"]][["markers_up"]])
neuron_down_genes <- rownames(markers_list[["Neurons"]][["markers_down"]])

neuronal_markers <- c(neuron_up_genes, neuron_down_genes)
cat("Number of significant neuronal mouse genes:", length(neuronal_markers), "\n")

# Map neuronal mouse gene symbols to human gene symbols using mapIt function
mapped_neurons <- mapIt(neuronal_markers, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
neuronal_human_genes <- unique(mapped_neurons$Homo_symbol)
neuronal_human_genes <- neuronal_human_genes[!is.na(neuronal_human_genes) & neuronal_human_genes != ""]
cat("Number of mapped human genes from neuronal DEGs:", length(neuronal_human_genes), "\n")

# Convert the human gene symbols to human Entrez IDs
neuronal_gene_ids <- bitr(neuronal_human_genes, 
                          fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)
if(is.null(neuronal_gene_ids) || nrow(neuronal_gene_ids) == 0){
  stop("No human Entrez IDs could be mapped from the neuronal DEGs.")
}
cat("Number of human Entrez IDs from neuronal DEGs:", nrow(neuronal_gene_ids), "\n")


# GO Enrichment Analysis (Neuron Upregulated DEGs)
neurons_up_BP <- enrichGO(gene          = neuron_up_genes,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

neurons_up_BP@result <- neurons_up_BP@result %>%
  mutate(
    gene_ratio_numeric = sapply(GeneRatio, function(x) {
      nums <- as.numeric(unlist(strsplit(x, "/")))
      nums[1] / nums[2]
    }),
    bg_ratio_numeric = sapply(BgRatio, function(x) {
      nums <- as.numeric(unlist(strsplit(x, "/")))
      nums[1] / nums[2]
    }),
    FoldEnrichment = gene_ratio_numeric / bg_ratio_numeric,
    q_log = -log10(qvalue),
    total_pathway_genes = sapply(BgRatio, function(x) {
      as.numeric(unlist(strsplit(x, "/")))[1]
    }),
    coverage = Count / total_pathway_genes
  )

neurons_up_BP_top_terms <- neurons_up_BP@result %>%
  arrange(p.adjust) %>%
  head(20)


# GO Enrichment Analysis (Neuron Downregulated DEGs)
neurons_down_BP <- enrichGO(gene        = neuron_down_genes,
                            OrgDb         = org.Mm.eg.db,
                            keyType       = "SYMBOL",  
                            ont           = "BP",     
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)     

neurons_down_BP@result <- neurons_down_BP@result %>%
  mutate(
    gene_ratio_numeric = sapply(GeneRatio, function(x) {
      nums <- as.numeric(unlist(strsplit(x, "/")))
      nums[1] / nums[2]
    }),
    bg_ratio_numeric = sapply(BgRatio, function(x) {
      nums <- as.numeric(unlist(strsplit(x, "/")))
      nums[1] / nums[2]
    }),
    FoldEnrichment = gene_ratio_numeric / bg_ratio_numeric,
    q_log = -log10(qvalue),
    total_pathway_genes = sapply(BgRatio, function(x) {
      as.numeric(unlist(strsplit(x, "/")))[1]
    }),
    coverage = Count / total_pathway_genes
  )

neurons_down_BP_top_terms <- neurons_down_BP@result %>% 
  arrange(p.adjust) %>% 
  head(20)




wb <- createWorkbook()

# GO BP - Upregulated
addWorksheet(wb, "Neuron GO BP UP")
writeData(wb, sheet = "Neuron GO BP UP",
          x     = neurons_up_BP_top_terms,
          rowNames = FALSE)

# GO BP - Downregulated
addWorksheet(wb, "Neuron GO BP DOWN")
writeData(wb, sheet = "Neuron GO BP DOWN",
          x     = neurons_down_BP_top_terms,
          rowNames = FALSE)

saveWorkbook(wb, "neuron_GO_BP_up_down.xlsx", overwrite = TRUE)
