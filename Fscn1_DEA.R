# Code by Jason Huang modified from code provided by Zack Gaertner and Amanda Schneeweis

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

# Load comparison lists
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

conversion <- function(gene_list) {
  sapply(gene_list, function(gene) {
    gsub("^(\\w)(\\w*)$", "\\U\\1\\L\\2", gene, perl = TRUE)
  })
}

inflam <- read_excel("C:/Users/weihu/OneDrive/Desktop/Ma Lab/Important Lists/jason_comprehensive_inflammation_list.xlsx")
schwaner_inflam <- inflam$`Inflam related genes (Schwaner)` %>% 
  conversion()
cdg_inflam <- inflam$`Human inflam genes (CDG)` %>% 
  conversion()
cytokines <- inflam$`Cytokines (KEGG)` %>% 
  conversion()
necroptosis <- inflam$`Necroptosis (KEGG)` %>% 
  conversion()
pyroptosis <- inflam$`Pyroptosis (MGI)` %>% 
  conversion()
combined_inflam <- c(schwaner_inflam, cdg_inflam, cytokines, necroptosis, pyroptosis)

mitocarta <- read_excel("C:/Users/weihu/OneDrive/Desktop/Ma Lab/Important Lists/MitoCarta3.0_mitochondrial_genes.xlsx")
mitocarta <- mitocarta$Symbol

AD_KEGG <- read_excel("C:/Users/weihu/OneDrive/Desktop/Ma Lab/Important Lists/AD_KEGG_gene_list.xlsx")
AD_KEGG_genes <- AD_KEGG$gene_symbol %>% 
  conversion()
AD_MalaCards <- read_csv("C:/Users/weihu/OneDrive/Desktop/Ma Lab/Important Lists/AD_MalaCards_gene_list.csv")
AD_MalaCards_genes <- AD_MalaCards$`Gene Symbol` %>% 
  conversion()
combined_AD <- c(AD_KEGG_genes, AD_MalaCards_genes)


# Load the integrated Seurat object
combined <- readRDS("combined_fascin_seurat_object.rds")

# Run MAST to generate DEGs across all cells
Idents(combined) <- "condition"

mast_res <- FindMarkers(
  object        = combined,
  ident.1       = "KO",
  ident.2       = "WT",
  test.use      = "MAST",
  latent.vars   = "sex",
  min.pct       = 0.10,
  logfc.threshold = 0.0
)

# Filter out mitochondrial & ribosomal genes
mast_res <- mast_res[
  !grepl("^mt-",    rownames(mast_res),      ignore.case = TRUE) &
    !grepl("^(Rps|Rpl)", rownames(mast_res),    ignore.case = TRUE),
]

# Generate annotations for DEGs
ann <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = rownames(mast_res),
  columns = "GENENAME",
  keytype = "SYMBOL"
)
rownames(ann) <- ann$SYMBOL
mast_res$annotation <- ann[rownames(mast_res),"GENENAME"]

mast_res <- mast_res[order(mast_res$avg_log2FC, decreasing = TRUE), ]
mast_up   <- subset(mast_res, avg_log2FC >  0 & p_val < 0.01)
mast_down <- subset(mast_res, avg_log2FC <  0 & p_val < 0.01)

wb <- createWorkbook()
addWorksheet(wb, "MAST_All_DEGs");  writeData(wb, "MAST_All_DEGs",  mast_res,    rowNames = TRUE)
addWorksheet(wb, "MAST_KO_UP");    writeData(wb, "MAST_KO_UP",    mast_up,     rowNames = TRUE)
addWorksheet(wb, "MAST_KO_DOWN");  writeData(wb, "MAST_KO_DOWN",  mast_down,   rowNames = TRUE)
saveWorkbook(wb, "KO_vs_WT_DGE_MAST.xlsx", overwrite = TRUE)

# Create a composite identifier combining condition & assigned cell type
DefaultAssay(combined) <- "SCT"
combined$composite_ident <- paste(combined$condition, tolower(combined$AssignedCellType), sep = " ")
Idents(combined) <- combined$composite_ident

# Identify composite IDs of the form 
cell_types      <- sort(unique(combined$AssignedCellType))
markers_list    <- list()

# Run MAST per celltype
for (ct in cell_types) {
  
  ident_wt <- paste("WT", tolower(ct))
  ident_ko <- paste("KO", tolower(ct))
  
  if (ident_wt %in% combined$composite_ident &&
      ident_ko %in% combined$composite_ident) {
    # Seurat's built-in MAST function with sex as latent variable
    de <- FindMarkers(
      combined,
      ident.1        = ident_ko,
      ident.2        = ident_wt,
      test.use       = "MAST",
      latent.vars    = c("sex"),   
      min.pct        = 0.10,
      logfc.threshold= 0,
      verbose        = FALSE
    )
    
    # Rename percentage columns
    colnames(de)[colnames(de) == "pct.1"] <- paste0("pct.1 (KO)")
    colnames(de)[colnames(de) == "pct.2"] <- paste0("pct.2 (WT)")
    
    # Filter out mitochondrial & ribosomal genes
    de <- de[
      !grepl("^mt-",    rownames(de),      ignore.case = TRUE) &
        !grepl("^(Rps|Rpl)", rownames(de),    ignore.case = TRUE),
    ]
    
    # Generation annotations for DEGs
    ann <- AnnotationDbi::select(
      org.Mm.eg.db,
      keys    = rownames(de),
      columns = "GENENAME",
      keytype = "SYMBOL")
    rownames(ann)     <- ann$SYMBOL
    de$annotation     <- ann[rownames(de), "GENENAME"]
    
    # Subsetting
    de <- de[order(de$avg_log2FC, decreasing = TRUE), ]
    
    de_up        <- subset(de, avg_log2FC >  0 & p_val < 0.01)
    de_down      <- subset(de, avg_log2FC <  0 & p_val < 0.01) %>%
      arrange(avg_log2FC)
    
    inflam_tbl   <- subset(de, rownames(de) %in% combined_inflam & p_val < 0.01) %>%
      arrange(avg_log2FC)
    mitocarta_tbl<- subset(de, rownames(de) %in% mitocarta       & p_val < 0.01) %>%
      arrange(desc(avg_log2FC))
    ad_tbl       <- subset(de, rownames(de) %in% combined_AD     & p_val < 0.01) %>%
      arrange(desc(avg_log2FC))
    
    markers_list[[ct]] <- list(
      markers      = de,
      markers_up   = de_up,
      markers_down = de_down,
      inflam       = inflam_tbl,
      mitocarta    = mitocarta_tbl,
      AD           = ad_tbl
    )
    
  } else {
    message("Skipping ", ct, " – missing WT or KO cells.")
  }
}

abbrs <- c("Immature Oligodendrocytes"="IO",
           "Mature Oligodendrocytes"  ="MO",
           "Astrocyte Precursor Cells"="APCs",
           "Ependymal Cells"          ="Ependymal")

wb <- createWorkbook()

for (ct in names(markers_list)) {
  short_ct  <- if (ct %in% names(abbrs)) abbrs[[ct]] else ct
  
  # All DE genes
  sheet_all  <- paste(short_ct, "DE")
  addWorksheet(wb, sheet_all)
  writeData(wb, sheet = sheet_all,
            x     = markers_list[[ct]]$markers,
            rowNames = TRUE)
  
  # KO-up
  sheet_up   <- paste(short_ct, "KO UP")
  addWorksheet(wb, sheet_up)
  writeData(wb, sheet = sheet_up,
            x     = markers_list[[ct]]$markers_up,
            rowNames = TRUE)
  
  # KO-down
  sheet_down <- paste(short_ct, "KO DOWN")
  addWorksheet(wb, sheet_down)
  writeData(wb, sheet = sheet_down,
            x     = markers_list[[ct]]$markers_down,
            rowNames = TRUE)
  
  # Inflammation panel
  sheet_infl <- paste(short_ct, "Inflam")
  addWorksheet(wb, sheet_infl)
  writeData(wb, sheet = sheet_infl,
            x     = markers_list[[ct]]$inflam,
            rowNames = TRUE)
  
  # Mitocarta panel
  sheet_mito <- paste(short_ct, "Mitocarta")
  addWorksheet(wb, sheet_mito)
  writeData(wb, sheet = sheet_mito,
            x     = markers_list[[ct]]$mitocarta,
            rowNames = TRUE)
  
  # AD panel
  sheet_ad   <- paste(short_ct, "AD")
  addWorksheet(wb, sheet_ad)
  writeData(wb, sheet = sheet_ad,
            x     = markers_list[[ct]]$AD,
            rowNames = TRUE)
}

saveWorkbook(wb, "celltype_DGE_MAST.xlsx", overwrite = TRUE)