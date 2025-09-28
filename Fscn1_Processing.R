# Code by Jason Huang with help from Zack Gaertner and Amanda Schneeweis
# Email: jasonhuang2024.2@u.northwestern.edu, zachary.gaertner@northwestern.edu, amanda.schneeweis@northwestern.edu

# Load libraries
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(ggbeeswarm)
library(knitr)
library(kableExtra)
library(gt)
library(patchwork)
library(tidyverse)
library(patchwork)
library(writexl)
library(ggrepel)
library(tools)
library(RColorBrewer)
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

theme_set(theme_minimal())

options(future.globals.maxSize = 1e9)  

set.seed(888888)

# Read in data
WT_F_data <- Read10X_h5("WT_F_filtered_feature_bc_matrix.h5")
KO_F_data <- Read10X_h5("KO_F_filtered_feature_bc_matrix.h5")
WT_M_data <- Read10X_h5("WT_M_filtered_feature_bc_matrix.h5")
KO_M_data <- Read10X_h5("KO_M_filtered_feature_bc_matrix.h5")

# Make Seurat objects
WT_F <- CreateSeuratObject(WT_F_data)
KO_F <- CreateSeuratObject(KO_F_data)
WT_M <- CreateSeuratObject(WT_M_data)
KO_M <- CreateSeuratObject(KO_M_data)

# Add metadata slots
WT_F[['orig.ident']] <- 'WT_F'
WT_F[['condition']] <- 'WT'
WT_F[['sex']] <- 'F'
WT_M[['orig.ident']] <- 'WT_M'
WT_M[['condition']] <- 'WT'
WT_M[['sex']] <- 'M'
KO_F[['orig.ident']] <- 'KO_F'
KO_F[['condition']] <- 'KO'
KO_F[['sex']] <- 'F'
KO_M[['orig.ident']] <- 'KO_M'
KO_M[['condition']] <- 'KO'
KO_M[['sex']] <- 'M'

WT_F <- PercentageFeatureSet(WT_F, pattern = "^mt-", col.name = "percent.mt")
WT_F <- PercentageFeatureSet(WT_F, pattern = "^Rp[sl]", col.name = "percent.ribo")
KO_F <- PercentageFeatureSet(KO_F, pattern = "^mt-", col.name = "percent.mt")
KO_F <- PercentageFeatureSet(KO_F, pattern = "^Rp[sl]", col.name = "percent.ribo")

WT_M <- PercentageFeatureSet(WT_M, pattern = "^mt-", col.name = "percent.mt")
WT_M <- PercentageFeatureSet(WT_M, pattern = "^Rp[sl]", col.name = "percent.ribo")
KO_M <- PercentageFeatureSet(KO_M, pattern = "^mt-", col.name = "percent.mt")
KO_M <- PercentageFeatureSet(KO_M, pattern = "^Rp[sl]", col.name = "percent.ribo")

# QC Filtering
VlnPlot(WT_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(WT_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(KO_F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)
VlnPlot(KO_M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

WT_F <- subset(WT_F, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 0.75 & percent.ribo < 0.75)
WT_M <- subset(WT_M, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 1 & percent.ribo < 0.5)
KO_F <- subset(KO_F, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 1.5 & percent.ribo < 0.75)
KO_M <- subset(KO_M, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 1 & percent.ribo < 0.5)

combined <- merge(WT_F, c(WT_M, KO_F, KO_M))

DefaultAssay(combined) <- "RNA"

# Normalization and integration
combined <- SCTransform(combined, vst.flavor = "v2", vars.to.regress = c('percent.mt', 'percent.ribo'))
combined <- RunPCA(combined)
combined <- IntegrateLayers(object = combined, method = CCAIntegration, orig.reduction = "pca", new.reduction = 'cca', assay = "SCT", normalization.method = "SCT")
combined <- IntegrateLayers(object = combined, method = RPCAIntegration, orig.reduction = "pca", new.reduction = 'rpca', assay = "SCT", normalization.method = "SCT")
combined <- FindNeighbors(combined, dims = 1:25, reduction = "cca")
combined <- FindClusters(combined, resolution = 0.2)
combined <- RunUMAP(combined, reduction = "cca", dims = 1:25, reduction.name = "umap.cca")

# Visualize clusters
DimPlot(combined, reduction = "umap.cca", label = TRUE)

DimPlot(combined, reduction = "umap.cca", split.by = "condition", label = TRUE) +
  labs(title = "UMAP by Condition")
DimPlot(combined, reduction = "umap.cca", split.by = "sex", label = TRUE) +
  labs(title = "UMAP by Sex")

# Distribution by label in each cluster
cluster_condition_df <- as.data.frame(table(Idents(combined), combined$condition))
colnames(cluster_condition_df) <- c("Cluster", "Condition", "CellCount")

cluster_sex_df <- as.data.frame(table(Idents(combined), combined$sex))
colnames(cluster_sex_df) <- c("Cluster", "Sex", "CellCount")

ggplot(cluster_condition_df, aes(x = Cluster, y = CellCount, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion", title = "Proportion of WT and KO in Each Cluster") +
  theme_minimal()

ggplot(cluster_sex_df, aes(x = Cluster, y = CellCount, fill = Sex)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion", title = "Proportion of M and F in Each Cluster") +
  theme_minimal()

# Prepare for marker detection
combined <- PrepSCTFindMarkers(combined)

# Finding all markers
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
markers

# Marker list for calculating module score
marker_list <- list(
  Neurons = c("Rbfox3", "Elavl2", "Elavl3", "Elavl4", "Syn1", "Syp", "Dcx",
              "Slc17a6", "Slc17a7", "Slc32a1", "Gad1", "Gad2", "Pvalb", "Sst", 
              "Vip", "Lhx6"),
  Microglia = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Itgam", "Hexb", "C1qa", "Csf1r", "Siglech"),
  OPCs = c("Pdgfra", "Cspg4", "Sox10", "Olig1", "Olig2", "Nkx2-2", "Tnr"),
  Oligodendrocytes = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Sox10", "Ermn", "Olig1", "Olig2", "Tnr", "Enpp6", "Gpr17"),
  Astrocytes = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "Gja1", "Itpr2", "Sparcl1", "Glul", "Sox9"),
  Fibroblasts = c("Dcn", "Lum", "Col1a1", "Col1a2", "Prrx1", "Thy1", "Pi16", "Pdgfrb"),
  Endothelial = c("Pecam1", "Cdh5", "Vwf", "Plvap", "Flt1", "Kdr", "Cldn5", "Tek", "Nos3")
)

# Calculate module scores
combined <- AddModuleScore(object = combined, features = marker_list, name = "CellTypeScore", ctrl = 50)

# Visualize module scores
FeaturePlot(combined, features = c("CellTypeScore1", "CellTypeScore2", "CellTypeScore3", 
                                   "CellTypeScore4", "CellTypeScore5", "CellTypeScore6",
                                   "CellTypeScore7"),
            cols = c("lightgrey", "blue"))

VlnPlot(combined, features = c("CellTypeScore1", "CellTypeScore2", "CellTypeScore3", 
                               "CellTypeScore4", "CellTypeScore5", "CellTypeScore6",
                               "CellTypeScore7"),
        group.by = "seurat_clusters", pt.size = 0)

# Compute per-cell Z-scores for each cell type
meta_data <- combined@meta.data

normalized_meta_data <- meta_data %>%
  mutate(across(starts_with("CellTypeScore"), ~ scale(.)))

score_summary <- normalized_meta_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Neurons = median(CellTypeScore1),
    Microglia = median(CellTypeScore2),
    OPCs = median(CellTypeScore3),
    Oligodendrocytes = median(CellTypeScore4),
    Astrocytes = median(CellTypeScore5),
    Fibroblasts = median(CellTypeScore6),
    Endothelial = median(CellTypeScore7)
  )

print(score_summary)

# Assign the cell type for each cluster.
score_summary <- score_summary %>%
  rowwise() %>%
  mutate(AssignedType = {
    scores <- c(Neurons, Microglia, OPCs, Oligodendrocytes, Astrocytes, Fibroblasts, Endothelial)
    # If the highest score is below 0.2, assign "Unassigned"
    if(max(scores, na.rm = TRUE) < 0.2) {
      "Unassigned"
    } else {
      cell_types <- c("Neurons", "Microglia", "OPCs", "Oligodendrocytes", "Astrocytes", "Fibroblasts", "Endothelial")
      cell_types[which.max(scores)]
    }
  }) %>%
  ungroup()

print(score_summary, n = 100)

# Create a lookup vector mapping cluster IDs to the assigned cell type.
cluster_to_type <- setNames(score_summary$AssignedType, as.character(score_summary$seurat_clusters))
print(cluster_to_type)

combined@meta.data$AssignedCellType <- cluster_to_type[as.character(combined@meta.data$seurat_clusters)]

# Check celltypes and manually assign unknown clusters 
Idents(combined) <- combined$seurat_clusters
DimPlot(combined, label = TRUE)

FindMarkers(combined, ident.1 = 0, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 1, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 2, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 3, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 4, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 5, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 6, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 7, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 8, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 9, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 10, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 11, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 12, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 13, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 14, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FindMarkers(combined, ident.1 = 15, test.use = "roc", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "0")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "1")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "2")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "3")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "4")] <- "Oligodendrocytes"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "5")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "6")] <- "Astrocytes"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "7")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "8")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "9")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "10")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "11")] <- "Microglia"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "12")] <- "OPCs"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "13")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "14")] <- "Neurons"
combined@meta.data$AssignedCellType[which(combined@meta.data$seurat_clusters == "15")] <- "Pericytes"

# Assign more detailed annotations if possible
combined@meta.data$SpecificCellType <- combined@meta.data$AssignedCellType
if("SpecificCellType" %in% colnames(combined@meta.data)){
  combined@meta.data$SpecificCellType <- as.character(combined@meta.data$SpecificCellType)
} else {
  combined@meta.data$SpecificCellType <- as.character(combined@meta.data$AssignedCellType)
}

combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "0")] <- "Cux2^{'+'}~L2/L3/L4~Excitatory~Neurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "1")] <- "Rorb^{'+'}~L4~Excitatory~Neurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "2")] <- "Foxp2^{'+'}~L5/L6~Excitatory~Neurons" 
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "3")] <- "Il1rapl2^{'+'}/Satb2^{'+'}~Excitatory~Neurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "5")] <- "Nxph1^{'+'}/Sox6^{'+'}~MGE-derived~Inhibitory~Interneurons" 
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "7")] <- "Adarb2^{'+'}~CGE-derived~Inhibitory~Interneurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "8")] <- "Tafa1^{'+'}/Tox^{'+'}~Excitatory~Neurons" 
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "9")] <- "Rgs9^{'+'}/Rarb^{'+'}~MSNs"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "10")] <- "Tshz2^{'+'}/Vwc2l^{'+'}~Excitatory~Neurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "13")] <- "Col11a1^{'+'}/Tmem163^{'+'}/Ntng2^{'+'}~Excitatory~Neurons"
combined@meta.data$SpecificCellType[which(combined@meta.data$seurat_clusters == "14")] <- "Tle4^{'+'}/Inpp4b^{'+'}~L6~Excitatory~Neurons"

combined@meta.data$SpecificCellType <- factor(combined@meta.data$SpecificCellType)

# Generate UMAP
DimPlot(combined, group.by = "SpecificCellType") + scale_color_discrete(labels = function(x) parse(text = x))
p <- DimPlot(combined, group.by = "SpecificCellType") + theme_void() + NoLegend()
p
ggsave(filename = "fascin_umap_celltype.pdf", plot = p, device = "pdf", width = 8, height = 6)

# Save processed Seurat object
saveRDS(combined, file = "combined_fascin_seurat_object.rds")

