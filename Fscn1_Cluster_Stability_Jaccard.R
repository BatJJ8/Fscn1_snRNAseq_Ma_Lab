# Load libraries
library(Seurat)
library(igraph)
library(ggplot2)
library(reshape2)

theme_set(theme_minimal())

set.seed(888888)

# Load the integrated Seurat object
combined <- readRDS("combined_fascin_seurat_object.rds")

Idents(combined) <- combined$seurat_clusters

# Original clustering settings
reduct_use <- "cca"
dims_use   <- 1:25
res_use    <- 0.2

# Bootstrap settings
prop_keep  <- 0.80
n_iter     <- 50

# Rerun and extract SNN
combined <- FindNeighbors(combined, reduction = reduct_use, dims = dims_use, verbose = FALSE)

snn_name <- if ("orig_snn" %in% names(combined@graphs)) "orig_snn" else "SCT_snn"
full_snn <- combined@graphs[[snn_name]]  # sparse adjacency matrix
stopifnot(inherits(full_snn, "dgCMatrix"))

# Cell index mapping
all_cells <- colnames(combined)
cell_index <- setNames(seq_along(all_cells), all_cells)  # map cellname -> column index

# igraph Louvain clustering
cluster_on_snn <- function(snn_sub) {
  # Convert SNN subgraph to igraph object
  # Treat matrix as weighted undirected graph
  g <- graph_from_adjacency_matrix(snn_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
  memb <- cluster_louvain(g, resolution = res_use)$membership
  factor(memb)
}

# Bootstrap stability (using Jaccard index as metric)
jac_list <- vector("list", n_iter)

for (b in seq_len(n_iter)) {
  # 1. Subsample cells
  keep <- sample(all_cells, size = floor(prop_keep * length(all_cells)))
  keep_idx <- cell_index[keep]
  
  # 2. Extract induced SNN subgraph for kept cells
  sub_snn <- full_snn[keep_idx, keep_idx, drop = FALSE]
  
  # 3. Re-cluster using Louvain on kept cells
  memb <- cluster_on_snn(sub_snn)
  names(memb) <- keep   # align membership vector to cell names
  
  # 4. Calculate Jaccard index between original and new clusters
  orig_lab <- Idents(combined)[keep]
  sub_lab  <- memb
  
  ## Contingency table 
  tab   <- table(orig_lab, sub_lab)
  rowN  <- rowSums(tab)
  colN <- colSums(tab)
  
  ## Jaccard calculation (|A ∩ B| / |A ∪ B|)
  union <- outer(rowN, colN, "+") - tab
  J     <- tab / union; J[is.na(J)] <- 0
  
  ## Best matching subcluster for each original cluster
  jac_list[[b]] <- apply(J, 1, max)
  
  rm(sub_snn)
  gc()
}

# Visualize as box and whisker plot
jac_mat <- do.call(rbind, jac_list)
df <- reshape2::melt(jac_mat, varnames = c("Iter","Cluster"), value.name = "Jaccard")
df$Cluster <- factor(df$Cluster, levels = levels(Idents(combined)))


theme_set(theme_bw(base_size = 14))
p <- ggplot(df, aes(Cluster, Jaccard)) +
  geom_boxplot(fill = "grey80", color = "black", 
               outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.6) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "Cluster Stability",
       y = "Normalized Jaccard Index", 
       x = "Cluster Number") + 
  theme(panel.grid = element_blank())
p

# Save plot
ggsave(filename = "jaccard_box_plot.pdf", plot = p, device = "pdf", width = 8, height = 6)