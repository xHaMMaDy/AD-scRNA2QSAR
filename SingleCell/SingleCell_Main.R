#=========================================================================================
# 🎯 Single-Cell RNA-seq Analysis Pipeline using Seurat
#=========================================================================================
# This script performs preprocessing, normalization, dimensional reduction, clustering,
# doublet removal, visualization, marker gene identification, and cluster annotation.
# Authors: Ibrahim Hammad - Mira Moheb - Abdulrahman Wagih - Reem Sharaf - Lorance Gergis
# Created: 12-03-2025
# Last Edited: 14-04-2025
#=========================================================================================
#=========================================================================================
# 📌 Load Required Libraries
#=========================================================================================
library(DoubletFinder)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(scater)
library(scran)
library(harmony)
library(batchelor)
library(SeuratWrappers)
library(glmGamPoi)
library(SingleR)
library(BiocParallel)
cat("[✔] All required libraries are  loaded.\n")
#=========================================================================================
# 📁 Define All Paths (Project structure)
#=========================================================================================
script_dir <- "C:/Users/HaMMaDy/Desktop/Grad"
input_rds_path <- file.path(script_dir, "All_Studies.rds")
qc_dir <- file.path(script_dir, "results", "QC")
int_dir <- file.path(script_dir, "results", "Integration")
doublet_input_dir <- file.path(script_dir, "results", "Ready4Doublet")
doublet_output_dir <- file.path(script_dir, "results", "DoubletDone")
doublet_plot_dir <- file.path(script_dir, "results", "DoubletResult")
integration_save_path <- file.path(script_dir, "results", "All_Studies_Harmony.rds")
final_save_path <- file.path(script_dir, "results", "All_Studies_Done.rds")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(int_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doublet_input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doublet_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doublet_plot_dir, showWarnings = FALSE, recursive = TRUE)
#=========================================================================================
# 🔹 1. Load Data
#=========================================================================================
cat("[1] Loading Seurat object...","\n")
start_time <- Sys.time()
seurat_obj <- readRDS(input_rds_path)
seurat_obj <- JoinLayers(seurat_obj)
end_time <- Sys.time()
cat("[✔]  Loaded in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 2. Quality Control & Cell Filtering
#=========================================================================================
cat("[2] Running quality control...","\n")

if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}
# Save QC violin and scatter plots
vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(qc_dir, "Violin_QC.png"), plot = vln_plot, width = 24, height = 12, dpi = 200)
scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(file.path(qc_dir, "Scatter_QC.png"), plot = scatter1 + scatter2, width = 24, height = 12, dpi = 200)
# Cell filtering
start_time <- Sys.time()
seurat_obj <- subset(seurat_obj, subset = (
  (orig.ident %in% c("Entorhinal Cortex", "Hippocampus", "Prefrontal cortex", "Superior frontal cortex")) &
    nFeature_RNA > 200 &
    ((orig.ident == "Entorhinal Cortex" & nFeature_RNA < 1600) |
       (orig.ident == "Hippocampus" & nFeature_RNA < 10000) |
       (orig.ident == "Prefrontal cortex" & nFeature_RNA < 10000) |
       (orig.ident == "Superior frontal cortex" & nFeature_RNA < 4500)) &
    percent.mt < 5
))
end_time <- Sys.time()
cat("[✔]  QC filtering completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 3. Data Normalization
#=========================================================================================
cat("[3] Normalizing data...","\n")
start_time <- Sys.time()
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
end_time <- Sys.time()
cat("[✔]  Normalized in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 4. Identify Highly Variable Features
#=========================================================================================
cat("[4] Identifying variable features...","\n")
start_time <- Sys.time()
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(file.path(qc_dir, "Top10.png"), plot = plot2, bg = 'white', width = 24, height = 12, dpi = 200)
end_time <- Sys.time()
cat("[✔]  Variable features identified in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 5. Data Scaling
#=========================================================================================
cat("[5] Scaling data...","\n")
start_time <- Sys.time()
seurat_obj <- ScaleData(seurat_obj)
end_time <- Sys.time()
cat("[✔]  Data scaled in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 6. Dimensional Reduction (PCA)
#=========================================================================================
cat("[6] Running PCA...","\n")
start_time <- Sys.time()
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
end_time <- Sys.time()
cat("[✔]  PCA completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(qc_dir, "ElbowPlot.png"), plot = elbow_plot, bg = 'white', width = 12, height = 6, dpi = 200)
num_pcs <- 35
#=========================================================================================
# 🔹 7. Clustering
#=========================================================================================
cat("[7] Finding neighbors and clustering...","\n")
start_time <- Sys.time()
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:num_pcs)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
end_time <- Sys.time()
cat("[✔]  Clustering completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 8. UMAP for Visualization
#=========================================================================================
cat("[8] Running UMAP...","\n")
start_time <- Sys.time()
seurat_obj <- RunUMAP(seurat_obj, dims = 1:num_pcs)
end_time <- Sys.time()
cat("[✔]  UMAP completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 9. Visualization of Clusters
#=========================================================================================
cat("[9] Generating UMAP plots by metadata groupings...","\n")
start_time <- Sys.time()
umap_plot1 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'orig.ident', raster = FALSE)
umap_plot2 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'condition', raster = FALSE)
umap_plot3 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'Tissue_Condition', raster = FALSE)
umap_combined <- grid.arrange(umap_plot1, umap_plot2, umap_plot3, ncol = 2, nrow = 2)
ggsave(file.path(qc_dir, "UMAP_Combined.png"), plot = umap_combined, bg = 'white', width = 24, height = 18, dpi = 200)
ggsave(file.path(qc_dir, "UMAP_Combined.svg"), plot = umap_combined, bg = 'white', width = 24, height = 18, dpi = 200)
end_time <- Sys.time()
cat("[✔]  UMAP visualization completed in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")
#=========================================================================================
# 🔹 10. Doublet Detection and Filtering
#=========================================================================================
cat("[10] Doublet Detection and Filtering...\n")

# Set the identity class based on original condition
metadata <- seurat_obj@meta.data
seurat_obj <- SetIdent(seurat_obj, value = metadata$orig.ident)
standard_rate <- 0.075
num_pcs <- 35

sweep_list <- paramSweep(seurat_obj, PCs = 1:num_pcs, sct = FALSE)
sweep_stats <- summarizeSweep(sweep_list, GT = FALSE)
find.pK(sweep_stats)

homotypic_prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)
nExp_poi <- round(standard_rate * nrow(seurat_obj@meta.data), 0)

pN_val <- 0.25
pK_val <- 0.09

seurat_obj <- doubletFinder(seurat_obj, PCs = 1:num_pcs, pN = pN_val, pK = pK_val,
                            nExp = nExp_poi, sct = FALSE)

col_poi <- paste0("DF.classifications_", pN_val, "_", pK_val, "_", nExp_poi)
seurat_obj@meta.data$DoubletStatus <- seurat_obj@meta.data[[col_poi]]

singlet_count <- sum(seurat_obj@meta.data$DoubletStatus == "Singlet")
doublet_count <- sum(seurat_obj@meta.data$DoubletStatus == "Doublet")

cat("[✔] Singlets: ", singlet_count, " | Doublets: ", doublet_count, "\n")

plot_title <- paste("Singlets =", singlet_count, "| Doublets =", doublet_count)
plot <- DimPlot(seurat_obj, reduction = "umap", group.by = col_poi) +
  scale_color_manual(values = c("Singlet" = "black", "Doublet" = "red")) +
  ggtitle(plot_title) +
  theme_minimal()

ggsave("Doublet_UMAP.png", plot, width = 12, height = 6, bg = "white")

seurat_obj <- subset(seurat_obj, subset = DoubletStatus == "Singlet")

cat("[✔] Doublets removed. Only Singlets are retained.\n")
#=========================================================================================
# 🔹 11. Integration Methods Comparison: PCA (Unintegrated), Harmony, FastMNN.
#=========================================================================================
# Split by condition
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$condition)

# Basic processing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_sct <- SCTransform(seurat_obj)

# PCA & Clustering (Unintegrated baseline)
seurat_sct <- RunPCA(seurat_sct)
seurat_sct <- FindNeighbors(seurat_sct, dims = 1:num_pcs, reduction = "pca")
seurat_sct <- FindClusters(seurat_sct, resolution = 0.4, cluster.name = "unintegrated_clusters")
seurat_sct <- RunUMAP(seurat_sct, dims = 1:num_pcs, reduction = "pca", reduction.name = "umap.unintegrated")

p_unintegrated <- DimPlot(seurat_sct, reduction = "umap.unintegrated", group.by = 'condition')
p_unintegrated
ggsave(file.path(int_dir, "p_unintegrated.png"), plot = p_unintegrated, bg = 'white', width = 24, height = 12, dpi = 200)
#=========================================================================================
# 🔹 11.1 Harmony Integration
#=========================================================================================
cat("[11.1] Running Harmony integration...\n")
seurat_harmony <- seurat_sct
start_time <- Sys.time()
#Harmony
seurat_harmony <- RunHarmony(seurat_harmony, group.by.vars = 'condition', plot_convergence = FALSE)

# Harmony reduction for UMAP & clustering
seurat_harmony <- seurat_harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:num_pcs, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:num_pcs) %>%
  FindClusters(resolution = 0.4)

p_harmony <- DimPlot(seurat_harmony, reduction = "umap.harmony", group.by = 'condition')
p_harmony
ggsave(file.path(int_dir, "p_harmony.png"), plot = p_harmony, bg = 'white', width = 24, height = 12, dpi = 200)
saveRDS(seurat_harmony, integration_save_path)
end_time <- Sys.time()
cat("[ ! ] Haromony integration completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
#=========================================================================================
# 🔹 11.2 FastMNN Integration
#=========================================================================================
cat("[11.2] Running FastMNN integration...\n")

start_time <- Sys.time()
# Split the SCT-processed object into a list by condition
seurat_list <- SplitObject(seurat_sct, split.by = "condition")
# FastMNN requires the same preprocessing for each split object
# So we apply SCTransform again per object (if not already done)
seurat_list <- lapply(seurat_list, SCTransform)

# Run FastMNN integration
seurat_mnn <- RunFastMNN(object.list = seurat_list)

# Run clustering and UMAP based on FastMNN
seurat_mnn <- RunUMAP(seurat_mnn, reduction = "mnn", dims = 1:20, reduction.name = "umap.mnn")
seurat_mnn <- FindNeighbors(seurat_mnn, reduction = "mnn", dims = 1:20)
seurat_mnn <- FindClusters(seurat_mnn, resolution = 0.4)

p_mnn <- DimPlot(seurat_mnn, reduction = "umap.mnn", group.by = 'condition')
p_mnn
ggsave(file.path(int_dir, "p_mnn.png"), plot = p_mnn, bg = 'white', width = 24, height = 12, dpi = 200)
saveRDS(seurat_mnn, integration_save_path)
end_time <- Sys.time()
cat("[ ! ] FastMNN integration completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
#=========================================================================================
# 🔹 11.3 Combined UMAPs Comparison
#=========================================================================================
# Combine all UMAPs (Unintegrated + Harmony + FastMNN)
p_unintegrated <- p_unintegrated + ggtitle("Unintegrated (PCA)")
p_harmony      <- p_harmony      + ggtitle("Harmony Integration")
p_mnn          <- p_mnn          + ggtitle("FastMNN Integration")
combined_umap_plot <- (p_unintegrated | p_harmony | p_mnn)
combined_umap_plot
ggsave(file.path(int_dir, "integration_umap_comparison_all.png"), plot = combined_umap_plot, bg = 'white', width = 36, height = 18, dpi = 200)
#=========================================================================================
# 🔹 12. Marker Gene Identification and DotPlot Visualization
#=========================================================================================
cat("[12] Identifying marker genes...\n")
start_time <- Sys.time()
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1.5)

top_marker <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 1, with_ties = FALSE) %>%
  ungroup()

top_marker_genes <- top_marker %>%
  arrange(cluster, p_val_adj) %>%
  pull(gene) %>%
  unique()

full_dotplot <- DotPlot(seurat_obj, features = top_marker_genes, group.by = "seurat_clusters") +
  labs(title = "Top Marker per Cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

print(full_dotplot)

ggsave("TopMarker_DotPlot_AllClusters.png", plot = full_dotplot, width = 16, height = 10, dpi = 200)

end_time <- Sys.time()
cat("[✔] Marker identification and DotPlot saved in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
#=========================================================================================
# 🔹 12.1 Cell Type Annotation using SingleR with SEA-AD Reference
#=========================================================================================
cat("[12.1] Annotating cell types with SingleR using SEA-AD reference...\n")
start_time <- Sys.time()

# Load reference object
ref_seurat <- readRDS("C:/Users/HaMMaDy/Desktop/Grad/AnnotationReference/Reference_RNAseq_all-nuclei.2022-06-07.rds")

# Build SingleCellExperiment for reference
ref_expr <- GetAssayData(ref_seurat, layer = "data")
ref_labels <- ref_seurat$subclass_label
ref_sce <- SingleCellExperiment(
  assays = list(logcounts = ref_expr),
  colData = DataFrame(label.main = ref_labels)
)

# Prepare dataset
seurat_expr <- GetAssayData(seurat_obj, layer = "data")
seurat_sce <- SingleCellExperiment(
  assays = list(logcounts = seurat_expr),
  colData = seurat_obj@meta.data
)

# SingleR annotation
singleR_results <- SingleR(
  test = seurat_sce,
  ref = ref_sce,
  labels = ref_sce$label.main,
)

#Add SingleR labels to Seurat object
seurat_obj$SingleR_Labels <- singleR_results$labels

#UMAP
p_singleR <- DimPlot(
  seurat_obj,
  group.by = "SingleR_Labels",
  label = TRUE,
  repel = TRUE,
  label.size = 8
) +
  ggtitle("SingleR Cell Type Annotation") +
  theme_minimal()

ggsave(filename = "results/QC/SingleR_UMAP_SEAAD.png", plot = p_singleR, width = 24, height = 12, dpi = 200)

write.csv(data.frame(
  Cell = colnames(seurat_obj),
  CellType = seurat_obj$SingleR_Labels
), "results/SingleR_CellTypeAnnotations_SEAAD.csv", row.names = FALSE)

#Done
end_time <- Sys.time()
cat("[✔] SingleR annotation completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
#=========================================================================================
# 🔹 13. Save Final Processed Object
#=========================================================================================
cat("[13] Saving final Seurat object...","\n")
start_time <- Sys.time()
saveRDS(singleR_results,"C:/Users/HaMMaDy/Documents/Backup/singleR_results_Alldata.rds")
saveRDS(markers, "C:/Users/HaMMaDy/Documents/Backup/AllStudies_Markers.rds")
saveRDS(singleR_results, "C:/Users/HaMMaDy/Documents/Backup/SingleR_Results.rds")
saveRDS(seurat_obj, "C:/Users/HaMMaDy/Desktop/Grad/AllStudies_SeuratDone_Annotated.rds")
end_time <- Sys.time()
cat("[✔] Final object saved in ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds.","\n")