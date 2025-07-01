# =======================
# [1] Load Required Libraries
# =======================
library(Seurat)
library(SingleCellExperiment)
library(miloR)
library(scater)
library(ggplot2)
library(dplyr)
library(patchwork)
library(BiocParallel)

# =======================
# [2] Load Seurat Object
# =======================
cat("[1] Loading Seurat object...\n")
seurat_obj <- readRDS("C:/Users/Hammady/Documents/Backup/AllStudies_Annotated.rds")
Reductions(seurat_obj)

# =======================
# [3] Prepare SingleCellExperiment with SCT assay
# =======================
cat("[2] Converting to SingleCellExperiment...\n")
sce <- SingleCellExperiment(
  assays = list(SCT = GetAssayData(seurat_obj, assay = "SCT", layer = "data")),
  colData = seurat_obj@meta.data
)

# Add UMAP from custom reduction (e.g. umap.harmony)
cat("[3] Adding UMAP coordinates...\n")
reducedDim(sce, "UMAP") <- Embeddings(seurat_obj, "umap.harmony")

cat("[4] Using Harmony PCs from Seurat...\n")
num_pcs <- 35 
reducedDim(sce, "PCA") <- Embeddings(seurat_obj, "harmony")[, 1:num_pcs]

# =======================
# [4] Build Milo Object and Graph
# =======================
cat("[5] Building Milo object and KNN graph...\n")
milo_obj <- Milo(sce)
milo_obj <- buildGraph(milo_obj, k = 30, d = num_pcs)
saveRDS(milo_obj,"Milo_Step5Done.rds")

# =======================
# [5] Define Neighborhoods
# =======================
cat("[6] Defining neighborhoods...\n")
start_time <- Sys.time()
milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = num_pcs, refined = TRUE)
end_time <- Sys.time()
cat("Neighborhoods built in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")

saveRDS(milo_obj, "C:/Users/Hammady/Documents/Backup/AllStudies_milo_obj_MakeNhoods_18_04.rds")
milo_obj <- readRDS("C:/Users/Hammady/Documents/Backup/AllStudies_milo_obj_MakeNhoods_18_04.rds")

# =======================
# [6] Plot Neighborhood Sizes
# =======================
cat("[7] Plotting neighborhood size distribution...\n")
milosizeplot <- plotNhoodSizeHist(milo_obj)
ggsave("results/Milo/Nhood.png", plot = milosizeplot, width = 16, height = 10, dpi = 300)
ggsave("results/Milo/Nhood.svg", plot = milosizeplot, width = 16, height = 10, dpi = 300)

# =======================
# [7] Fix Sample Labels
# =======================
cat("[8] Resolving NA sample labels...\n")
cell_names <- colnames(milo_obj)
sample_ids <- milo_obj$sample
na_cells <- is.na(sample_ids)

sample_ids[na_cells & grepl("_AD7_AD8", cell_names)] <- "AD_ERC_1"
sample_ids[na_cells & grepl("_AD5_AD6", cell_names)] <- "AD_ERC_2"
sample_ids[na_cells & grepl("_AD3_AD4", cell_names)] <- "AD_ERC_3"
sample_ids[na_cells & grepl("_AD1_AD2", cell_names)] <- "AD_ERC_4"
sample_ids[na_cells & grepl("_Ct5_Ct6", cell_names)] <- "CT_ERC_1"
sample_ids[na_cells & grepl("_Ct3_Ct4", cell_names)] <- "CT_ERC_2"
sample_ids[na_cells & grepl("_Ct1_Ct2", cell_names)] <- "CT_ERC_3"
sample_ids[na_cells & grepl("_Ct7_Ct8", cell_names)] <- "CT_ERC_4"

milo_obj$sample <- sample_ids
cat("Sample label counts:\n")
print(table(milo_obj$sample, useNA = "ifany"))

# =======================
# [8] Count Cells Per Sample
# =======================
cat("[9] Counting cells per sample...\n")
milo_obj <- countCells(milo_obj, meta.data = data.frame(colData(milo_obj)), samples = "sample")
cell_counts = nhoodCounts(milo_obj)
cellcount_data = as.data.frame(cell_counts)
write.csv(cellcount_data, file= "results/Milo/cellcount_data.csv", row.names = FALSE)

# =======================
# [9] Design Matrix
# =======================
cat("[10] Creating design matrix...\n")
colData(milo_obj)$condition <- factor(colData(milo_obj)$condition)
colData(milo_obj)$condition <- relevel(colData(milo_obj)$condition, ref = "Control")
design_df <- distinct(data.frame(colData(milo_obj))[, c("sample", "condition")])
rownames(design_df) <- design_df$sample
design_df <- design_df[colnames(nhoodCounts(milo_obj)), , drop = FALSE]

saveRDS(milo_obj,"Milo_Step10_Done_DesignMatrix.rds")
# =======================
# [10] Calculate Distances + Test for DA
# =======================
cat("[11] Calculating distances and testing for differential abundance...\n")
milo_obj <- calcNhoodDistance(milo_obj, d = num_pcs)
saveRDS(milo_obj,"Milo_Step11_Done.rds")

bp <- SnowParam(workers = 6)  # Parallelization
da_results <- testNhoods(milo_obj, design = ~condition, design.df = design_df, BPPARAM = bp)

saveRDS(da_results, "results/Milo/DA_results.rds")
saveRDS(milo_obj, "results/Milo/MiloObject.rds")

# =======================
# [11] Annotate DA Results
# =======================
da_results <- annotateNhoods(milo_obj, da_results, coldata_col = "SingleR_Labels")
head(da_results)

# =======================
# [12] UMAP and DA Plotting
# =======================
cat("[12] Visualizing UMAP and differential abundance results...\n")

umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap.harmony",   # use your Harmony UMAP embedding
  group.by = "SingleR_Labels",  # color by SingleR cell type annotations
  label = TRUE,                 # add text labels
  repel = TRUE,                 # avoid overlapping labels
  label.size = 6                # bigger text
) +
  ggtitle(" ") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


ggsave("results/Milo/umap_plot.png", plot = umap_plot, width = 16, height = 10, dpi = 300)
ggsave("results/Milo/umap_plot.svg", plot = umap_plot, width = 16, height = 10, dpi = 300)

MiloR_DA_beeswarm_plot <- plotDAbeeswarm(da_results, group.by = "SingleR_Labels")
ggsave("results/Milo/MiloR_DA_beeswarm_plot.png", plot = MiloR_DA_beeswarm_plot, width = 16, height = 10, dpi = 100)

# Differential abundance graph
milo_obj <- buildNhoodGraph(milo_obj)
da_plot <- plotNhoodGraphDA(milo_obj, da_results, alpha = 0.05)
ggsave("results/Milo/da_plot.png", plot = da_plot, width = 12, height = 8, dpi = 300)
ggsave("results/Milo/da_plot.svg", plot = da_plot, width = 12, height = 8, dpi = 300)

# Save for future plotting
saveRDS(da_results, "results/Milo/DA_results_ready4plot.rds")
saveRDS(milo_obj, "results/Milo/MiloObject_ready4plot.rds")
saveRDS(umap_plot, "results/Milo/umap_plot.rds")
saveRDS(da_plot, "results/Milo/da_plot.rds")


# Combined DA + UMAP
combined_plot <- umap_plot + da_plot + plot_layout(guides = "collect")
ggsave("results/Milo/combined_plot_annotated.png", plot = combined_plot, width = 16, height = 10, dpi = 300)
ggsave("results/Milo/combined_plot_annotated.svg", plot = combined_plot, width = 16, height = 10, dpi = 300)
# =======================
# [✔] Done
# =======================
end_time <- Sys.time()
cat("[✔] Milo differential abundance testing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
