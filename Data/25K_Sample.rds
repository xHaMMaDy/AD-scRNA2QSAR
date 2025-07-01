#=========================================================================================
# Script: CellChat Analysis & Comparison (Control vs AD) 
# Purpose: Comprehensive cell-cell communication analysis comparing control and AD conditions
# Authors: Dr Amr ElHefnawy - Ibrahim Hammad - Mira Moheb - Abdulrahman Wagih - Reem Sharaf - Lorance Gergis
# Created: 22-05-2025
#=========================================================================================

#=========================================================================================
# 📌 Load Required Libraries
# Load all necessary packages for single-cell analysis and visualization
#=========================================================================================
library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(ggplotify)
library(future)
library(future.apply)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(grid)
#=========================================================================================
# 📁 Define All Paths (Project structure) + Parallel Setup
#=========================================================================================
seurat_rds_path <- "C:/Users/HaMMaDy/Desktop/Grad/All_Studies.rds"
output_directory <- "C:/Users/HaMMaDy/Desktop/Grad/results/CellChat"
setwd(output_directory)
figure_directory <- "figure_objects_new"
image_directory <- file.path(figure_directory, "images")
output_control_directory <- file.path(output_directory, "Control")
output_alzheimers_directory <- file.path(output_directory, "AD")
output_merged_directory <- file.path(output_directory, "Merged")

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(image_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(output_control_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(output_alzheimers_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(output_merged_directory, recursive = TRUE, showWarnings = FALSE)

plan(multisession, workers = 6)
options(future.seed = TRUE)
options(future.globals.maxSize = 168192 * 1024^2)
options(ggrepel.max.overlaps = 100)
#=========================================================================================
# 🔹 1. Load Data and Split Control vs AD
# Import Seurat object and separate samples by experimental condition
#=========================================================================================
cat("[ 1 ] Loading Seurat object...\n")
start_time <- Sys.time()
seurat_object <- readRDS(seurat_rds_path)
seurat_object@meta.data$samples <- as.factor(seurat_object@meta.data$sample)
cat("[ 1.1 ] Splitting into Control vs AD based on condition metadata...\n")
# Extract control & Alzheimer's disease based on condition 
seurat_control <- subset(seurat_object, subset = condition == "Control")
seurat_alzheimers <- subset(seurat_object, subset = condition == "AD")

# Store cell type labels for control & AD samples for later reference
cell_labels_control <- Idents(seurat_control)
cell_labels_alzheimers <- Idents(seurat_alzheimers)

end_time <- Sys.time()
cat("[ ✔ ] Data loaded and split in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")
#=========================================================================================
# 🔹 2. Configure CellChatDB
#=========================================================================================
cat("[ 2 ] Setting up CellChatDB for human ligand-receptor interactions...\n")
# Load human ligand-receptor interaction database for analysis
cellchat_database <- CellChatDB.human
#=========================================================================================
# 🔹 3. Build CellChat Object for Control & AD Condition
#=========================================================================================
cat("[ 3 ] Running CellChat analysis for CONTROL & AD conditions...\n")

# Initialize CellChat object for control & AD samples using SingleR cell type annotations
cellchat_control <- createCellChat(
  object = seurat_control, 
  group.by = "SingleR_Labels", 
  assay = "RNA"
)
cellchat_alzheimers <- createCellChat(
  object = seurat_alzheimers, 
  group.by = "SingleR_Labels", 
  assay = "RNA"
)

start_time <- Sys.time()
# Filter database to focus on secreted signaling pathways only
cellchat_database_filtered <- subsetDB(cellchat_database, search = "Secreted Signaling")
cellchat_control@DB <- cellchat_database_filtered

# Subset expression data to include only genes present in the interaction database
cellchat_control <- subsetData(cellchat_control)
cellchat_alzheimers <- subsetData(cellchat_alzheimers)

# Identify significantly over-expressed genes in each cell type
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_alzheimers <- identifyOverExpressedGenes(cellchat_alzheimers)

# Identify over-expressed ligand-receptor interactions based on gene expression
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_alzheimers <- identifyOverExpressedInteractions(cellchat_alzheimers)

# Extract significant ligand-receptor pairs for downstream analysis
ligand_receptor_interactions_control <- cellchat_control@LR$LRsig
ligand_receptor_interactions_alzheimers <- cellchat_alzheimers@LR$LRsig

# Calculate communication probabilities between all cell type pairs
cellchat_control <- computeCommunProb(cellchat_control)
cellchat_alzheimers <- computeCommunProb(cellchat_alzheimers)

# Filter communications to include only those with sufficient cell numbers (min 10 cells)
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10)
cellchat_alzheimers <- filterCommunication(cellchat_alzheimers, min.cells = 10)

# Extract detailed communication network data for analysis
communication_network_control <- subsetCommunication(cellchat_control)
communication_network_alzheimers <- subsetCommunication(cellchat_alzheimers)

# Extract pathway-level communication probabilities for pathway analysis
pathway_network_control <- subsetCommunication(cellchat_control, slot.name = "netP")
pathway_network_alzheimers <- subsetCommunication(cellchat_alzheimers, slot.name = "netP")

# Calculate communication probabilities at the pathway level for pathway analysis
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_alzheimers <- computeCommunProbPathway(cellchat_alzheimers)

# Aggregate cell-cell communication network for visualization
cellchat_control <- aggregateNet(cellchat_control)
cellchat_alzheimers <- aggregateNet(cellchat_alzheimers)

end_time <- Sys.time()
cat("[ ✔ ] Control processing completed in",round(difftime(end_time, start_time, units = "secs"), 2), "seconds.\n")

# Save processed CellChat object for control condition for future use
saveRDS(cellchat_control,file.path(output_control_directory, "cellchat_control_processed.rds"))

# loading option for control data
cellchat_control <- readRDS("C:/Users/Hammady/Desktop/Grad/results/CellChat/Control/cellchat_control.rds")
cellchat_alzheimers <- readRDS("C:/Users/Hammady/Desktop/Grad/results/CellChat/AD/cellchat_AD.rds")
#=========================================================================================
# 🔹 4. Merge CellChat Objects for Comparative Analysis
#=========================================================================================
cat("[ 5 ] Merging CellChat objects for comparative analysis...\n")

# Create named list containing both CellChat objects for comparison
cellchat_object_list <- list(
  Control = cellchat_control, 
  AD = cellchat_alzheimers
)

# Merge CellChat objects for comparative analysis with proper naming
cellchat_merged <- mergeCellChat(
  cellchat_object_list, 
  add.names = names(cellchat_object_list), 
  cell.prefix = TRUE
)

# Save merged CellChat object for future comparative analyses
saveRDS(
  cellchat_merged, 
  file.path(output_directory, "cellchat_merged_analysis.rds")
)

# Alternative loading option for pre-processed merged data (commented out for safety)
cellchat_merged <- readRDS("C:/Users/Hammady/Desktop/Grad/results/CellChat/Merged/cellchat_merged.rds")
#=========================================================================================
# 🔹 6. Comprehensive Visualization Suite
# Generate complete set of comparative visualizations
#=========================================================================================
cat("[ 6 ] Generating visualizations\n")
#-----------------------------------------------------------------------------------------
# 6.1 Interaction Count and Strength Comparisons
# Compare total interactions and their strengths between conditions
#-----------------------------------------------------------------------------------------
cat("[ 6.1 ] Creating interaction comparison plots...\n")
# Generate comparison plot for total number of interactions between conditions
interaction_count_comparison <- compareInteractions(
  cellchat_merged, 
  show.legend = FALSE, 
  group = c(1, 2)
)

# Generate comparison plot for interaction strength between conditions
interaction_weight_comparison <- compareInteractions(
  cellchat_merged, 
  show.legend = FALSE, 
  group = c(1, 2), 
  measure = "weight"
)

# Combined plot layout using grid arrangement for side-by-side comparison
interaction_comparison_combined <- grid.arrange(
  as_grob(interaction_count_comparison), 
  as_grob(interaction_weight_comparison), 
  ncol = 2
)

ggsave(
  file.path(image_directory, "interaction_comparison_combined.png"), 
  plot = interaction_comparison_combined, 
  width = 12, 
  height = 6, 
  dpi = 300
)
saveRDS(
  interaction_comparison_combined,
  file = file.path(figure_directory, "interaction_comparison_combined.rds")
)
#-----------------------------------------------------------------------------------------
# 6.2 Circle Network Visualizations
#-----------------------------------------------------------------------------------------
cat("[ 6.2 ] Creating circular network visualizations...\n")

# Calculate group sizes for control & AD  condition (number of cells per cell type)
group_size_control <- as.numeric(table(cellchat_control@idents))
group_size_alzheimers <- as.numeric(table(cellchat_alzheimers@idents))

# Set up 2x2 plot layout for circle visualizations
par(mfrow = c(2, 2), xpd = TRUE) # Allow titles to extend outside plot area

# Circle plot showing interaction counts for control condition
circle_plot_control_count <- netVisual_circle(
  cellchat_control@net$count, 
  vertex.weight = group_size_control, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  title.name = "Interaction Count - Control"
)

# Circle plot showing interaction strength for control condition
circle_plot_control_weight <- netVisual_circle(
  cellchat_control@net$weight, 
  vertex.weight = group_size_control, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  title.name = "Interaction Strength - Control"
)

# Circle plot showing interaction counts for AD condition
circle_plot_alzheimers_count <- netVisual_circle(
  cellchat_alzheimers@net$count, 
  vertex.weight = group_size_alzheimers, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  title.name = "Interaction Count - AD"
)

# Generate circle plot showing interaction strength for AD condition
circle_plot_alzheimers_weight <- netVisual_circle(
  cellchat_alzheimers@net$weight, 
  vertex.weight = group_size_alzheimers, 
  weight.scale = TRUE, 
  label.edge = FALSE, 
  title.name = "Interaction Strength - AD"
)

# Convert final circle plot to grob for saving
circle_plots_combined <- as_grob(circle_plot_alzheimers_weight)

ggsave(
  file.path(image_directory, "circle_plots_comparison.png"), 
  plot = circle_plots_combined, 
  width = 18, 
  height = 12, 
  dpi = 300
)
saveRDS(
  circle_plots_combined,
  file = file.path(figure_directory, "circle_plots_comparison.rds")
)
#-----------------------------------------------------------------------------------------
# 6.3 Heatmap Visualizations
#-----------------------------------------------------------------------------------------
cat("[ 6.3 ] Creating comprehensive heatmap visualizations...\n")

# Generate heatmap showing interaction counts for control condition
heatmap_control_count <- netVisual_heatmap(
  cellchat_control, 
  color.heatmap = c('Reds'), 
  title.name = "Interaction Count (Control)"
)
# Generate heatmap showing interaction strength for control condition
heatmap_control_weight <- netVisual_heatmap(
  cellchat_control, 
  measure = "weight", 
  color.heatmap = c('Reds'), 
  title.name = "Interaction Strength (Control)"
)

# Combine control condition heatmaps for side-by-side comparison
heatmap_control_combined <- heatmap_control_count + heatmap_control_weight

png(
  file.path(image_directory, "heatmap_control_combined.png"), 
  width = 12, 
  height = 8, 
  units = "in", 
  res = 300, 
  bg = "white"
)
ComplexHeatmap::draw(heatmap_control_combined)
dev.off()

# Generate heatmap showing interaction counts for AD condition
heatmap_alzheimers_count <- netVisual_heatmap(
  cellchat_alzheimers, 
  color.heatmap = c("Reds"), 
  title.name = "Interaction Count (AD)"
)
# Generate heatmap showing interaction strength for AD condition
heatmap_alzheimers_weight <- netVisual_heatmap(
  cellchat_alzheimers, 
  measure = "weight", 
  color.heatmap = c("Reds"), 
  title.name = "Interaction Strength (AD)"
)
png(
  file.path(image_directory, "heatmap_alzheimers_combined.png"), 
  width = 12, 
  height = 8, 
  units = "in", 
  res = 300, 
  bg = "white"
)
ComplexHeatmap::draw(heatmap_alzheimers_count + heatmap_alzheimers_weight)
dev.off()

png(
  file.path(image_directory, "heatmap_count_condition_comparison.png"), 
  width = 10, 
  height = 12, 
  units = "in", 
  res = 300, 
  bg = "white"
)
ComplexHeatmap::draw(heatmap_control_count %v% heatmap_alzheimers_count)
dev.off()

png(
  file.path(image_directory, "heatmap_weight_condition_comparison.png"), 
  width = 10, 
  height = 12, 
  units = "in", 
  res = 300, 
  bg = "white"
)
ComplexHeatmap::draw(heatmap_control_weight %v% heatmap_alzheimers_weight)
dev.off()

#-----------------------------------------------------------------------------------------
# 6.4 Cell Type Aggregation Analysis
#-----------------------------------------------------------------------------------------
cat("[ 6.4 ] Performing cell type aggregation analysis...\n")

grouped_cell_types <- c(
  rep("Oligodendrocyte", 2), 
  rep("Astrocyte", 2), 
  rep("Microglia-PVM", 2), 
  rep("L6 CT", 2), 
  rep("L6 IT", 2)
)

grouped_cell_types <- factor(
  grouped_cell_types, 
  levels = c("Oligodendrocyte", "Astrocyte", "Microglia-PVM", "L6 CT", "L6 IT")
)

# Merge interactions based on defined cell type groups for simplified analysis
cellchat_object_list <- lapply(cellchat_object_list, function(cellchat_obj) {
  mergeInteractions(cellchat_obj, grouped_cell_types)
})
# Merge updated objects for comparative analysis
cellchat_grouped <- mergeCellChat(
  cellchat_object_list, 
  add.names = names(cellchat_object_list)
)
# Calculate maximum weights across datasets 
weight_maximum <- getMaxWeight(
  cellchat_object_list, 
  slot.name = c("idents", "net", "net"), 
  attribute = c("idents", "count", "count.merged")
)
#-----------------------------------------------------------------------------------------
# 6.5 Grouped Circle Visualizations
#-----------------------------------------------------------------------------------------
cat("[ 6.5 ] Creating circle visualizations for aggregated cell groups...\n")

png(
  file.path(image_directory, "grouped_circle_interactions_comparison.png"), 
  width = 12, 
  height = 6, 
  units = "in", 
  res = 300
)
par(mfrow = c(1, 2), xpd = TRUE)
for (condition_index in 1:length(cellchat_object_list)) {
  netVisual_circle(
    cellchat_object_list[[condition_index]]@net$count.merged, 
    weight.scale = TRUE,
    label.edge = TRUE, 
    edge.weight.max = weight_maximum[3],
    edge.width.max = 12,
    title.name = paste0(
      "Number of interactions - ", 
      names(cellchat_object_list)[condition_index]
    )
  )
}
dev.off()

#-----------------------------------------------------------------------------------------
# 6.6 Signaling Role Analysis
#-----------------------------------------------------------------------------------------
cat("[ 6.6 ] Performing comprehensive signaling role analysis...\n")
# Extract detailed communication networks from merged object for analysis
detailed_communication_network <- subsetCommunication(cellchat_merged)
pathway_communication_network <- subsetCommunication(
  cellchat_merged, 
  slot.name = "netP"
)
# Calculate total number of interactions per cell type across datasets
total_interactions_per_celltype <- sapply(cellchat_object_list, function(cellchat_obj) {
  rowSums(cellchat_obj@net$count) + 
    colSums(cellchat_obj@net$count) - 
    diag(cellchat_obj@net$count)
})
# Define weight range for consistent dot sizing across plots
weight_min_max_range <- c(
  min(total_interactions_per_celltype), 
  max(total_interactions_per_celltype)
)

future::plan(sequential)
# Compute network centrality measures and generate signaling role scatter plots
for (condition_index in seq_along(cellchat_object_list)) {
  # Calculate centrality measures for each cell type in current condition
  cellchat_object_list[[condition_index]] <- netAnalysis_computeCentrality(
    cellchat_object_list[[condition_index]], 
    slot.name = 'netP'
  )
  # Generate scatter plot showing signaling roles (outgoing vs incoming)
  signaling_role_scatter_plot <- netAnalysis_signalingRole_scatter(
    cellchat_object_list[[condition_index]],
    title = names(cellchat_object_list)[condition_index],
    weight.MinMax = weight_min_max_range
  )
  ggsave(
    file.path(
      image_directory, 
      paste0(
        "signaling_role_scatter_", 
        names(cellchat_object_list)[condition_index], 
        ".png"
      )
    ),
    plot = signaling_role_scatter_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

#-----------------------------------------------------------------------------------------
# 6.7 Combined Signaling Role Scatter Plots
#-----------------------------------------------------------------------------------------
cat("[ 6.7 ] Generating combined signaling role scatter plots...\n")

signaling_scatter_plots <- list()

# Generate signaling role scatter plots for both conditions
for (condition_index in 1:length(cellchat_object_list)) {
  signaling_scatter_plots[[condition_index]] <- netAnalysis_signalingRole_scatter(
    cellchat_object_list[[condition_index]],
    title = names(cellchat_object_list)[condition_index],
    weight.MinMax = weight_min_max_range
  )
}

# Combine all scatter plots
signaling_roles_combined_plot <- patchwork::wrap_plots(plots = signaling_scatter_plots)
ggsave(
  file.path(image_directory, "signaling_roles_combined_comparison.png"), 
  plot = signaling_roles_combined_plot, 
  width = 12, 
  height = 6, 
  dpi = 300
)
saveRDS(
  signaling_roles_combined_plot,
  file = file.path(figure_directory, "signaling_roles_combined_comparison.rds")
)

#-----------------------------------------------------------------------------------------
# 6.8 Cell Type-Specific Signaling Changes Analysis
# Analyze how specific cell types change their signaling patterns
#-----------------------------------------------------------------------------------------
cat("[ 6.8 ] Analyzing cell type-specific signaling changes...\n")

# Generate signaling change scatter plots for specific cell types of interest
endothelial_signaling_changes <- netAnalysis_signalingChanges_scatter(
  cellchat_object_list, 
  idents.use = "Endothelial"
)
layer4_it_signaling_changes <- netAnalysis_signalingChanges_scatter(
  cellchat_object_list, 
  idents.use = "L4 IT"
)
microglia_signaling_changes <- netAnalysis_signalingChanges_scatter(
  cellchat_object_list, 
  idents.use = "Microglia-PVM"
)
oligodendrocyte_precursor_signaling_changes <- netAnalysis_signalingChanges_scatter(
  cellchat_object_list, 
  idents.use = "OPC"
)

# Combine all cell type-specific plots into 2x2 grid layout
celltype_signaling_changes_combined <- patchwork::wrap_plots(
  plots = list(
    endothelial_signaling_changes, 
    layer4_it_signaling_changes, 
    microglia_signaling_changes, 
    oligodendrocyte_precursor_signaling_changes
  ), 
  ncol = 2
)

ggsave(
  file.path(image_directory, "celltype_signaling_changes_comprehensive.png"), 
  plot = celltype_signaling_changes_combined,
  width = 18, 
  height = 12, 
  dpi = 300
)
#-----------------------------------------------------------------------------------------
# 6.9 Pathway-Level Signaling Pattern Analysis
# Compare signaling patterns at the pathway level across conditions
#-----------------------------------------------------------------------------------------
cat("[ 6.9 ] Analyzing pathway-level signaling patterns across conditions...\n")

# Identify union of all pathways present in both conditions
pathway_union <- union(
  cellchat_object_list[[1]]@netP$pathways, 
  cellchat_object_list[[2]]@netP$pathways
)

# Generate heatmaps for outgoing signaling patterns
outgoing_heatmap_control <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[1]], 
  pattern = "outgoing", 
  signaling = pathway_union,
  title = "Control - Outgoing Signals",
  width = 8, 
  height = 10
)
outgoing_heatmap_alzheimers <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[2]], 
  pattern = "outgoing", 
  signaling = pathway_union,
  title = names(cellchat_object_list)[2],
  width = 8, 
  height = 10
)

# Generate heatmaps for incoming signaling patterns
incoming_heatmap_control <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[1]], 
  pattern = "incoming", 
  signaling = pathway_union,
  title = "Control - Incoming Signals",
  width = 8, 
  height = 10,
  color.heatmap = "GnBu"
)
incoming_heatmap_alzheimers <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[2]], 
  pattern = "incoming", 
  signaling = pathway_union,
  title = names(cellchat_object_list)[2],
  width = 8, 
  height = 10,
  color.heatmap = "GnBu"
)

# Assign unique names to heatmap objects for legend handling
outgoing_heatmap_control@name <- "OutgoingSignaling_Control"
outgoing_heatmap_alzheimers@name <- "OutgoingSignaling_AD"
incoming_heatmap_control@name <- "IncomingSignaling_Control"
incoming_heatmap_alzheimers@name <- "IncomingSignaling_AD"

# Combine outgoing and incoming heatmaps for each condition
control_signaling_patterns <- outgoing_heatmap_control + incoming_heatmap_control
alzheimers_signaling_patterns <- outgoing_heatmap_alzheimers + incoming_heatmap_alzheimers

control_signaling_grob <- grid.grabExpr(
  draw(control_signaling_patterns, ht_gap = unit(4, "mm"), merge_legends = TRUE), 
  newpage = FALSE
)
alzheimers_signaling_grob <- grid.grabExpr(
  draw(alzheimers_signaling_patterns, ht_gap = unit(4, "mm"), merge_legends = TRUE), 
  newpage = FALSE
)

png(
  file.path(image_directory, "comprehensive_signaling_patterns_comparison.png"), 
  width = 14, 
  height = 12, 
  units = "in", 
  res = 300, 
  bg = "white"
)
grid.arrange(control_signaling_grob, alzheimers_signaling_grob, nrow = 2)
dev.off()

# Generate heatmaps showing all signaling patterns (combined outgoing and incoming)
all_signaling_heatmap_control <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[1]],
  pattern = "all",
  signaling = pathway_union,
  title = names(cellchat_object_list)[1],
  width = 6,
  height = 8,
  color.heatmap = "OrRd"
)
all_signaling_heatmap_alzheimers <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[2]],
  pattern = "all",
  signaling = pathway_union,
  title = names(cellchat_object_list)[2],
  width = 6,
  height = 8,
  color.heatmap = "OrRd"
)
png(
  file.path(image_directory, "all_signaling_patterns_comparison.png"), 
  width = 8, 
  height = 6, 
  units = "in", 
  res = 300, 
  bg = "white"
)
ComplexHeatmap::draw(all_signaling_heatmap_control + all_signaling_heatmap_alzheimers)
dev.off()

# ============================ 
# 25-02-2025
# ============================ 

# Compute network similarity for functional analysis
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
#> Compute signaling network similarity for datasets 1 2

# Perform manifold learning and network clustering for functional networks
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat_merged <- netClustering(cellchat_merged, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2

# Visualization in 2D-space for functional networks
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

# Compute network similarity for structural analysis
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural")

# Visualization in 2D-space for structural networks
netVisual_embeddingPairwise(cellchat_merged, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat_merged, type = "structural", nCol = 2)

### Part III: Identify the up-regulated and down-regulated signaling ligand-receptor pairs ###

# Specific cell type interaction analyses
netVisual_bubble(cellchat_merged, 
                 sources.use = "Astrocyte", 
                 targets.use = c("Microglia-PVM"), 
                 comparison = c(1, 2), 
                 angle.x = 45)

netVisual_bubble(cellchat_merged, 
                 sources.use = "Microglia-PVM", 
                 targets.use = c("Astrocyte","Endothelial","L6 IT","L2/3 IT","L4 IT","L5 IT","L6 cT"), 
                 comparison = c(1, 2), 
                 angle.x = 45)

netVisual_bubble(cellchat_merged, 
                 sources.use = "VLMC", 
                 targets.use = c("OPC"), 
                 comparison = c(1, 2), 
                 angle.x = 45, 
                 line.on = FALSE)

netVisual_bubble(cellchat_merged, 
                 sources.use = "Sst Chodl", 
                 targets.use = c("Sst Chodl"), 
                 comparison = c(1, 2), 
                 angle.x = 45, 
                 line.on = FALSE)

netVisual_bubble(cellchat_merged, 
                 sources.use = "L6b", 
                 targets.use = c("Endothelial"), 
                 comparison = c(1, 2), 
                 angle.x = 45, 
                 line.on = FALSE)

netVisual_bubble(cellchat_merged, 
                 sources.use = "L6b", 
                 targets.use = c("L2/3 IT"), 
                 comparison = c(1, 2), 
                 angle.x = 45, 
                 line.on = FALSE)

# Generate comprehensive bubble plot for all pathways
allpathwayz <- netVisual_bubble(cellchat_merged, 
                                sources.use = c(1:24), 
                                targets.use = c(1:24), 
                                comparison = c(1, 2), 
                                angle.x = 45)

# Generate individual bubble plots for all source-target combinations
for (s in 1:24) {
  for (t in 1:24) {
    
    p <- netVisual_bubble(
      cellchat_merged,
      sources.use = s,
      targets.use = t,
      comparison  = c(1, 2),
      angle.x     = 45,
      line.on     = FALSE        # avoids the seq() error
    )
    
    ggsave(
      sprintf("%s/bubble_plot_source_%02d_target_%02d.png",
              Visualbubble, s, t),
      plot   = p,
      width  = 8,
      height = 6,
      dpi    = 300
    )
  }
}

# Comparative analysis for specific neuronal interactions
gg1 <- netVisual_bubble(cellchat_merged, 
                        sources.use = "Astrocyte", 
                        targets.use = c(6,7,8), 
                        comparison = c(1, 2), 
                        max.dataset = 2, 
                        title.name = "Neurons signaling pathways", 
                        angle.x = 45, 
                        remove.isolate = T)

gg2 <- netVisual_bubble(cellchat_merged, 
                        sources.use = 1, 
                        targets.use = c(6,7,8), 
                        comparison = c(1, 2), 
                        max.dataset = 1, 
                        title.name = "macrophage signaling pathways", 
                        angle.x = 45, 
                        remove.isolate = T)

gg1 + gg2

## Identify the upgulated and down-regulated signaling ligand-receptor pairs based on differential expression analysis (DEA)
pos.dataset = "AD"
features.name = paste0(pos.dataset, ".merged")

# Identify over-expressed genes
cellchat_merged <- identifyOverExpressedGenes(cellchat_merged, 
                                              group.dataset = "datasets", 
                                              pos.dataset = pos.dataset, 
                                              features.name = features.name, 
                                              only.pos = FALSE, 
                                              thresh.pc = 0.1, 
                                              thresh.fc = 0.05,
                                              thresh.p = 0.05, 
                                              group.DE.combined = FALSE) 

# Map differentially expressed genes to signaling networks
net <- netMappingDEG(cellchat_merged, features.name = features.name, variable.all = TRUE)

# Extract up-regulated and down-regulated communications
net.up <- subsetCommunication(cellchat_merged, 
                              net = net, 
                              datasets = "AD",
                              ligand.logFC = 0.05, 
                              receptor.logFC = NULL)

net.down <- subsetCommunication(cellchat_merged, 
                                net = net, 
                                datasets = "Control",
                                ligand.logFC = -0.05, 
                                receptor.logFC = NULL)

# Extract gene subsets from interaction pairs
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merged)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merged)

# 
# pairLR.use.up = net.up[, "interaction_name", drop = F]
# gg1 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.up, sources.use = 5, targets.use = c(6,7,8), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = "Neurons signaling pathways")
# #> Comparing communications on a merged object
# pairLR.use.down = net.down[, "interaction_name", drop = F]
# gg2 <- netVisual_bubble(cellchat_merged, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(6,7,8), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = "macrophage signaling pathways")
# #> Comparing communications on a merged object
# gg1 + gg2
#

# Visualize the enriched ligands in the first condition using wordcloud
computeEnrichmentScore(net.down, species = 'human', variable.both = TRUE)

# Visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'human', variable.both = TRUE)

# 
# ## Compare the signaling gene expression distribution between different datasets
# cellchat_merged@meta$datasets = factor(cellchat_merged@meta$datasets, levels = c("Control", "AD")) # set factor level
# plotGeneExpression(cellchat_merged, signaling = "PTN", split.by = "datasets", colors.ggplot = T)
#
#=========================================================================================
# 📊 Session Documentation and Analysis Completion
#=========================================================================================
session <- capture.output(sessionInfo())
writeLines(session, file.path(output_directory, "session_info.txt"))
cat("[ ✔ ] CellChat analysis complete! All visualizations saved to:", figure_directory, "directory.\n")