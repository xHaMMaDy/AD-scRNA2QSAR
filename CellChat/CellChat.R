# ================================================================================
# 💬 CELLCHAT CELL-CELL COMMUNICATION ANALYSIS PIPELINE
# ================================================================================
# 📋 Project: Alzheimer's Disease Cell-Cell Communication Analysis
# 👨‍🔬 Authors: Dr Amr ElHefnawy, Ibrahim Hammad, Mira Moheb, Abdulrahman Wagih, Reem Sharaf, Lorance Gergis
# 📅 Created: 2025-05-01
# 🔄 Last Modified: 2025-05-22
# 🎯 Objective: Comprehensive comparative analysis of cell-cell communication patterns
# 🔬 Method: CellChat analysis comparing Control vs AD conditions
# 📊 Analysis: Ligand-receptor interactions, signaling pathways, and network topology
# 🧠 Focus: Identifying dysregulated communication in Alzheimer's Disease
# ================================================================================
# 📦 LIBRARY LOADING
# ================================================================================
cat("🔄 Loading required libraries...\n")
start_time <- Sys.time()

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
library(dplyr)

end_time <- Sys.time()
cat("✔ Libraries loaded successfully -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# 📁 PATH DEFINITIONS & PARALLEL SETUP
# ================================================================================
input_path <- "C:/Users/HaMMaDy/Desktop/Grad/All_Studies.rds"
output_directory <- "C:/Users/HaMMaDy/Desktop/Grad/results/CellChat"
figure_directory <- file.path(output_directory, "figures")
image_directory <- file.path(figure_directory, "images")

# Create directory structure
directories <- c(
  output_directory,
  figure_directory,
  image_directory,
  file.path(output_directory, "Control"),
  file.path(output_directory, "AD"),
  file.path(output_directory, "Merged")
)

for (dir in directories) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat("📁 Created directory:", dir, "\n")
  }
}

# 🔧 ANALYSIS PARAMETERS
# ================================================================================
analysis_params <- list(
  workers = 6,                    # Parallel processing workers
  min_cells = 10,                 # Minimum cells for communication filtering
  secreted_only = TRUE,           # Focus on secreted signaling only
  thresh_pc = 0.1,                # Threshold for percentage expression
  thresh_fc = 0.05,               # Threshold for fold change
  thresh_p = 0.05,                # Threshold for p-value
  alpha = 0.05                    # Significance level
)

# Configure parallel processing
plan(multisession, workers = analysis_params$workers)
options(future.seed = TRUE)
options(future.globals.maxSize = 168192 * 1024^2)
options(ggrepel.max.overlaps = 100)

cat("📊 Analysis Parameters:\n")
cat("   - Parallel workers:", analysis_params$workers, "\n")
cat("   - Minimum cells:", analysis_params$min_cells, "\n")
cat("   - Secreted signaling only:", analysis_params$secreted_only, "\n")
cat("   - Expression threshold:", analysis_params$thresh_pc, "\n")
cat("   - Fold change threshold:", analysis_params$thresh_fc, "\n")
cat("   - P-value threshold:", analysis_params$thresh_p, "\n\n")

# ================================================================================
# 1️⃣ DATA LOADING AND PREPARATION
# ================================================================================
cat("🎯 Step 1: Loading and Preparing Seurat Object\n")
start_time <- Sys.time()

# 📌 Load Seurat Object
seurat_object <- readRDS(input_path)
seurat_object@meta.data$samples <- as.factor(seurat_object@meta.data$sample)

# 🔹 Split by Condition
cat("🔹 Splitting data by condition...\n")
seurat_control <- subset(seurat_object, subset = condition == "Control")
seurat_alzheimers <- subset(seurat_object, subset = condition == "AD")

# 📊 Data Summary
cat("📊 Data Summary:\n")
cat("   - Total cells:", ncol(seurat_object), "\n")
cat("   - Control cells:", ncol(seurat_control), "\n")
cat("   - AD cells:", ncol(seurat_alzheimers), "\n")
cat("   - Cell types:", length(unique(seurat_object$SingleR_Labels)), "\n")

end_time <- Sys.time()
cat("✔ Step 1 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 2️⃣ CELLCHAT DATABASE SETUP
# ================================================================================
cat("🎯 Step 2: Setting up CellChat Database\n")
start_time <- Sys.time()

# 📌 Load Human Ligand-Receptor Database
cellchat_database <- CellChatDB.human

# 🔹 Filter for Secreted Signaling
if (analysis_params$secreted_only) {
  cellchat_database_filtered <- subsetDB(cellchat_database, search = "Secreted Signaling")
  cat("🔹 Filtered to secreted signaling pathways\n")
} else {
  cellchat_database_filtered <- cellchat_database
}

cat("📊 Database Summary:\n")
cat("   - Total interactions:", nrow(cellchat_database_filtered$interaction), "\n")
cat("   - Signaling pathways:", length(unique(cellchat_database_filtered$interaction$pathway_name)), "\n")

end_time <- Sys.time()
cat("✔ Step 2 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 3️⃣ CELLCHAT OBJECT CONSTRUCTION - CONTROL
# ================================================================================
cat("🎯 Step 3: Building CellChat Object - Control\n")
start_time <- Sys.time()

# 📌 Create CellChat Object
cellchat_control <- createCellChat(
  object = seurat_control,
  group.by = "SingleR_Labels",
  assay = "RNA"
)

# 🔹 Set Database
cellchat_control@DB <- cellchat_database_filtered

# 🔹 Preprocessing Pipeline
cat("🔹 Preprocessing control data...\n")
cellchat_control <- subsetData(cellchat_control)
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)

# 🔹 Communication Analysis
cat("🔹 Computing communication probabilities...\n")
cellchat_control <- computeCommunProb(cellchat_control)
cellchat_control <- filterCommunication(cellchat_control, min.cells = analysis_params$min_cells)

# 🔹 Pathway Analysis
cat("🔹 Analyzing signaling pathways...\n")
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)

# 📊 Control Analysis Summary
cat("📊 Control Analysis Summary:\n")
cat("   - Cell types analyzed:", length(levels(cellchat_control@idents)), "\n")
cat("   - Significant interactions:", nrow(cellchat_control@net$count), "\n")
cat("   - Active pathways:", length(cellchat_control@netP$pathways), "\n")

# 💾 Save Control Object
saveRDS(cellchat_control, file.path(output_directory, "Control", "cellchat_control.rds"))

end_time <- Sys.time()
cat("✔ Step 3 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 4️⃣ CELLCHAT OBJECT CONSTRUCTION - AD
# ================================================================================
cat("🎯 Step 4: Building CellChat Object - AD\n")
start_time <- Sys.time()

# 📌 Create CellChat Object
cellchat_alzheimers <- createCellChat(
  object = seurat_alzheimers,
  group.by = "SingleR_Labels",
  assay = "RNA"
)

# 🔹 Set Database
cellchat_alzheimers@DB <- cellchat_database_filtered

# 🔹 Preprocessing Pipeline
cat("🔹 Preprocessing AD data...\n")
cellchat_alzheimers <- subsetData(cellchat_alzheimers)
cellchat_alzheimers <- identifyOverExpressedGenes(cellchat_alzheimers)
cellchat_alzheimers <- identifyOverExpressedInteractions(cellchat_alzheimers)

# 🔹 Communication Analysis
cat("🔹 Computing communication probabilities...\n")
cellchat_alzheimers <- computeCommunProb(cellchat_alzheimers)
cellchat_alzheimers <- filterCommunication(cellchat_alzheimers, min.cells = analysis_params$min_cells)

# 🔹 Pathway Analysis
cat("🔹 Analyzing signaling pathways...\n")
cellchat_alzheimers <- computeCommunProbPathway(cellchat_alzheimers)
cellchat_alzheimers <- aggregateNet(cellchat_alzheimers)

# 📊 AD Analysis Summary
cat("📊 AD Analysis Summary:\n")
cat("   - Cell types analyzed:", length(levels(cellchat_alzheimers@idents)), "\n")
cat("   - Significant interactions:", nrow(cellchat_alzheimers@net$count), "\n")
cat("   - Active pathways:", length(cellchat_alzheimers@netP$pathways), "\n")

# 💾 Save AD Object
saveRDS(cellchat_alzheimers, file.path(output_directory, "AD", "cellchat_AD.rds"))

end_time <- Sys.time()
cat("✔ Step 4 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 5️⃣ COMPARATIVE ANALYSIS SETUP
# ================================================================================
cat("🎯 Step 5: Setting up Comparative Analysis\n")
start_time <- Sys.time()

# 📌 Create Object List
cellchat_object_list <- list(
  Control = cellchat_control,
  AD = cellchat_alzheimers
)

# 🔹 Merge Objects
cat("🔹 Merging CellChat objects...\n")
cellchat_merged <- mergeCellChat(
  cellchat_object_list,
  add.names = names(cellchat_object_list),
  cell.prefix = TRUE
)

# 💾 Save Merged Object
saveRDS(cellchat_merged, file.path(output_directory, "Merged", "cellchat_merged.rds"))

end_time <- Sys.time()
cat("✔ Step 5 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 6️⃣ BASIC COMPARATIVE VISUALIZATIONS
# ================================================================================
cat("🎯 Step 6: Generating Basic Comparative Visualizations\n")
start_time <- Sys.time()

# 📌 Interaction Count Comparison
cat("🔹 Creating interaction comparison plots...\n")
interaction_count_plot <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2)
)

interaction_weight_plot <- compareInteractions(
  cellchat_merged,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight"
)

# 🔹 Combine and Save
interaction_comparison_combined <- grid.arrange(
  as_grob(interaction_count_plot),
  as_grob(interaction_weight_plot),
  ncol = 2
)

ggsave(
  file.path(image_directory, "interaction_comparison_combined.png"),
  plot = interaction_comparison_combined,
  width = 12, height = 6, dpi = 300
)

saveRDS(interaction_comparison_combined, 
        file.path(figure_directory, "interaction_comparison_combined.rds"))

end_time <- Sys.time()
cat("✔ Step 6 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 7️⃣ NETWORK VISUALIZATIONS
# ================================================================================
cat("🎯 Step 7: Creating Network Visualizations\n")
start_time <- Sys.time()

# 📌 Calculate Group Sizes
group_size_control <- as.numeric(table(cellchat_control@idents))
group_size_alzheimers <- as.numeric(table(cellchat_alzheimers@idents))

# 🔹 Create Circle Plots
cat("🔹 Generating circle network plots...\n")
png(file.path(image_directory, "circle_plots_comparison.png"), 
    width = 18, height = 12, units = "in", res = 300)

par(mfrow = c(2, 2), xpd = TRUE)

netVisual_circle(
  cellchat_control@net$count,
  vertex.weight = group_size_control,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction Count - Control"
)

netVisual_circle(
  cellchat_control@net$weight,
  vertex.weight = group_size_control,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction Strength - Control"
)

netVisual_circle(
  cellchat_alzheimers@net$count,
  vertex.weight = group_size_alzheimers,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction Count - AD"
)

netVisual_circle(
  cellchat_alzheimers@net$weight,
  vertex.weight = group_size_alzheimers,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction Strength - AD"
)

dev.off()

end_time <- Sys.time()
cat("✔ Step 7 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 8️⃣ HEATMAP VISUALIZATIONS
# ================================================================================
cat("🎯 Step 8: Creating Heatmap Visualizations\n")
start_time <- Sys.time()

# 📌 Generate Heatmaps
cat("🔹 Creating comprehensive heatmaps...\n")

# 🔹 Control Heatmaps
heatmap_control_count <- netVisual_heatmap(
  cellchat_control,
  color.heatmap = "Reds",
  title.name = "Interaction Count (Control)"
)

heatmap_control_weight <- netVisual_heatmap(
  cellchat_control,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Interaction Strength (Control)"
)

# 🔹 AD Heatmaps
heatmap_alzheimers_count <- netVisual_heatmap(
  cellchat_alzheimers,
  color.heatmap = "Reds",
  title.name = "Interaction Count (AD)"
)

heatmap_alzheimers_weight <- netVisual_heatmap(
  cellchat_alzheimers,
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Interaction Strength (AD)"
)

# 🔹 Save Combined Heatmaps
png(file.path(image_directory, "heatmap_count_comparison.png"),
    width = 10, height = 12, units = "in", res = 300, bg = "white")
ComplexHeatmap::draw(heatmap_control_count %v% heatmap_alzheimers_count)
dev.off()

png(file.path(image_directory, "heatmap_weight_comparison.png"),
    width = 10, height = 12, units = "in", res = 300, bg = "white")
ComplexHeatmap::draw(heatmap_control_weight %v% heatmap_alzheimers_weight)
dev.off()

end_time <- Sys.time()
cat("✔ Step 8 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 9️⃣ SIGNALING ROLE ANALYSIS
# ================================================================================
cat("🎯 Step 9: Performing Signaling Role Analysis\n")
start_time <- Sys.time()

# 📌 Calculate Network Centrality
cat("🔹 Computing network centrality measures...\n")
future::plan(sequential)

# 🔹 Calculate Weight Range
total_interactions_per_celltype <- sapply(cellchat_object_list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})

weight_range <- c(min(total_interactions_per_celltype), max(total_interactions_per_celltype))

# 🔹 Generate Signaling Role Plots
signaling_plots <- list()
for (i in seq_along(cellchat_object_list)) {
  cellchat_object_list[[i]] <- netAnalysis_computeCentrality(
    cellchat_object_list[[i]], 
    slot.name = 'netP'
  )
  
  signaling_plots[[i]] <- netAnalysis_signalingRole_scatter(
    cellchat_object_list[[i]],
    title = names(cellchat_object_list)[i],
    weight.MinMax = weight_range
  )
}

# 🔹 Combine and Save
signaling_combined <- patchwork::wrap_plots(plots = signaling_plots)
ggsave(
  file.path(image_directory, "signaling_roles_comparison.png"),
  plot = signaling_combined,
  width = 12, height = 6, dpi = 300
)

saveRDS(signaling_combined, 
        file.path(figure_directory, "signaling_roles_comparison.rds"))

end_time <- Sys.time()
cat("✔ Step 9 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 🔟 PATHWAY ANALYSIS
# ================================================================================
cat("🎯 Step 10: Pathway-Level Analysis\n")
start_time <- Sys.time()

# 📌 Get Union of Pathways
pathway_union <- union(
  cellchat_object_list[[1]]@netP$pathways,
  cellchat_object_list[[2]]@netP$pathways
)

cat("🔹 Analyzing", length(pathway_union), "unique pathways...\n")

# 🔹 Generate Pathway Heatmaps
outgoing_heatmap_control <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[1]],
  pattern = "outgoing",
  signaling = pathway_union,
  title = "Control - Outgoing Signals",
  width = 8, height = 10
)

outgoing_heatmap_ad <- netAnalysis_signalingRole_heatmap(
  cellchat_object_list[[2]],
  pattern = "outgoing",
  signaling = pathway_union,
  title = "AD - Outgoing Signals",
  width = 8, height = 10
)

# 🔹 Save Pathway Analysis
png(file.path(image_directory, "pathway_signaling_comparison.png"),
    width = 16, height = 10, units = "in", res = 300, bg = "white")
ComplexHeatmap::draw(outgoing_heatmap_control + outgoing_heatmap_ad)
dev.off()

end_time <- Sys.time()
cat("✔ Step 10 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 1️⃣1️⃣ DIFFERENTIAL EXPRESSION ANALYSIS
# ================================================================================
cat("🎯 Step 11: Differential Expression Analysis\n")
start_time <- Sys.time()

# 📌 Identify Over-expressed Genes
pos_dataset <- "AD"
features_name <- paste0(pos_dataset, ".merged")

cat("🔹 Identifying differentially expressed genes...\n")
cellchat_merged <- identifyOverExpressedGenes(
  cellchat_merged,
  group.dataset = "datasets",
  pos.dataset = pos_dataset,
  features.name = features_name,
  only.pos = FALSE,
  thresh.pc = analysis_params$thresh_pc,
  thresh.fc = analysis_params$thresh_fc,
  thresh.p = analysis_params$thresh_p,
  group.DE.combined = FALSE
)

# 🔹 Map DE Genes to Networks
cat("🔹 Mapping differentially expressed genes to signaling networks...\n")
net <- netMappingDEG(cellchat_merged, features.name = features_name, variable.all = TRUE)

# 🔹 Extract Up/Down-regulated Communications
net_up <- subsetCommunication(
  cellchat_merged,
  net = net,
  datasets = "AD",
  ligand.logFC = analysis_params$thresh_fc,
  receptor.logFC = NULL
)

net_down <- subsetCommunication(
  cellchat_merged,
  net = net,
  datasets = "Control",
  ligand.logFC = -analysis_params$thresh_fc,
  receptor.logFC = NULL
)

# 📊 DE Analysis Summary
cat("📊 Differential Expression Summary:\n")
cat("   - Up-regulated interactions:", nrow(net_up), "\n")
cat("   - Down-regulated interactions:", nrow(net_down), "\n")

# 💾 Save DE Results
write.csv(net_up, file.path(output_directory, "upregulated_interactions.csv"))
write.csv(net_down, file.path(output_directory, "downregulated_interactions.csv"))

end_time <- Sys.time()
cat("✔ Step 11 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 1️⃣2️⃣ BUBBLE PLOT VISUALIZATIONS
# ================================================================================
cat("🎯 Step 12: Creating Bubble Plot Visualizations\n")
start_time <- Sys.time()

# 📌 Key Cell Type Interactions
cat("🔹 Generating bubble plots for key interactions...\n")

# 🔹 Astrocyte-Microglia Interactions
astrocyte_microglia_plot <- netVisual_bubble(
  cellchat_merged,
  sources.use = "Astrocyte",
  targets.use = "Microglia-PVM",
  comparison = c(1, 2),
  angle.x = 45
)

# 🔹 Microglia to Multiple Targets
microglia_multi_plot <- netVisual_bubble(
  cellchat_merged,
  sources.use = "Microglia-PVM",
  targets.use = c("Astrocyte", "Endothelial", "L6 IT", "L2/3 IT", "L4 IT"),
  comparison = c(1, 2),
  angle.x = 45
)

# 🔹 Save Bubble Plots
ggsave(
  file.path(image_directory, "astrocyte_microglia_interactions.png"),
  plot = astrocyte_microglia_plot,
  width = 10, height = 8, dpi = 300
)

ggsave(
  file.path(image_directory, "microglia_multi_interactions.png"),
  plot = microglia_multi_plot,
  width = 14, height = 10, dpi = 300
)

end_time <- Sys.time()
cat("✔ Step 12 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 1️⃣3️⃣ NETWORK SIMILARITY ANALYSIS
# ================================================================================
cat("🎯 Step 13: Network Similarity Analysis\n")
start_time <- Sys.time()

# 📌 Compute Network Similarity
cat("🔹 Computing functional network similarity...\n")
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional")

# 🔹 Visualize Network Embedding
functional_embedding_plot <- netVisual_embeddingPairwise(
  cellchat_merged, 
  type = "functional", 
  label.size = 3.5
)

# 🔹 Structural Network Analysis
cat("🔹 Computing structural network similarity...\n")
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "structural")
cellchat_merged <- netEmbedding(cellchat_merged, type = "structural")
cellchat_merged <- netClustering(cellchat_merged, type = "structural")

# 🔹 Visualize Structural Embedding
structural_embedding_plot <- netVisual_embeddingPairwise(
  cellchat_merged, 
  type = "structural", 
  label.size = 3.5
)

# 💾 Save Final Merged Object
saveRDS(cellchat_merged, file.path(output_directory, "cellchat_merged_final.rds"))

end_time <- Sys.time()
cat("✔ Step 13 completed -", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n\n")

# ================================================================================
# 📊 SESSION DOCUMENTATION
# ================================================================================
cat("📊 Documenting session information...\n")
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(output_directory, "session_info.txt"))

# 🎉 PIPELINE COMPLETION
# ================================================================================
total_end_time <- Sys.time()
cat("🎉 CELLCHAT ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("📁 All results saved to:", output_directory, "\n")
cat("📊 Key outputs:\n")
cat("   - Control CellChat object: Control/cellchat_control.rds\n")
cat("   - AD CellChat object: AD/cellchat_AD.rds\n")
cat("   - Merged analysis: cellchat_merged_final.rds\n")
cat("   - Visualizations: figures/images/\n")
cat("   - DE interactions: upregulated_interactions.csv, downregulated_interactions.csv\n")
cat("⏱️  Total analysis time:", round(difftime(total_end_time, start_time, units = "mins"), 2), "minutes\n")
cat("================================================================================\n")