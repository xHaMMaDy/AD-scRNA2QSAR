library(dplyr)
library(Seurat)
library(SeuratObject)
library(readxl)
library(patchwork)
library(tibble)
library(utils)
library(readr)
library(hdf5r)
library(base)
library(devtools)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridExtra)

##################### Study 1 (GSE138852) ##################### 


# Study GSE138852 Import , Region : Entorhinal Cortex
GSE138852_counts <- read.csv("C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE138852/GSE138852_counts.csv", row.names = 1)
GSE138852_covariates <- read.csv("C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE138852/GSE138852_covariates.csv")

# Create a Seurat object
seurat_obj <- CreateSeuratObject(counts = GSE138852_counts, project = "Entorhinal Cortex")
metadata = seurat_obj@meta.data

# condition & Tissue Condition
cell_names <- colnames(seurat_obj)
condition <- rep(NA, length(cell_names))
ad_pattern <- grepl("_AD[0-9]+(_|$)", cell_names)
condition[ad_pattern] <- "AD"
ct_pattern <- grepl("_Ct[0-9]+(_|$)", cell_names)
condition[ct_pattern] <- "Control"
seurat_obj$condition <- condition
 seurat_obj$Tissue_Condition <- paste0(
  seurat_obj$orig.ident, 
  "_", 
  seurat_obj$condition
)
head(seurat_obj@meta.data)

# SAVE OBJ
saveRDS(seurat_obj, "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE138852_obj_seurat.rds")


##################### Study 2 (GSE157827) ##################### 

data_dirs <- list(
  "AD1" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD1",
  "AD2" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD2",
  "AD4" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD4",
  "AD5" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD5",
  "AD6" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD6",
  "AD8" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD8",
  "AD9" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD9",
  "AD10" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD10",
  "AD13" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD13",
  "AD19" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD19",
  "AD20" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD20",
  "AD21" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/AD21",
  "NC3" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC3",
  "NC7" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC7",
  "NC11" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC11",
  "NC12" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC12",
  "NC14" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC14",
  "NC15" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC15",
  "NC16" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC16",
  "NC17" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC17",
  "NC18" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE157827/NC18"
)

seurat_list <- lapply(names(data_dirs), function(sample_name) {
  data <- Read10X(data_dirs[[sample_name]])
  seurat_obj <- CreateSeuratObject(counts = data, project = "Prefrontal cortex")
  seurat_obj$sample <- sample_name  # Add metadata to distinguish samples
  return(seurat_obj)
})
names(seurat_list) <- names(data_dirs)  # Assign sample names

combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list)
)

# condition & Tissue Condition
cell_names <- colnames(combined_seurat)
condition <- rep(NA, length(cell_names))
ad_pattern <- grepl("^AD[0-9]+_", cell_names)
condition[ad_pattern] <- "AD"
nc_pattern <- grepl("^NC[0-9]+_", cell_names)
condition[nc_pattern] <- "Control"
combined_seurat$condition <- condition
combined_seurat$Tissue_Condition <- paste0(
  combined_seurat$orig.ident, 
  "_", 
  combined_seurat$condition
)
head(combined_seurat@meta.data)

# SAVE OBJ -- GSE157827 --
saveRDS(combined_seurat, "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE157827_obj_seurat.rds")


##################### Study 3 (GSE175814) ##################### 

data_dirs <- list(
  "A1_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE175814/A1_AD",
  "A2_Control" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE175814/A2_Control",
  "A3_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE175814/A3_AD",
  "A4_Control" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE175814/A4_Control"

)

seurat_list <- lapply(names(data_dirs), function(sample_name) {
  data <- Read10X(data_dirs[[sample_name]])
  seurat_obj <- CreateSeuratObject(counts = data, project = "Hippocampus")
  seurat_obj$sample <- sample_name  # Add metadata to distinguish samples
  return(seurat_obj)
})
names(seurat_list) <- names(data_dirs)  # Assign sample names

combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list)
)


# condition & Tissue Condition
cell_names <- colnames(combined_seurat)
condition <- rep(NA, length(cell_names))
ad_pattern <- grepl("_AD_", cell_names)
condition[ad_pattern] <- "AD"
control_pattern <- grepl("_Control_", cell_names)
condition[control_pattern] <- "Control"
combined_seurat$condition <- condition
combined_seurat$Tissue_Condition <- paste0(
  combined_seurat$orig.ident, 
  "_", 
  combined_seurat$condition
)
head(combined_seurat@meta.data)



saveRDS(combined_seurat, "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE175814_obj_seurat.rds")

##################### Study 4_1 (GSE163577) ##################### 
#Brain region:
#Mixed Between : Hippocampus and Superior frontal cortex


data_dirs <- list(
  "01_05_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/01_05_AD",
  "01_06_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/01_06_C",
  "01_09_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/01_09_C",
  "01_10_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/01_10_C",
  "02_01_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/02_01_C",
  "02_04_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/02_04_AD",
  "02_11_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/02_11_C",
  "03_02_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/03_02_AD",
  "03_03_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/03_03_C",
  "03_09_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/03_09_C",
  "04_09_C" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/04_09_C",
  "04_11_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/04_11_AD",
  "04_13_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/04_13_AD",
  "04_15_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/04_15_AD",
  "05_04_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/05_04_AD",
  "05_10_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/05_10_AD",
  "H2004_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/H2004_AD"
)


seurat_list <- lapply(names(data_dirs), function(sample_name) {
  data <- Read10X(data_dirs[[sample_name]])
  seurat_obj <- CreateSeuratObject(counts = data, project = "Hippocampus")
  seurat_obj$sample <- sample_name  # Add metadata to distinguish samples
  return(seurat_obj)
})
names(seurat_list) <- names(data_dirs)  # Assign sample names

combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list)
)

# condition & Tissue Condition
cell_names <- colnames(combined_seurat)
condition <- rep(NA, length(cell_names))
ad_pattern <- grepl("_AD_", cell_names)
condition[ad_pattern] <- "AD"
c_pattern <- grepl("_C_", cell_names)
condition[c_pattern] <- "Control"
combined_seurat$condition <- condition
combined_seurat$Tissue_Condition <- paste0(
  combined_seurat$orig.ident, 
  "_", 
  combined_seurat$condition
)
head(combined_seurat@meta.data)


# SAVE OBJ -- GSE163577 --
saveRDS(combined_seurat, "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE163577_HPC_obj_seurat.rds")

##################### Study 4_2 (GSE163577) ##################### 
#Brain region:
#Mixed Between : Hippocampus and Superior frontal cortex


data_dirs <- list(
  "O1_1OC_Ctx" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O1_1OC_Ctx",
  "O3_O3C_Ctx" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O3_O3C_Ctx",
  "04_15_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O3_O9C_Ctx",
  "05_04_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O4_11AD_Ctx",
  "05_10_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O4_13AD_Ctx",
  "H2004_AD" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O4_15AD_Ctx",
  "O1_1OC_Ctx" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O4_O9C_Ctx",
  "O3_O3C_Ctx" = "C:/Users/HaMMaDy/Downloads/Compressed/Grad Data/GSE163577/O5_1OAD_Ctx"
)


seurat_list <- lapply(names(data_dirs), function(sample_name) {
  data <- Read10X(data_dirs[[sample_name]])
  seurat_obj <- CreateSeuratObject(counts = data, project = "Superior frontal cortex")
  seurat_obj$sample <- sample_name  # Add metadata to distinguish samples
  return(seurat_obj)
})
names(seurat_list) <- names(data_dirs)  # Assign sample names

combined_seurat <- merge(
  x = seurat_list[[1]], 
  y = seurat_list[-1], 
  add.cell.ids = names(seurat_list)
)

# condition & Tissue Condition
cell_names <- colnames(combined_seurat)
condition <- rep(NA, length(cell_names))
ad_pattern <- grepl("_AD_", cell_names)
condition[ad_pattern] <- "AD"
control_pattern <- grepl("C_", cell_names)
condition[control_pattern] <- "Control"
combined_seurat$condition <- condition
combined_seurat$Tissue_Condition <- paste0(
  combined_seurat$orig.ident, 
  "_", 
  combined_seurat$condition
)
head(combined_seurat@meta.data)

# SAVE OBJ -- GSE163577 --
saveRDS(combined_seurat, "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE163577_CTX_obj_seurat.rds")


##################### Reading After Merge  ##################### 
Entorhinal_Cortex_Obj <- readRDS(file = "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE138852_obj_seurat.rds")
Prefrontal_Cortex_Obj <- readRDS(file = "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE157827_obj_seurat.rds")
Hippocampus_Obj <- readRDS(file = "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE175814_obj_seurat.rds")
HPC_Obj  <- readRDS(file = "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE163577_HPC_obj_seurat.rds")
CTX_Obj  <- readRDS(file = "C:/Users/HaMMaDy/Desktop/Grad_Obj/GSE163577_CTX_obj_seurat.rds")


seurat_list <- list(Entorhinal_Cortex_Obj, Prefrontal_Cortex_Obj, Hippocampus_Obj, HPC_Obj, CTX_Obj)

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],  # Merge all remaining studies
  add.cell.ids = names(seurat_list)  # Prefix cell names with study name
)

merged_metadata = merged_seurat@meta.data
metadata = merged_seurat@meta.data
new_idents = unique(metadata$orig.ident)
merged_seurat <- SetIdent(merged_seurat, value = metadata$orig.ident)

saveRDS(merged_seurat, "C:/Users/HaMMaDy/Desktop/Grad_Obj/All_Studies.rds")
