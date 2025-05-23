source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

################################################################################
# Remove cell line samples
################################################################################
blood_samples <- as.data.frame(sce@colData) %>%
  dplyr::filter(sample_type == "blood") %>%
  pull(collection_id) %>%
  as.character()

sce <- sce[,sce$collection_id %in% blood_samples]

################################################################################

markers <- as.data.frame(rowData(sce)) %>%
  dplyr::filter(marker_class == "state") %>%
  pull(marker_name)

temp <- sce[markers,sce$collection_id == "NJH29-1"]

p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")

p1

################################################################################
# Perform batch correction

batch <- as.factor(colData(sce)$experiment_id)

design <- model.matrix(~ 0 + factor(colData(sce)$condition))  # one-hot encoding of conditions
colnames(design) <- levels(factor(colData(sce)$condition))

corrected_exprs <- removeBatchEffect(assay(sce, "exprs"), batch = batch, design = design)

assay(sce, "exprs") <- corrected_exprs

################################################################################

temp <- sce[markers,sce$collection_id == "NJH29-1"]

p1 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# temp <- sce_corrected[markers,sce_corrected$collection_id == "NJH29-1"]
# 
# p2 <- plotExprs(temp, color_by = "experiment_id", assay = "exprs")
# 
# p1+ggtitle("NJH29 (no batch correction)")
# p2+ggtitle("NJH29 (batch corrected)")
# 
# ################################################################################

# Get samples that are run in multiple experiments
to_test <- as.data.frame(sce@colData) %>%
  dplyr::select(experiment_id,collection_id) %>%
  distinct() %>%
dplyr::count(collection_id) %>%
  arrange(desc(n)) %>%
  dplyr::filter(n > 1) %>%
  pull(collection_id) %>%
  as.character()


sort(to_test)

# Checking SC454-1
temp <- sce[markers,sce$collection_id == "SC454-1"]

plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce_corrected[markers,sce_corrected$collection_id == "SC414-1"]
plotExprs(temp, color_by = "experiment_id",assay = "exprs")

################################################################################
# Remove outlier experiments
################################################################################
# sce <- sce[,sce$experiment_id != "531050"]
# sce <- sce[,sce$experiment_id != "508095"]
# sce <- sce[,sce$experiment_id != "515600"]

################################################################################
# Remove blood bank samples
################################################################################
# blood_bank_samples <- paste0("NORMAL", 7:20)
# 
# length(unique(sce$patient_id))
# 
# sce <- sce[,!sce$patient_id %in% blood_bank_samples]

################################################################################
# remove collections with < 30 cells
################################################################################
samples_to_remove <- as.data.frame(sce@colData) %>% 
  dplyr::count(collection_id) %>% 
  dplyr::filter(n<30) %>% 
  pull(collection_id) %>% 
  as.character()

sce <- sce[,!sce$collection_id %in% samples_to_remove]

################################################################################

saveRDS(sce, "data/cytof_objects/sclc_all_samples_object.rds")

