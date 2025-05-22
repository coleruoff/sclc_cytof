source("source/sclc_cytof_functions.R")
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_object_no_qc.rds")

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
temp <- sce[,sce$collection_id == "SC454-1"]

plotExprs(temp, color_by = "experiment_id", assay = "exprs")

temp <- sce[,sce$collection_id == "SC414-1"]
plotExprs(temp, color_by = "experiment_id",assay = "exprs")

temp <- sce[,sce$collection_id == "SC443-1"]

plotExprs(temp, color_by = "experiment_id", assay = "exprs")

################################################################################
#remove experiment 531050
################################################################################
# sce <- sce[,sce$experiment_id != "531050"]

# sce <- sce[,sce$patient_id != "SC443"]

################################################################################
# remove collections with < 30 cells
################################################################################
samples_to_remove <- as.data.frame(sce@colData) %>% 
  dplyr::count(collection_id) %>% 
  dplyr::filter(n<30) %>% 
  pull(collection_id) %>% 
  as.character()

sce <- sce[,!sce$collection_id %in% samples_to_remove]


saveRDS(sce, "data/cytof_objects/sclc_all_samples_object.rds")

