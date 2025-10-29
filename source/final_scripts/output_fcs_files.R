source("source/sclc_cytof_functions.R")
################################################################################
# Output all CTCs
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

fcs <- CATALYST::sce2fcs(ctcs, split_by = "collection_id")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Subtypes
################################################################################
fcs <- CATALYST::sce2fcs(ctcs, split_by = "subtype")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Naive and Treated files
################################################################################

# Subset to remove post-tarla cells
treatment_status_ctcs <- ctcs[,is.na(ctcs$tarla) | ctcs$tarla != "post"]

fcs <- CATALYST::sce2fcs(treatment_status_ctcs, split_by = "treatment_status")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Output Pre- and Post-Tarla Files
################################################################################

# Subset to only tarla cells
tarla_ctcs <- ctcs[,!is.na(ctcs$tarla)]

fcs <- CATALYST::sce2fcs(tarla_ctcs, split_by = "tarla")

write.flowSet(fcs, outdir = "data/output_fcs_files")

################################################################################
# Initial Clusters 1-5
################################################################################

sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

clusters_1_to_5 <- sce[,sce$cluster_id %in% c(1,2,3,4,5)]

fcs <- CATALYST::sce2fcs(clusters_1_to_5)

write.FCS(fcs, filename = "data/output_fcs_files/clusters_1_to_5.fcs")
################################################################################
# Second clustering 1 and 3
################################################################################

sce <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

clusters_1_and_3 <- sce[,sce$cluster_id %in% c(1,3)]

fcs <- CATALYST::sce2fcs(clusters_1_and_3)

write.FCS(fcs, filename = "data/output_fcs_files/clusters_1_and_3.fcs")

################################################################################
# Second clustering 2,4,5,6,7,8
################################################################################

sce <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

ctc_clusters <- sce[,!sce$cluster_id %in% c(1,3)]

fcs <- CATALYST::sce2fcs(ctc_clusters)

write.FCS(fcs, filename = "data/output_fcs_files/ctc_clusters")
