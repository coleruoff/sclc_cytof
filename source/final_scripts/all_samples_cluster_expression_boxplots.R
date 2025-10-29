################################################################################
# This script generates boxplots for selected proteins displaying the 
# distribution of expression across all cells in the given comparison groups
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

cancer_enriched_clusters <- readRDS("data/cancer_enriched_clusters.rds")
################################################################################
# Create cluster boxplots using SCLC subtype TFs.
################################################################################

sce$cancer_enriched <- ifelse(sce$new_clusters %in% cancer_enriched_clusters, "cancer_enriched", "non-cancer_enriched")

markers_to_use <- c("NeuroD1","ASCL1","POU2F3","p-Rb")

p1 <- create_marker_boxplots(sce, markers_to_use, group="new_clusters",fill="cancer_enriched")


p1
p1 <- p1 +
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  rremove("legend")

################################################################################
# Create cancer-enriched vs normal enriched boxplots using SCLC subtype TFs 
################################################################################
p2 <- create_marker_boxplots(sce, markers_to_use, "cancer_enriched","cancer_enriched")

p2 <- p2+
  stat_compare_means(label = "p.signif", tip.length = 0, comparisons = list(c("cancer_enriched","non-cancer_enriched")))+
  labs(y="Expression",
       x="",
       fill="")+
  ylim(0,10)+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  scale_x_discrete(labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Normal Enriched"))+
  rremove("legend")

################################################################################
# Create cluster boxplots using selected markers
################################################################################
markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

p3 <- create_marker_boxplots(sce, markers_to_use, "new_clusters","cancer_enriched")

p3 <- p3+
  labs(y="Expression",
       x="",
       fill="")+
  scale_fill_manual(
    values = c("cancer_enriched" = "#E57373", "non-cancer_enriched" = "#64B5F6"),
    labels = c("cancer_enriched" = "Cancer Enriched", "non-cancer_enriched" = "Non-Cancer Enriched"))+
  rremove("legend")

################################################################################
# Save figures
################################################################################
tiff(glue("figures/all_samples_cluster_tf_expression_boxplot.tiff"), width=160,height=160, units = "mm", res=1000)
print(p1)
dev.off()

tiff(glue("figures/all_samples_ce_vs_nonce_tf_expression_boxplot.tiff"), width=160,height=160, units = "mm", res=1000)
print(p2)
dev.off()

tiff(glue("figures/all_samples_cluster_expression_boxplot.tiff"), width=320,height=200, units = "mm", res=1000)
print(p3)
dev.off()
