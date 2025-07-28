source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")



markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam", "E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist", "SLUG", "PD-L1", "p-YAP", "CD44", "CD24")

scaled_heatmap <- create_expression_heatmap(ctcs, "subtype", markers_to_use,"", scale = T)


col_fun = colorRamp2(c(-3, -1, 0, 1, 3), 
                     c("#313695",  # deep blue
                       "#74add1",  # light blue
                       "#f7f7f7",  # white (center, 0)
                       "#f46d43",  # light red
                       "#a50026"))

ht <- Heatmap(t(scaled_heatmap),column_names_rot = 0,col = col_fun,
              cluster_columns = F, cluster_rows=F, column_title = "",
              row_names_gp = gpar(fontsize = 20),column_names_gp = gpar(fontsize = 20),
              heatmap_legend_param = list(
                title = "   Scaled\nExpression",      
                title_gp = gpar(fontsize = 18), 
                labels_gp = gpar(fontsize = 14),
                legend_height = unit(3, "cm"),
                grid_width = unit(.5,"cm")))

print(ht)

tiff("figures/ctcs_subtype_expression_heatmap.tiff", width=200,height=100, units = "mm", res=600)
print(ht)
dev.off()

# 
# # assay(ctcs, "exprs") <- t(scale(t(assay(ctcs, "exprs"))))
# 
# all_markers <- readRDS("data/state_markers.rds")
# 
# naive <- ctcs[,ctcs$treatment_status == "naive"]
# naive_ht <- create_expression_heatmap(naive, "subtype", all_markers)
# 
# treated <- ctcs[,ctcs$treatment_status == "treated"]
# treated_ht <- create_expression_heatmap(treated, "subtype", all_markers)
# 
# naive_ht+treated_ht
# 
# subtypes <- c("A","N","P","I")
# heatmaps <- list()
# for(curr_subtype in subtypes){
#   curr_sce <- ctcs[,ctcs$subtype == curr_subtype]
#   
#   # curr_sce <- curr_sce[-which(rownames(curr_sce) %in% c("CD24","Vimentin"))]
#   
#   curr_ht <- create_expression_heatmap(curr_sce, "treatment_status", rownames(curr_sce),curr_subtype,
#                                        scale = F)
#  
#   heatmaps <- append(heatmaps, list(curr_ht))
#    
# }
# 
# 
# 
# heatmaps[[1]]+heatmaps[[2]]+heatmaps[[3]]+heatmaps[[4]]
# # 
# # jpeg("figures/ctcs_subtype_expression_heatmap.jpg", width=300,height=100, units = "mm", res=1000)
# # print(ht)
# # dev.off()
# 
# ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")
# assay(ctcs, "exprs") <- t(scale(t(assay(ctcs, "exprs"))))
# curr_ht <- create_expression_heatmap(ctcs, "subtype", rownames(curr_sce), "all CTCs", scale=F)
# curr_ht
