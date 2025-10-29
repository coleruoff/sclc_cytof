
################################################################################
# This script plots violin plots of the expression of all protein markers between
# naive CTCs and CTCs treated with SOC
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")


curr_data <- ctcs[,ctcs$patient_id == "SC547"]
################################################################################
# Create plot dataframe
################################################################################
# markers_to_use <- c("ASCL1", "NeuroD1", "POU2F3", "DLL3", "Alcam","SLUG", "PD-L1", "p-YAP", "CD44", "CD24","E-Cad", "EpCAM", "MUC-1", "Vimentin", "Twist")
markers_to_use <- c("DLL3")

y <- assay(curr_data, "exprs")

df <- data.frame(t(y), colData(curr_data), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(curr_data)))

plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$antigen <- factor(plot_df$antigen, levels = markers_to_use)


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","CTX + ICI")

plot_df$treatment_status <- ifelse(is.na(plot_df$tarla), plot_df$treatment_status,
                                   ifelse(plot_df$tarla == "pre", plot_df$treatment_status, "Tarla"))

# plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive", "SOC")

plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive", "CTX + ICI","Tarla"))

# Remove post-tarla samples
# plot_df <- plot_df %>%
#   filter(tarla != "post" | is.na(tarla))

################################################################################
# Plot violin plots
################################################################################

plot_df$collection_id
stat.test <- plot_df %>%
  wilcox_test(expression ~ collection_id, comparisons = list(c("SC547-1","SC547-2"),c("SC547-1","SC547-3"),c("SC547-2","SC547-3"))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "collection_id")

stat.test

ggviolin(plot_df, x="collection_id" ,y="expression", fill="#457B9D", lwd=.3, outlier.size = .1,draw_quantiles =0.5,)+
  stat_pvalue_manual(stat.test, y.position = c(6.2,7,6.5),label = "p.adj.signif",size=5,tip.length = 0)+
  facet_wrap(~antigen,nrow=2)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  # scale_fill_manual(name = "Treatment",values=c("#457B9D","#A8DADC"))+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        axis.text.y = element_text(size=16),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank(),
        legend.position = "none")





