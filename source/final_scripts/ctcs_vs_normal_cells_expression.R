source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

sce <- readRDS("data/cytof_objects/sclc_all_samples_object.rds")

sce$ctc <- ifelse(sce$cell_id %in% ctcs$cell_id, "ctc","normal")

sce <- sce[,sce$condition == "normal" | sce$ctc == "ctc"]

################################################################################
# Plot CTCs vs normal cells

markers_to_use <- c("p-Rb")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
               id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$ctc <- ifelse(plot_df$ctc == "normal", "Normal","CTCs")
plot_df$ctc <- factor(plot_df$ctc, levels=c("CTCs","Normal"))

p <- ggviolin(plot_df, x="ctc" ,y="expression", fill="ctc",draw_quantiles = 0.5)+
  stat_compare_means(comparisons = list(c("CTCs","Normal")))+
  facet_wrap(~antigen)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

p

plot_df %>% 
  filter(ctc == "CTCs") %>% 
  pull(expression) %>% 
  summary()


jpeg(glue("figures/ctcs_vs_normal_rb_violin_plot.jpg"), width=100,height=100, units = "mm", res=1000)
print(p)
dev.off()

################################################################################
# Plot CTCs vs non-CTCs in cancer enriched clusters
sce <- readRDS("data/cytof_objects/cancer_enriched_with_clusters.rds")

sce$ctc <- ifelse(sce$new_clusters %in% c(2,4,5,6,7,8), "CTCs","Non-CTCs")

sce@colData %>% 
  as.data.frame() %>% 
  count(ctc)

markers_to_use <- c("p-Rb")

y <- assay(sce, "exprs")

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))


plot_df <- gg_df %>%
  dplyr::filter(antigen %in% markers_to_use)

plot_df$ctc <- factor(plot_df$ctc, levels=c("CTCs","Non-CTCs"))


p <- ggviolin(plot_df, x="ctc" ,y="expression", fill="ctc",draw_quantiles = 0.5)+
  stat_compare_means(comparisons = list(c("CTCs","Non-CTCs")))+
  facet_wrap(~antigen)+
  ylim(0,NA)+
  labs(y="Expression",
       x= "")+
  scale_fill_manual(name = "Subtype",values=c("#E63946","#457B9D"))+
  theme(axis.title = element_text(size=14),
        axis.text.x = element_text(size=10,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=12), 
        strip.background = element_blank())+
  rremove("legend")

plot_df %>% 
  filter(ctc == "CTCs") %>% 
  pull(expression) %>% 
  summary()

p


jpeg(glue("figures/ctcs_vs_nonctcs_rb_violin_plot.jpg"), width=100,height=100, units = "mm", res=1000)
print(p)
dev.off()
################################################################################
################################################################################
# 
# ggboxplot(plot_df, x="ctc",y="expression",fill="ctc",add = "jitter", add.params = list(size = 1, alpha = 0.025))+
#   facet_wrap(~antigen)+
#   stat_compare_means(label = "p.signif", label.y = c(7,11,7,9), tip.length = 0, size=8, comparisons = list(c("normal","ctc")))+
#   labs(y=glue("Expression"),
#        x="")+
#   ylim(c(0,12))+
#   rremove("legend")+
#   theme(axis.text = element_text(size=20),
#         axis.title = element_text(size=20))
# 
# 
# p1 <- p + facet_wrap(~antigen, scales = "free_y")+
#   stat_compare_means()
# 
# p1
# 
# 
# 
# 
# rows_to_keep <- intersect(rownames(sce),rownames(ctcs))
# sce <- sce[rows_to_keep,]
# 
# 
# ctcs$cell_id
# sce$cell_id
# sce$sample_id
# 
# 
# 
# ctcs <- ctcs[,ctcs$treatment_status == "naive"]
# normal_cells <- sce[,sce$condition == "normal"]
# 
# 
# curr_marker <- "Alcam"
# 
# 
# 
# plot_df <- as.data.frame(rbind(cbind(ctcs@assays@data$exprs[curr_marker,],"CTCs"),cbind(normal_cells@assays@data$exprs[curr_marker,],"Normal Cells")))
# 
# 
# colnames(plot_df) <- c("value","class")
# 
# plot_df$value <- as.numeric(plot_df$value)
# 
# ggviolin(plot_df, x="class",y="value",fill="class",add = "none", add.params = list(size = 1, alpha = 0.3))+
#   stat_compare_means(label = "p.signif", label.y = 6, tip.length = 0, size=8, comparisons = list(c("CTCs","Normal Cells")))+
#   labs(y=glue("{curr_marker} Expression"),
#        x="")+
#   ylim(c(0,6.5))+
#   rremove("legend")+
#   theme(axis.text = element_text(size=20),
#         axis.title = element_text(size=20))
# 
# 
# 
# 
#   geom_beeswarm(binaxis = "y",
#               binwidth = .05,
#               stackdir = "center",
#               dotsize = 0.01)
# 
# 
# 
# 
#   
#   
#   y <- assay(ctcs, "exprs")
#   
#   df <- data.frame(t(y), colData(ctcs), check.names = FALSE)
#   
#   value <- ifelse("exprs" == "exprs", "expression", "exprs")
#   
#   gg_df2 <- melt(df, value.name = "expression", variable.name = "antigen", 
#                  id.vars = names(colData(ctcs)))
#   
#   dim(gg_df1)
#   dim(gg_df2)
#   
#   
#   
#   ################################################################################
#   
#   plot_df <- gg_df %>%
#     dplyr::filter(antigen %in% markers_to_use)
#   
#   
#   if(is.null(fill)){
#     p <- ggboxplot(plot_df, x=group ,y="expression", fill=group)
#   } else{
#     p <- ggboxplot(plot_df, x=group ,y="expression", fill=fill)
#   }
#   
#   
#   p1 <- p + facet_wrap(~antigen, scales = "free_y")+
#     stat_compare_means()
#   
# 
# 
