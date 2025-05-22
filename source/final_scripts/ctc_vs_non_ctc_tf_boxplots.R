


sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

sce$ctc <- ifelse(sce$new_clusters %in% c(8) & sce$condition == "cancer", "ctc", "non-ctc")

y <- assay(sce, "exprs")

dim(y)

df <- data.frame(t(y), colData(sce), check.names = FALSE)

value <- ifelse("exprs" == "exprs", "expression", "exprs")

gg_df <- melt(df, value.name = "expression", variable.name = "antigen", 
              id.vars = names(colData(sce)))

sclc_tfs <- c("NeuroD1","ASCL1","POU2F3")

temp <- gg_df %>% 
  dplyr::filter(antigen %in% sclc_tfs)

p <- ggboxplot(temp, x="ctc",y="expression",fill="ctc")

p + facet_wrap(~antigen, scales = "free_y")+
  stat_compare_means()
