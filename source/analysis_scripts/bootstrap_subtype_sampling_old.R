source("source/sclc_cytof_functions.R")

set.seed(42)

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

# cluster_colors <- c("#b54a3a","#E6EAE3","#8FB6B8","#3D617C")


subtype_status_colors <- c("#b54a3a","#dd4b33","#E6EAE3", "#F1FAEE","#8FB6B8","#A8DADC","#3D617C","#3D617C")



all_data <- ctcs@colData %>% 
  as.data.frame()

all_data$treatment_status <- ifelse(all_data$treatment_status == "naive","Naive","SOC")

all_data$treatment_status <- ifelse(is.na(all_data$tarla), all_data$treatment_status,
                                    ifelse(all_data$tarla == "pre", all_data$treatment_status, "Tarla"))


################################################################################
# Calculate mean subtype proportions
# Naive vs SOC
################################################################################
n_cells <- 30

data_df <- all_data %>% 
  filter(treatment_status %in% c("Naive","SOC"))

# Repeat 1000 times
resamples <- map(1:1000, ~ {
  all_data %>% 
    group_by(collection_id) %>%
    filter(n() >= n_cells) %>%
    slice_sample(n = n_cells) %>%
    count(subtype,treatment_status) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

df_resampled$treatment_status <- factor(df_resampled$treatment_status, levels = c("Naive","SOC"))
df_resampled$subtype <- factor(df_resampled$subtype, levels = c("A","N","P","I"))
df_resampled$treatment_subtype <- paste0(df_resampled$treatment_status,"_",df_resampled$subtype)
df_resampled$treatment_subtype <- factor(df_resampled$treatment_subtype, levels=c("Naive_A","SOC_A","Naive_N","SOC_N","Naive_P","SOC_P","Naive_I","SOC_I"))


stat.test <- df_resampled %>%
  group_by(subtype) %>%
  wilcox_test(prop ~ treatment_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "treatment_status")

p <- ggviolin(df_resampled, x="treatment_status",y="prop",fill="treatment_subtype",draw_quantiles = 0.5)+
  facet_wrap(~subtype, scales="free_x",nrow=1)+
  stat_pvalue_manual(stat.test, label = "p.adj.signif",size=5,y.position = 1,tip.length = .01)+
  scale_fill_manual(values = subtype_status_colors)+
  labs(y="Proportion",
       x= "Treatment Status")+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank())+
  rremove("legend")
  
p

tiff(glue("figures/resampling_naive_soc_proportion_violinplot.tiff"), width=300,height=200, units = "mm", res=600)
print(p)
dev.off()

naive_soc_pvals <- stat.test$p.adj
sprintf("%.4f",stat.test$p.adj)

# Global consensus with variability
global_consensus <- df_resampled %>%
  group_by(subtype,treatment_status) %>%
  summarise(
    mean_prop = mean(prop),
    sd_prop   = sd(prop),
    lower_ci  = quantile(prop, 0.025),
    upper_ci  = quantile(prop, 0.975),
    .groups = "drop"
  )

global_consensus

################################################################################
# Naive vs Tarla
################################################################################
data_df <- all_data %>% 
  filter(treatment_status %in% c("Naive","Tarla"))

# Repeat 1000 times
resamples <- map(1:1000, ~ {
  data_df %>% 
    group_by(collection_id) %>%
    filter(n() >= n_cells) %>%
    slice_sample(n = n_cells) %>%
    count(subtype,treatment_status) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

df_resampled$treatment_status <- factor(df_resampled$treatment_status, levels = c("Naive","Tarla"))
df_resampled$subtype <- factor(df_resampled$subtype, levels = c("A","N","P","I"))
df_resampled$treatment_subtype <- paste0(df_resampled$treatment_status,"_",df_resampled$subtype)
df_resampled$treatment_subtype <- factor(df_resampled$treatment_subtype, levels=c("Naive_A","Tarla_A","Naive_N","Tarla_N","Naive_P","Tarla_P","Naive_I","Tarla_I"))


stat.test <- df_resampled %>%
  group_by(subtype) %>%
  wilcox_test(prop ~ treatment_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "treatment_status")

p <- ggviolin(df_resampled, x="treatment_status",y="prop",fill="treatment_subtype",draw_quantiles = 0.5)+
  facet_wrap(~subtype, scales="free_x",nrow=1)+
  stat_pvalue_manual(stat.test, label = "p.adj.signif",size=5,y.position = 1,tip.length = .01)+
  scale_fill_manual(values = subtype_status_colors)+
  labs(y="Proportion",
       x= "Treatment Status")+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank())+
  rremove("legend")

p

tiff(glue("figures/resampling_naive_tarla_proportion_violinplot.tiff"), width=300,height=200, units = "mm", res=600)
print(p)
dev.off()

naive_tarla_pvals <- stat.test$p.adj


sprintf("%.4f",stat.test$p.adj)
# Global consensus with variability
# global_consensus <- df_resampled %>%
#   group_by(subtype,tarla) %>%
#   summarise(
#     mean_prop = mean(prop),
#     sd_prop   = sd(prop),
#     lower_ci  = quantile(prop, 0.025),
#     upper_ci  = quantile(prop, 0.975),
#     .groups = "drop"
#   )
# 
# global_consensus
###########
################################################################################
# Naive vs Tarla
################################################################################
data_df <- all_data %>% 
  filter(treatment_status %in% c("SOC","Tarla"))

# Repeat 1000 times
resamples <- map(1:1000, ~ {
  data_df %>% 
    group_by(collection_id) %>%
    filter(n() >= n_cells) %>%
    slice_sample(n = n_cells) %>%
    count(subtype,treatment_status) %>%                   # count subtypes across ALL patients
    mutate(prop = n / sum(n), 
           iteration = .x) %>% 
    ungroup()
})

# Combine
df_resampled <- bind_rows(resamples)

df_resampled$treatment_status <- factor(df_resampled$treatment_status, levels = c("SOC","Tarla"))
df_resampled$subtype <- factor(df_resampled$subtype, levels = c("A","N","P","I"))
df_resampled$treatment_subtype <- paste0(df_resampled$treatment_status,"_",df_resampled$subtype)
df_resampled$treatment_subtype <- factor(df_resampled$treatment_subtype, levels=c("SOC_A","Tarla_A","SOC_N","Tarla_N","SOC_P","Tarla_P","SOC_I","Tarla_I"))


stat.test <- df_resampled %>%
  group_by(subtype) %>%
  wilcox_test(prop ~ treatment_status) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>% 
  add_xy_position(x = "treatment_status")

p <- ggviolin(df_resampled, x="treatment_status",y="prop",fill="treatment_subtype",draw_quantiles = 0.5)+
  facet_wrap(~subtype, scales="free_x",nrow=1)+
  stat_pvalue_manual(stat.test, label = "p.adj.signif",size=5,y.position = 1,tip.length = .01)+
  scale_fill_manual(values = subtype_status_colors)+
  labs(y="Proportion",
       x= "Treatment Status")+
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=16,angle=45,hjust=1,vjust=1),
        strip.text = element_text(face = "bold", size=16), 
        strip.background = element_blank())+
  rremove("legend")

p

tiff(glue("figures/resampling_soc_tarla_proportion_violinplot.tiff"), width=300,height=200, units = "mm", res=600)
print(p)
dev.off()

soc_tarla_pvals <- stat.test$p.adj


# sprintf("%.4f",stat.test$p.adj)
# # Global consensus with variability
# global_consensus <- df_resampled %>%
#   group_by(subtype,tarla) %>%
#   summarise(
#     mean_prop = mean(prop),
#     sd_prop   = sd(prop),
#     lower_ci  = quantile(prop, 0.025),
#     upper_ci  = quantile(prop, 0.975),
#     .groups = "drop"
#   )
# 
# global_consensus


naive_soc_pvals
naive_tarla_pvals
soc_tarla_pvals
