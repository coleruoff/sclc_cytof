source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")
################################################################################
# Cell level
################################################################################

plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1")



plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))


to_lodes_form(data.frame(plot_df),
              key = "category", value = "group", id = "cohort",
              axes = 1:2) %>% 
  add_count(group)


plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 







p1 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=10) +
  geom_text(stat = "stratum", aes(label = glue("\U00A0\n\nn={nn} ({freq}%)")),size=5) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p1

tiff("figures/cell_level_alluvial_plot.tiff", width=300,height=400, units = "mm", res=600)
print(p1)
dev.off()
################################################################################
# Sample level
################################################################################
samples_to_use <- ctcs@colData %>% 
  as.data.frame() %>% 
  count(collection_id) %>% 
  filter(n > 10) %>% 
  pull(collection_id) %>% 
  as.character() %>% 
  unique()


sample_level_subtype <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(collection_id,treatment_status,subtype) %>% 
  group_by(treatment_status) %>% 
  count(collection_id,subtype) %>% 
  group_by(collection_id,treatment_status) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) 



plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% samples_to_use) %>% 
  select(collection_id,treatment_status,tarla) %>%
  distinct() %>% 
  merge(.,sample_level_subtype,by="collection_id") %>% 
  select(treatment_status.x,subtype,tarla) %>% 
  rename("treatment_status" = "treatment_status.x") %>% 
  cbind(1) %>% 
  rename("n" = "1") 


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))



plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) %>% 
  mutate(stats = ifelse(group == "A", "", glue("\U00A0\n\n\nn={nn} ({freq}%)"))) %>% 
  mutate(a_stats = ifelse(group == "A",glue("n={nn} ({freq}%)"),""))



p2 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total, label=stats)) +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=10) +
  geom_text(stat = "stratum", aes(label = stats),size=5) +
  # ggfittext::geom_fit_text(stat = "stratum", width = .5, min.size = 10) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(group == "A", as.character(a_stats), NA)),
    stat = "stratum", size = 5, direction = "y", nudge_x = .3) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p2

tiff("figures/sample_level_alluvial_plot.tiff", width=300,height=500, units = "mm", res=600)
print(p2)
dev.off()

################################################################################
# Sample level (non-longitudinal patient)
################################################################################
samples_to_use <- ctcs@colData %>% 
  as.data.frame() %>% 
  count(collection_id) %>% 
  filter(n > 10) %>% 
  pull(collection_id) %>% 
  as.character() %>% 
  unique()

non_long_patients <- ctcs@colData %>%
  as.data.frame() %>% 
  filter(collection_id %in% samples_to_use) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  count(patient_id) %>% 
  filter(n == 1) %>% 
  pull(patient_id) %>% 
  as.character() 


patient_level_subtype <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(patient_id %in% non_long_patients & collection_id %in% samples_to_use) %>% 
  select(collection_id,subtype) %>% 
  count(collection_id,subtype) %>% 
  group_by(collection_id) %>% 
  filter(n == max(n)) %>% 
  select(collection_id,subtype) %>% 
  distinct()


tied_patients <- patient_level_subtype %>% 
  count(collection_id) %>% 
  filter(n > 1)



plot_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  select(collection_id,treatment_status,tarla) %>%
  distinct() %>% 
  merge(.,patient_level_subtype,by="collection_id") %>% 
  select(treatment_status,subtype,tarla) %>% 
  cbind(1) %>% 
  rename("n" = "1") 


plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive","Naive","SOC")
plot_df$tarla <- ifelse(plot_df$tarla == "pre","Pre-Tarla","Tarla")

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "Naive", plot_df$treatment_status, ifelse(plot_df$tarla == "Pre-Tarla" | is.na(plot_df$tarla),plot_df$treatment_status,"Tarla"))

plot_df$treatment_status <- factor(plot_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

plot_df$subtype <- factor(plot_df$subtype, levels=c("A","N","P",'I'))



plot_df_long <- to_lodes_form(data.frame(plot_df),
                              key = "category", value = "group", id = "cohort",
                              axes = 1:2) %>% 
  add_count(group) %>% 
  group_by(category) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=sprintf("%.1f",(nn/total)*100)) 


p3 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=10) +
  geom_text(stat = "stratum", aes(label = glue("\U00A0\n\nn={nn} ({freq}%)")),size=5) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p3

tiff("figures/patient_level_alluvial_plot.tiff", width=300,height=500, units = "mm", res=600)
print(p3)
dev.off()

