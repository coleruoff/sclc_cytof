################################################################################
# This script downsamples cells from patients to n_cells and then plots the 
# proportions of subtypes and treatment status as an alluvial plot.
################################################################################
source("source/sclc_cytof_functions.R")

set.seed(42)
################################################################################
# Read in data
################################################################################
ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

cluster_colors <- c("#dd4b33", "#D1DACF", "#A8DADC", "#457B9D")

################################################################################
# # Sample n cells from each patient
################################################################################
n_cells <- 50

sampled_data <- ctcs@colData %>% 
  as.data.frame() %>%
  group_by(patient_id) %>%
  filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
  slice_sample(n = n_cells) %>%
  ungroup()

length(unique(sampled_data$patient_id))

################################################################################
# Create plot dataframe
################################################################################
plot_df <- sampled_data %>% 
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


plot_df_long_left <- subset(plot_df_long, group %in% c("Naive","SOC","Tarla"))
plot_df_long_right <- subset(plot_df_long, !group %in% c("Naive","SOC","Tarla"))

################################################################################
# Plot alluvial 
################################################################################
p1 <- ggplot(data = plot_df_long,
             aes(x = category, stratum = group, alluvium = cohort, y = total)) +
  coord_flip() +
  scale_y_reverse() +
  geom_flow(aes(fill=group),width=.3, aes.flow = "backward") +
  geom_stratum(aes(fill=group),width=.3) +
  geom_text(stat = "stratum", aes(label = glue("{group}")),size=12) +
  ggrepel::geom_text_repel(data=plot_df_long_left,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = -.3,
                           size=10,segment.color = NA) +
  ggrepel::geom_text_repel(data=plot_df_long_right,stat = "stratum", aes(label = glue("n={nn}\n({freq}%)")),nudge_x = .3,
                           size=10,segment.color = NA) +
  scale_fill_manual(name = "group",values=c("gray90","gray8 0","gray70",cluster_colors))+
  theme_void() +
  rremove("legend")

p1
################################################################################
# Save figure
################################################################################

tiff(glue("figures/cell_level_alluvial_plot_downsampled_{n_cells}.tiff"), width=500,height=300, units = "mm", res=600)
print(p1)
dev.off()

################################################################################
# Calculate percentage of subtypes in each treatment status
################################################################################
plot_df %>% 
  count(subtype,treatment_status) %>% 
  group_by(treatment_status) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  select(treatment_status,subtype,freq) %>% 
  arrange(treatment_status)

