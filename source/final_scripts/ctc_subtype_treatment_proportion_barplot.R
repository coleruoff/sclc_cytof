source("source/sclc_cytof_functions.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
cluster_colors <- c("#dd4b33", "#F1FAEE", "#A8DADC", "#457B9D")

ctcs <- readRDS("data/cytof_objects/ctcs_with_subtype.rds")

plot_df <-  ctcs@colData %>% 
  as.data.frame()  %>% 
  count(collection_id,subtype,treatment_status) %>% 
  group_by(collection_id) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = (n/total)*100) %>% 
  group_by(collection_id)


plot_df <- plot_df %>% 
  filter(total > 10)

sample_order <- plot_df %>% 
  filter(subtype == "I") %>% 
  arrange(desc(freq)) %>% 
  pull(collection_id)

plot_df$collection_id <- factor(plot_df$collection_id, levels = sample_order)

plot_df$subtype <- factor(plot_df$subtype, levels = c("A","N","P","I"))

plot_df$treatment_status <- ifelse(plot_df$treatment_status == "naive", "Naive","Treated")
plot_df$treatment_status <- factor(plot_df$treatment_status, levels = c("Naive","Treated"))

ggplot(plot_df)+
  geom_col(aes(x=collection_id,y=freq,fill=subtype))+
  geom_text(aes(label=total,x=collection_id,y=101.5),size=2)+
  # facet_wrap(~treatment_status,scales="free")+
  facet_grid(. ~ treatment_status, scales = "free", space = "free") +
  scale_fill_manual(values = cluster_colors)+
  theme_classic()+
  labs(x="Sample",
       y="Percentage",
       fill="Subtype")+
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
  


length(unique(plot_df$collection_id))
