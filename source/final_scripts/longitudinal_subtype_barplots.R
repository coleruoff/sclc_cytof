source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################

ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")


long_patients <- as.data.frame(ctcs@colData) %>% 
  select(patient_id,sample_num) %>% 
  distinct() %>% 
  dplyr::count(patient_id) %>% 
  dplyr::filter(n > 1) %>% 
  pull(patient_id) %>% 
  as.character()


long_data <- ctcs[,ctcs$patient_id %in% long_patients]


plot_df <- as.data.frame(long_data@colData) %>% 
  select(patient_id,sample_num,subtype) %>% 
  dplyr::count(patient_id,sample_num,subtype) %>% 
  group_by(patient_id,sample_num) %>% 
  mutate(total = sum(n)) %>% 
  mutate(freq = n/total)

plot_df$freq <- plot_df$freq*100


cluster_colors <- c(
  "#E57373",  # muted red
  "#64B5F6",  # muted blue
  "#81C784",  # muted green
  "#FFB74D"   # muted orange
)

p <- ggplot(plot_df)+
  geom_col(aes(x=sample_num,y=freq,fill=subtype))+
  facet_wrap(~patient_id, scales="free",ncol=5)+
  scale_fill_manual(name = "Subtype",values=cluster_colors)+
  labs(x="Sample Number",
       y="Percentage",
       fill="Subtype")+
  theme_classic()

p

jpeg("figures/longitudinal_subtype_barplots.jpg", width=180,height=100, units = "mm", res=1000)
print(p)
dev.off()


