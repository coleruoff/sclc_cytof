################################################################################
# Sample level (non-longitudinal patients) alluvial
################################################################################
# Find patient level subtypes
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

# Count cells
samples_used <- as.character(patient_level_subtype$collection_id)

count_df <- ctcs@colData %>% 
  as.data.frame() %>% 
  filter(collection_id %in% samples_used) 

count_df$treatment_status <- ifelse(count_df$treatment_status == "naive","Naive","SOC")
count_df$tarla <- ifelse(count_df$tarla == "pre","Pre-Tarla","Tarla")

count_df$treatment_status <- ifelse(count_df$treatment_status == "Naive", count_df$treatment_status, ifelse(count_df$tarla == "Pre-Tarla" | is.na(count_df$tarla),count_df$treatment_status,"Tarla"))

count_df$treatment_status <- factor(count_df$treatment_status, levels=c("Naive","SOC","Pre-Tarla","Tarla"))

count_df %>% 
  count(treatment_status)

count_df %>% 
  count(subtype)

#######

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


contin_table <- plot_df %>% 
  count(subtype,treatment_status) %>% 
  arrange(subtype) %>% 
  pivot_wider(names_from = treatment_status,values_from = n) %>% 
  column_to_rownames("subtype") 


contin_table[is.na(contin_table)] <- 0

fisher.test(contin_table)


treatment_table <- contin_table[,1:2]

subtypes <- rownames(treatment_table)

results_list <- list()

for(i in 1:4){
  
  contin_table <- rbind(treatment_table[i,],colSums(treatment_table[-i,]))
  
  fisher_res <- fisher.test((contin_table+1))
  
  results_list[["subtype"]] <- append(results_list[["subtype"]], subtypes[i])
  results_list[["or"]] <- append(results_list[["or"]], fisher_res$estimate)
  results_list[["pval"]] <- append(results_list[["pval"]], fisher_res$p.value)
  results_list[["up_or"]] <- append(results_list[["up_or"]], fisher_res$conf.int[1])
  results_list[["low_or"]] <- append(results_list[["low_or"]], fisher_res$conf.int[2])
  
}


plot_df <- as.data.frame(results_list)

plot_df <- plot_df %>% 
  mutate(padj = p.adjust(pval), method = "BH") %>% 
  mutate(signif = ifelse(padj < 0.05, "s","ns")) %>% 
  mutate(log_or = log(or)) %>% 
  mutate(log_upper_or = log(up_or)) %>% 
  mutate(log_lower_or = log(low_or)) %>% 
  mutate(star_height = ifelse(log_or > 0, log_or+.1,log_or-.1))

ggplot(plot_df,aes(x=log_or,y=fct_rev(subtype),color=subtype))+
  geom_point(aes(shape = factor(signif)),size=12,fill="white",show.legend = F, stroke=3)+
  scale_shape_manual(values = c("ns" = 1, "s" = 16)) +
  geom_errorbarh(aes(xmin = log_lower_or, xmax = log_upper_or), height = 0.1, linewidth = 1,show.legend = F)+
  geom_vline(xintercept = 0, linetype = 2)+
  scale_color_manual(values = cluster_colors)+
  xlim(-5.75,5.75)+
  labs(y="Subtype",
       x="log(OR)")+
  theme_classic()+
  annotate("text", x=-.85, y=4.5, label = "SOC", angle=0,size=8) +
  annotate("text", x=.85, y=4.5, label = "Naive", angle=0,size=8) +
  theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
        axis.title = element_text(size=24),
        axis.text.x = element_text(angle = 0, hjust = .5))



