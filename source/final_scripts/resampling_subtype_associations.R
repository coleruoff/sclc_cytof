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

for(n_cells in c(20,50)){
  
  treatment_results <- list()
  tarla_results <- list()
  for(i in 1:1000){
    sampled_data <- ctcs@colData %>% 
      as.data.frame() %>%
      group_by(patient_id) %>%
      filter(n() >= n_cells) %>%        # keep only patients with ≥ n_cells cells
      slice_sample(n = n_cells) %>%
      ungroup()
    
    ################################################################################
    # Treatment Status Association
    ################################################################################
    data_df <- sampled_data
    
    data_df$treatment_status <- factor(data_df$treatment_status, levels = c("naive","treated"))
    
    results_list <- list()
    for(curr_subtype in c("A","N","P","I")){
      data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
      
      formula_str <- glue("curr_subtype ~ treatment_status + (1 | patient_id)")
      
      model <- glmer(
        formula = as.formula(formula_str),
        family = binomial(link = "logit"),
        data = data_df)
      
      or <- exp(fixef(model)[2])
      
      tidy_out <- tidy(model,effects='fixed')
      curr_pval <- tidy_out$p.value[2]
      
      lower_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][1])
      upper_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][2])
      
      res <- data.frame("subtype"=curr_subtype,"or"=or,"pval"=curr_pval,up_or=upper_or,low_or=lower_or)
      
      results_list <- append(results_list, list(res))
      
    }
    
    # Combine into one data frame
    all_results <- bind_rows(results_list) %>% 
      mutate(iter = i)
    
    treatment_results <- append(treatment_results,list(all_results))
    
    ################################################################################
    # Tarla Status Association
    ################################################################################
    data_df <- sampled_data
    
    data_df$tarla <- factor(data_df$tarla, levels = c("pre","post"))
    
    results_list <- list()
    for(curr_subtype in c("A","N","P","I")){
      data_df$curr_subtype <- as.factor(as.integer(data_df$subtype == curr_subtype))
      
      formula_str <- glue("curr_subtype ~ tarla + (1 | patient_id)")
      
      model <- glmer(
        formula = as.formula(formula_str),
        family = binomial(link = "logit"),
        data = data_df)
      
      or <- exp(fixef(model)[2])
      
      tidy_out <- tidy(model,effects='fixed')
      curr_pval <- tidy_out$p.value[2]
      
      lower_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][1])
      upper_or <- exp(confint(model, parm = "beta_", method = "Wald")[2,][2])
      
      res <- data.frame("subtype"=curr_subtype,"or"=or,"pval"=curr_pval,up_or=upper_or,low_or=lower_or)
      
      results_list <- append(results_list, list(res))
      
    }
    
    # Combine into one data frame
    all_results <- bind_rows(results_list)
    
    # Combine into one data frame
    all_results <- bind_rows(results_list) %>% 
      mutate(iter = i)
    
    tarla_results <- append(tarla_results,list(all_results))
    
  }
  
  all_res <- bind_rows(treatment_results)
  
  wilcox_summary <- all_res %>%
    group_by(subtype) %>%
    summarise(
      median_logOR = median(or, na.rm = TRUE),
      lower_CI = quantile(or, 0.025, na.rm = TRUE),
      upper_CI = quantile(or, 0.975, na.rm = TRUE),
      pval = wilcox.test(log(or), mu = 0)$p.value
    )
  
  wilcox_summary <- wilcox_summary %>%
    mutate(padj = p.adjust(pval, method = "BH"))
  
  summary_res <-all_res %>%
    group_by(subtype) %>%
    summarise(
      median_logOR = median(or, na.rm = TRUE),
      lower_CI = quantile(or, 0.025, na.rm = TRUE),
      upper_CI = quantile(or, 0.975, na.rm = TRUE),
      prop_sig = mean(pval < 0.05, na.rm = TRUE)
    )
  
  all_res$subtype <- factor(all_res$subtype,levels = c("A","N","P",'I'))
  
  p1 <- ggplot(all_res, aes(x = log(or), y = fct_rev(subtype))) +
    geom_violin(aes(fill=subtype), alpha = .8, trim = FALSE, draw_quantiles = 0.5) +
    # geom_jitter(width = 0.2, alpha = 0.3, size = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    xlim(-2,2)+
    annotate("text", x=-1.5, y=4.5, label = "Naive", angle=0,size=6) +
    annotate("text", x=1.5, y=4.5, label = "SOC", angle=0,size=6) +
    theme_classic()+
    scale_fill_manual(values = cluster_colors)+
    labs(
      y = "Subtype",
      x = "log(OR)")+
    theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
          axis.title = element_text(size=24),
          axis.text.x = element_text(angle = 0, hjust = .5))+
    rremove("legend")
  
  all_res <- bind_rows(tarla_results)
  
  summary_res <-all_res %>%
    group_by(subtype) %>%
    summarise(
      median_logOR = median(or, na.rm = TRUE),
      lower_CI = quantile(or, 0.025, na.rm = TRUE),
      upper_CI = quantile(or, 0.975, na.rm = TRUE),
      prop_sig = mean(pval < 0.05, na.rm = TRUE)
    )
  
  all_res$subtype <- factor(all_res$subtype,levels = c("A","N","P",'I'))
  
  p2 <- ggplot(all_res, aes(x = log(or), y = fct_rev(subtype)))+
    geom_violin(aes(fill=subtype), alpha = .8, trim = FALSE, draw_quantiles = 0.5) +
    # geom_jitter(width = 0.2, alpha = 0.3, size = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    xlim(-3,3)+
    annotate("text", x=-1.5, y=4.5, label = "Pre-Tarlatamab", angle=0,size=6) +
    annotate("text", x=1.5, y=4.5, label = "Post-Tarlatamab", angle=0,size=6) +
    theme_classic()+
    scale_fill_manual(values = cluster_colors)+
    labs(
      y = "Subtype",
      x = "log(OR)")+
    theme(axis.text = element_text(size=22,angle = 0, hjust = 1),
          axis.title = element_text(size=24),
          axis.text.x = element_text(angle = 0, hjust = .5))+
    rremove("legend")
  
  # p2
  
  ################################################################################
  # Save figures
  ################################################################################
  
  
  tiff(glue("figures/downsampled_subtype_treatment_status_or_results_{n_cells}.tiff"), width=200,height=200, units = "mm", res=600)
  print(p1)
  dev.off()
  
  tiff(glue("figures/downsampled_subtype_tarla_or_results_{n_cells}.tiff"), width=200, height=200, units = "mm", res=600)
  print(p2)
  dev.off()
}









