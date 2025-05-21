source("source/cytof_de_function.R")

script_seed <- 42
set.seed(script_seed)
################################################################################
ctcs <- readRDS("data/cytof_objects/all_samples_ctcs_with_subtype.rds")

################################################################################
# Treatment status and subtype association
################################################################################
contin_table <- table(ctcs$treatment_status, ctcs$subtype) %>% 
  as.data.frame.matrix()

# colnames(contin_table) <- c("ASCL1","Inflamed","NeuroD1","POU2F3")


contin_table <- contin_table[c(1,3),]

test_res <- chisq.test(contin_table)

residuals <- test_res$stdres

# Calculate two-sided p-values from standard normal distribution
pvals <- 2 * (1 - pnorm(abs(residuals)))

# Adjust using Benjamini-Hochberg (FDR)
pvals_adj <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                    nrow = nrow(pvals),
                    dimnames = dimnames(pvals))

# Convert to long format
res_df <- as.data.frame(as.table(residuals)) %>%
  mutate(p_value = as.vector(pvals_adj),
         stars = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ ""))

colnames(res_df) <- c("status","subtype","residual","adj_pval","stars")
res_df$status <- ifelse(res_df$status == "naive", "Naive","Treated")
res_df$status <- factor(res_df$status, levels=c("Treated","Naive"))

status_plot <- ggplot(res_df, aes(x = subtype, y = status)) +
  geom_point(aes(size = abs(residual), fill = residual), shape = 21, color = "black", stroke = 0.6) +
  geom_text(aes(label = stars), vjust = -1.5, size = 4) +  # Stars above dots
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(title = "",
       size = "|Residual|", fill = "Residual",
       x = "Subtype",
       y = "Treatment Status")

status_plot

################################################################################
# Tarla status and subtype association
################################################################################

contin_table <- table(ctcs$tarla, ctcs$subtype) %>% 
  as.data.frame.matrix()

contin_table <- contin_table[c(1,2),]

test_res <- chisq.test(contin_table)

residuals <- test_res$stdres

# Calculate two-sided p-values from standard normal distribution
pvals <- 2 * (1 - pnorm(abs(residuals)))

# Adjust using Benjamini-Hochberg (FDR)
pvals_adj <- matrix(p.adjust(as.vector(pvals), method = "BH"),
                    nrow = nrow(pvals),
                    dimnames = dimnames(pvals))

# Convert to long format
res_df <- as.data.frame(as.table(residuals)) %>%
  mutate(p_value = as.vector(pvals_adj),
         stars = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ ""))

colnames(res_df) <- c("status","subtype","residual","adj_pval","stars")
res_df$status <- ifelse(res_df$status == "T", "Post-Tarlatamab","Pre-Tarlatamab")

tarla_plot <- ggplot(res_df, aes(x = subtype, y = status)) +
  geom_point(aes(size = abs(residual), fill = residual), shape = 21, color = "black", stroke = 0.6) +
  geom_text(aes(label = stars), vjust = -1.5, size = 4) +  # Stars above dots
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(title = "",
       size = "|Residual|", fill = "Residual",
       x = "Subtype",
       y = "Treatment Status")

tarla_plot

################################################################################
# SAVE FIGURES
################################################################################

jpeg("figures/all_samples_treatment_subtype_chi_results.jpg", width=180,height=100, units = "mm", res=1000)
print(status_plot)
dev.off()

jpeg("figures/all_samples_tarla_subtype_chi_results.jpg", width=180,height=100, units = "mm", res=1000)
print(tarla_plot)
dev.off()
