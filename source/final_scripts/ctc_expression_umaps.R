
# Save data with cluster assignments
sce <- readRDS("data/cytof_objects/sclc_all_samples_with_clusters.rds")

markers_to_use <- rowData(sce)$marker_name[rowData(sce)$marker_class == "state"]
y <- assay(sce, "exprs")
y <- t(y[markers_to_use,])

y <- as.data.frame(apply(y, 2, function(x) (x - min(x)) / (max(x) - min(x))))


xy <- reducedDim(sce, "UMAP")[, 1:2]
colnames(xy) <- c("x", "y")

df <- data.frame(colData(sce), xy,y, check.names = FALSE)

dim(y)

markers <- colnames(y)

# curr_marker <- markers[6]

for(curr_marker in markers){
   cat(curr_marker,"\n")
  
  # Plot UMAP
  curr_plot <- ggplot(df)+
    geom_point(aes(x=x, y=y, color=!!sym(curr_marker)),size=.1)+
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    labs(color = curr_marker)+
    scale_color_gradientn(colors = c("skyblue", "white", "red"))+
    theme_classic() +
    theme(panel.grid.minor = element_blank(), 
           strip.text = element_text(face = "bold", size=8), 
          axis.text = element_text(color = "black", size=8),
          axis.title = element_text(size=8),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8))
  
  
  
  jpeg(glue("figures/marker_expression_umaps/all_samples_{curr_marker}.jpg"), width=120,height=100, units = "mm", res=1000)
  print(curr_plot)
  dev.off()
}

