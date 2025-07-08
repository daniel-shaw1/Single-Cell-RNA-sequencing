# Set induction threshold
threshold <- 1.5

cell_types <- c("Pre.meiotic", "Meiotic", "Post.meiotic", "Somatic")

# Add category labels first (make sure your gene name columns match exactly)
combined_data$Category <- "Other"
combined_data$Category[combined_data$GeneName %in% Xgametolog$Genename] <- "Xgametolog"
combined_data$Category[combined_data$GeneName %in% Xonly$Genename] <- "Xonly"
combined_data$Category[combined_data$GeneName %in% Yonly$Genename] <- "Yonly"

induced_genes <- list()
dnds_values <- list()

for (ct in cell_types) {
  other_cts <- setdiff(cell_types, ct)
  medians <- apply(combined_data[, other_cts], 1, median, na.rm = TRUE)
  
  induced <- combined_data[, ct] > threshold * medians & !is.na(combined_data$dnds)
  
  induced_genes[[ct]] <- combined_data[induced, ]
  dnds_values[[ct]] <- combined_data$dnds[induced]
}

# Build the grouped list for plotting
dnds_grouped <- list()

for (stage in names(dnds_values)) {
  genes <- induced_genes[[stage]]
  categories <- genes$Category
  dnds <- genes$dnds
  
  for (cat in unique(categories)) {
    group_label <- paste(cat, stage, sep = "_")
    dnds_grouped[[group_label]] <- dnds[categories == cat]
  }
}

# Now create a combined data frame for ggplot2 plotting
library(ggplot2)

plot_df <- do.call(rbind, lapply(names(dnds_grouped), function(name) {
  data.frame(
    dnds = dnds_grouped[[name]],
    group = name,
    stringsAsFactors = FALSE
  )
}))

plot_df$Category <- sub("^(.*?)_.*$", "\\1", plot_df$group)
plot_df$Stage <- sub("^.*?_(.*)$", "\\1", plot_df$group)

plot_df$Category <- factor(plot_df$Category, levels = c("Xgametolog", "Xonly", "Yonly", "Other"))
plot_df$Stage <- factor(plot_df$Stage, levels = cell_types)

ggplot(plot_df, aes(x = Stage, y = dnds, fill = Category)) +
  geom_boxplot(outlier.size = 0.8) +
  theme_bw() +
  labs(title = "dN/dS by Gene Category and Induction Stage (Threshold = 1.5x Median)",
       x = "Induction Stage",
       y = "dN/dS") +
  scale_fill_manual(values = c("Xgametolog" = "#1b9e77", "Xonly" = "#d95f02", "Yonly" = "#7570b3", "Other" = "grey70")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
