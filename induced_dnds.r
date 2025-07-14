library(dplyr)
library(tidyr)

stage_cols <- c("Somatic", "Pre.meiotic", "Meiotic", "Post.meiotic")

induced_list <- list()

for (stage in stage_cols) {
  other_stages <- setdiff(stage_cols, stage)
  
  induced_genes <- counts1_cpm %>%
    filter(!is.na(dnds)) %>%
    rowwise() %>%
    mutate(
      median_other = median(c_across(all_of(other_stages)), na.rm = TRUE),
      induced = .data[[stage]] > 2 * median_other
    ) %>%
    ungroup() %>%
    filter(induced) %>%
    select(ensembl_gene_id, dnds, Chr) %>%
    mutate(Stage = stage)
  
  induced_list[[stage]] <- induced_genes
}

# Combine all stages
induced_df <- bind_rows(induced_list)
induced_df$Stage <- factor(induced_df$Stage, levels = stage_cols)
induced_df$Chr <- factor(induced_df$Chr)


library(ggplot2)

ggplot(induced_df, aes(x = Stage, y = dnds, fill = Chr)) +
  geom_boxplot(outlier.size = 0.8, position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(title = "dN/dS of Induced Genes (2x Median Threshold)",
       x = "Stage",
       y = "dN/dS") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")
