library(ggplot2)
library (dplyr)
library(stats)

mps1i <- c("Reversine", "AZ3146", "MPI.0479605", "CFI.402257")

gene_ttest_rank <- read.csv("treatment_vs_gene_corr_ttest_ranking_neg_effect_only_PRISM_23.csv")
gene_ttest_rank$treatment <- gsub('\\.\\.BRD\\.BRD.*','',gene_ttest_rank$treatment)

mps1i_data <- dplyr::filter(gene_ttest_rank, tolower(treatment) %in% tolower(mps1i))
gene_ttest_rank_filter <- gene_ttest_rank[gene_ttest_rank$rank != -1,]

avg_CDC20_rank <- round(mean(gene_ttest_rank_filter$precent))
avg_CDC20_rank_mps1i <- round(mean(mps1i_data$precent))


density_data <- density(gene_ttest_rank_filter$precent)


density_df <- data.frame(
  x = density_data$x,
  y = density_data$y
)

mps1i_density <- density_df %>%
  dplyr::filter(x %in% mps1i_data$precent)


avg_CDC20_rank_mps1i_density <- density_df %>%
  filter(x == avg_CDC20_rank_mps1i)
avg_CDC20_rank_density <- density_df %>%
  filter(x == avg_CDC20_rank)

density_data <- density(gene_ttest_rank_filter$precent)

density_df <- data.frame(
  x = density_data$x,
  y = density_data$y
)

plot1 <- ggplot(data = gene_ttest_rank_filter, aes(x = precent)) +
  geom_density() +
  geom_point(data = mps1i_data, aes(x = precent, y = approx(density_df$x, density_df$y, precent)$y), color = "red") + 
  geom_vline(xintercept = avg_CDC20_rank_mps1i, color = 'red', linetype = "dashed") +
  geom_vline(xintercept = avg_CDC20_rank, linetype = "dashed") +
  annotate("text", x = avg_CDC20_rank_mps1i - 3, y = approx(density_df$x, density_df$y, avg_CDC20_rank_mps1i-0.2)$y, label = "mps1i", angle = 90, color = "red", size = 3) +
  annotate("text", x = avg_CDC20_rank - 2, y = approx(density_df$x, density_df$y, avg_CDC20_rank)$y + 0.0005, label = "all drugs", angle = 90, size = 3) +
  labs(x = "Percentile of ranking (top to bottom)", y = "Density", title = "CDC20 gene ranking distribution by p-value", subtitle ="(unadjusted)")

ggsave("CDC20_gene_ranking_density_plot.png", plot1, width = 768 * 2, height = 498 * 2, units = "px",
       dpi = 300)


