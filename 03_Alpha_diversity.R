
# Title: Alpha diversity analysis
# Description: Calculates alpha diversity indices (Observed, Shannon, Simpson), generates boxplots, 
#              and performs statistical tests (Wilcoxon, Kruskal-Wallis, pairwise Wilcoxon).
# Author: Blanca Ruiz



library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(rstatix)
library(writexl)

# Ensure output folder exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# 1) Load phyloseq object
soil_ps <- readRDS("data/soil_16S_clean.rds")  # change if needed

# 2) Compute alpha diversity indices
alpha_div <- estimate_richness(
  soil_ps,
  measures = c("Observed", "Shannon", "Simpson")
) %>% as.data.frame()

# Add metadata (aligned to alpha_div row order for safety)
metadata <- as.data.frame(sample_data(soil_ps))
if (!all(rownames(alpha_div) %in% rownames(metadata))) {
  stop("Sample names in metadata and alpha diversity do not match.")
}
metadata <- metadata[rownames(alpha_div), , drop = FALSE]
alpha_div <- cbind(metadata, alpha_div)

# Add Inverse Simpson
alpha_div <- alpha_div %>% mutate(InvSimpson = 1 / Simpson)

# 3) Plots: Alpha diversity by management
alpha_long_mgmt <- alpha_div %>%
  select(Management, Observed, Shannon, InvSimpson) %>%
  pivot_longer(cols = c("Observed", "Shannon", "InvSimpson"),
               names_to = "Index", values_to = "Value")

alpha_long_mgmt$Index <- factor(alpha_long_mgmt$Index,
                                levels = c("Observed", "Shannon", "InvSimpson"))

fill_colors <- c("Organic" = "#9ACD32", "Conventional" = "#EE9A00")

p_mgmt <- ggplot(alpha_long_mgmt, aes(x = Management, y = Value, fill = Management)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.1, size = 2, color = "black", alpha = 0.7) +
  facet_wrap(~ Index, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = fill_colors) +
  labs(y = "Alpha diversity", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave("results/AlphaDiversity_byManagement.pdf", p_mgmt, width = 12, height = 8)

# 4) Plots: Alpha diversity by sampling time
time_colors <- c(
  "February" = "#8EE5EE",
  "May"      = "#43CD80",
  "July"     = "#FFD700",
  "October"  = "#DA70D6"
)

alpha_div$Sampling.time <- factor(alpha_div$Sampling.time,
                                  levels = c("February", "May", "July", "October"))

# Dynamic y-positions for significance asterisks
y_max <- max(alpha_div$Observed, na.rm = TRUE)
y_positions <- y_max * c(1.05, 1.10, 1.15, 1.20)

p_observed <- ggplot(alpha_div, aes(x = Sampling.time, y = Observed, fill = Sampling.time)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  scale_fill_manual(values = time_colors) +
  labs(y = "Observed richness") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_signif(
    comparisons = list(
      c("February", "July"),
      c("February", "October"),
      c("May", "July"),
      c("May", "October")
    ),
    annotations = c("**", "**", "**", "*"),
    y_position  = y_positions,
    tip_length  = 0.02,
    textsize    = 5
  )

p_shannon <- ggplot(alpha_div, aes(x = Sampling.time, y = Shannon, fill = Sampling.time)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  scale_fill_manual(values = time_colors) +
  labs(y = "Shannon index") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

p_simpson <- ggplot(alpha_div, aes(x = Sampling.time, y = Simpson, fill = Sampling.time)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +
  scale_fill_manual(values = time_colors) +
  labs(y = "Simpson index") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

final_plot <- p_observed / p_shannon / p_simpson + plot_layout(ncol = 1)
ggsave("results/AlphaDiversity_bySamplingTime.pdf", final_plot, width = 12, height = 9)

# 5) Statistics
# Wilcoxon tests (management per sampling time)
observed_stats <- alpha_div %>%
  group_by(Sampling.time) %>%
  wilcox_test(Observed ~ Management) %>%
  adjust_pvalue(method = "fdr")

shannon_stats <- alpha_div %>%
  group_by(Sampling.time) %>%
  wilcox_test(Shannon ~ Management) %>%
  adjust_pvalue(method = "fdr")

simpson_stats <- alpha_div %>%
  group_by(Sampling.time) %>%
  wilcox_test(Simpson ~ Management) %>%
  adjust_pvalue(method = "fdr")

combined_stats <- bind_rows(
  observed_stats %>% mutate(Index = "Observed"),
  shannon_stats  %>% mutate(Index = "Shannon"),
  simpson_stats  %>% mutate(Index = "Simpson")
)

write_xlsx(combined_stats, "results/AlphaDiversity_Wilcoxon.xlsx")

# Kruskal-Wallis tests (effect of sampling time)
kruskal_observed <- alpha_div %>% kruskal_test(Observed ~ Sampling.time)
kruskal_shannon  <- alpha_div %>% kruskal_test(Shannon  ~ Sampling.time)
kruskal_simpson  <- alpha_div %>% kruskal_test(Simpson  ~ Sampling.time)

# Pairwise Wilcoxon by sampling time
pairwise_observed <- alpha_div %>%
  pairwise_wilcox_test(Observed ~ Sampling.time, p.adjust.method = "fdr")
pairwise_shannon <- alpha_div %>%
  pairwise_wilcox_test(Shannon  ~ Sampling.time, p.adjust.method = "fdr")
pairwise_simpson <- alpha_div %>%
  pairwise_wilcox_test(Simpson  ~ Sampling.time, p.adjust.method = "fdr")

write_xlsx(
  list(
    "Observed" = pairwise_observed,
    "Shannon"  = pairwise_shannon,
    "Simpson"  = pairwise_simpson
  ),
  "results/AlphaDiversity_PairwiseWilcoxon.xlsx"
)

# 6) Export indices and summaries
sheet1 <- alpha_div %>%
  select(Sampling.time, Management, Observed, Shannon, InvSimpson)

sheet2 <- sheet1 %>%
  group_by(Sampling.time, Management) %>%
  summarise(across(c(Observed, Shannon, InvSimpson), list(mean = mean, sd = sd)), .groups = "drop")

write_xlsx(
  list(
    "Raw_Indices" = sheet1,
    "Summary_by_Time_and_Management" = sheet2
  ),
  "results/AlphaDiversity_summary.xlsx"
)