
# Title: db-RDA and Spearman correlations (Family level)
# Description: Aggregates data at the Family level, applies Hellinger 
#              transformation, performs Mantel test, db-RDA with environmental 
#              variables, and Spearman correlations. Generates heatmap and biplot.
# Author: Blanca Ruiz

# Load packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(tidyr)
library(dplyr)
library(tibble)
library(openxlsx)
library(reshape2)
library(RColorBrewer)
library(writexl)
library(grid)   # for unit() used in arrow()

# Ensure output folder exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Load input data (change if needed)
soil_16S_hellinger <- readRDS("data/soil_16S_hellinger.rds")
distance_matrix    <- readRDS("data/distance_matrix_hellinger.rds")
soil_16S_clean     <- readRDS("data/soil_16S_clean.rds")


# 1) Aggregate abundances at Family level

soil_16S_family <- tax_glom(soil_16S_clean, taxrank = "Family")
family_abundance <- psmelt(soil_16S_family)

# Keep only required columns
family_abundance <- family_abundance %>%
  select(Sample, Abundance, Replicate, Sampling.time, Management, 
         pH, Available_Water_Percent, Available_Phosphorus, C.N, 
         Organic_matter, Nitrate, Family)

# Top 25 families + group others as "Other"
top_families <- family_abundance %>%
  group_by(Family) %>%
  summarize(TotalAbundance = sum(Abundance)) %>%
  top_n(25, TotalAbundance) %>%
  pull(Family)

family_abundance_filtered <- family_abundance %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%
  group_by(Sample, Family) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# Convert to matrix
family_matrix <- family_abundance_filtered %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = "Sample")

# Remove empty samples (all zeros) to avoid NaN in transforms
family_matrix <- family_matrix[rowSums(family_matrix) > 0, ]

# Apply Hellinger transformation (sqrt of relative abundance)
family_matrix_hellinger <- decostand(family_matrix, method = "hellinger")

# 2) Environmental data

environmental_data <- sample_data(soil_16S_clean)[
  rownames(family_matrix_hellinger),  # align samples
  c("pH", "Available_Water_Percent", "Available_Phosphorus", "C.N", "Organic_matter", "Nitrate")
]
environmental_data <- as.data.frame(environmental_data)

# Scale environmental variables (unit variance, zero mean)
environmental_data_scaled <- as.data.frame(scale(environmental_data))


# 3) Mantel test (Bray vs Euclidean)

env_dist_matrix    <- dist(environmental_data_scaled, method = "euclidean")
family_dist_matrix <- vegdist(family_matrix_hellinger, method = "bray")

mantel_result_families <- mantel(
  family_dist_matrix, env_dist_matrix,
  method = "spearman", permutations = 9999
)
print(mantel_result_families)


# 4) db-RDA (capscale) with environmental variables

db_rda_families <- capscale(
  family_matrix_hellinger ~ pH + Available_Water_Percent + 
    Available_Phosphorus + C.N + Organic_matter + Nitrate,
  data = environmental_data_scaled, distance = "bray"
)

# Permutational ANOVA: global and by terms (with FDR)
anova_global_families <- anova(db_rda_families, permutations = 9999)
anova_terms_families  <- anova(db_rda_families, by = "terms", permutations = 9999)

p_values     <- anova_terms_families$`Pr(>F)`
p_values_fdr <- p.adjust(p_values, method = "fdr")
anova_terms_families_fdr <- cbind(anova_terms_families, p_value_fdr = p_values_fdr)

# Save ANOVA-by-terms results
write_xlsx(anova_terms_families_fdr, "results/dbRDA_anova_terms_families.xlsx")


# 5) Spearman correlations (Families vs Environmental variables) + FDR

colnames(environmental_data_scaled) <- c("pH", "AWP", "P", "C:N", "OM", "N")

cor_results <- data.frame()
for (family in colnames(family_matrix_hellinger)) {
  for (variable in colnames(environmental_data_scaled)) {
    cor_test <- cor.test(
      family_matrix_hellinger[, family],
      environmental_data_scaled[, variable],
      method = "spearman"
    )
    cor_results <- rbind(cor_results, data.frame(
      Family  = family,
      Variable = variable,
      Rho      = cor_test$estimate,
      P_value  = cor_test$p.value
    ))
  }
}
cor_results$P_value_FDR <- p.adjust(cor_results$P_value, method = "fdr")

# Keep only significant pairs (FDR < 0.05)
cor_filtered <- cor_results %>%
  filter(P_value_FDR < 0.05) %>%
  mutate(Significance = case_when(
    P_value_FDR < 0.001 ~ "***",
    P_value_FDR < 0.01  ~ "**",
    P_value_FDR < 0.05  ~ "*",
    TRUE ~ ""
  ))

# Save correlation table
write.xlsx(cor_filtered, file = "results/spearman_correlations.xlsx", rowNames = FALSE)


# 6) Heatmap of significant correlations

if (nrow(cor_filtered) > 0) {
  cor_matrix <- cor_filtered %>%
    select(Family, Variable, Rho) %>%
    pivot_wider(names_from = Variable, values_from = Rho, values_fill = 0) %>%
    column_to_rownames("Family")
  
  family_dendro <- hclust(dist(cor_matrix), method = "average")
  ordered_families <- family_dendro$labels[family_dendro$order]
  cor_filtered$Family <- factor(cor_filtered$Family, levels = ordered_families)
  
  p <- ggplot(cor_filtered, aes(x = Variable, y = Family, fill = Rho)) +
    geom_tile(width = 1, height = 1) +
    scale_fill_gradient2(low = "#27408B", mid = "white", high = "#8B1A1A", midpoint = 0) +
    geom_text(aes(label = Significance), color = "black", size = 5, fontface = "bold") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, face = "italic"),
      legend.position = "bottom"
    )
  ggsave("results/heatmap_spearman_families.pdf", plot = p, width = 6, height = 7)
}


# 7) db-RDA biplot

species_coords <- scores(db_rda_families, display = "species", scaling = 2)
env_coords     <- scores(db_rda_families, display = "bp",      scaling = 2)
sites_df       <- as.data.frame(scores(db_rda_families, display = "sites", scaling = 2))

# Append minimal metadata for aesthetics
sites_df <- cbind(sites_df, as.data.frame(sample_data(soil_16S_clean)[rownames(sites_df), ]))

# Optional relabel for vectors in the plot
rownames(env_coords) <- c("pH", "AWP", "P", "C:N", "OM", "N")

p <- ggplot(sites_df, aes(x = CAP1, y = CAP2, color = Management, shape = Sampling.time)) +
  geom_point(size = 4) +
  geom_segment(data = as.data.frame(species_coords), inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#EEAEEE", linewidth = 0.8) +
  geom_segment(data = as.data.frame(env_coords), inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#53868B", linewidth = 0.8) +
  theme_minimal(base_size = 12) +
  labs(x = "CAP1", y = "CAP2") +
  scale_color_manual(values = c("Organic" = "#9ACD32", "Conventional" = "#EE9A00")) +
  scale_shape_manual(values = c("February" = 16, "May" = 15, "July" = 17, "October" = 18))

ggsave("results/dbRDA_biplot_families.pdf", plot = p, width = 8, height = 6)

