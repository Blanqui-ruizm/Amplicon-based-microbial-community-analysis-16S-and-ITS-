
# Title: Beta diversity analysis
# Description: Performs Hellinger transformation, computes Bray–Curtis distances,
#              checks homogeneity of dispersion (PERMDISP), 
#              runs PERMANOVA (management, sampling time, interaction),
#              and generates PCoA plots with ellipses.
# Author: Blanca Ruiz


# ------------------------------------------------------------------------------
# Title: Beta diversity analysis
# Description: Performs Hellinger transformation, computes Bray–Curtis distances,
#              checks homogeneity of dispersion (PERMDISP),
#              runs PERMANOVA (management, sampling time, interaction),
#              and generates PCoA plots with ellipses.
# Author: Blanca Ruiz
# ------------------------------------------------------------------------------

library(phyloseq)
library(vegan)
library(ggplot2)
library(writexl)

# Ensure output folder exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# 1) Load phyloseq object
soil_ps <- readRDS("data/soil_16S_clean.rds")   # change filename if needed

# 2) Hellinger transform (rows = samples) + Bray–Curtis
M <- as(otu_table(soil_ps), "matrix")
if (taxa_are_rows(soil_ps)) M <- t(M)                 # decostand needs rows = samples
M <- M[rowSums(M) > 0, , drop = FALSE]                # drop empty samples (avoids NaN)
M_hell <- decostand(M, method = "hellinger")
bray   <- vegdist(M_hell, method = "bray")

# Rebuild a phyloseq with transformed abundances
soil_ps_hellinger <- phyloseq(
  otu_table(M_hell, taxa_are_rows = FALSE),
  sample_data(sample_data(soil_ps)[rownames(M_hell), , drop = FALSE]),
  tax_table(soil_ps)
)

# Save intermediates for downstream db-RDA
saveRDS(soil_ps_hellinger, file = "data/soil_16S_hellinger.rds")
saveRDS(bray,                file = "data/distance_matrix_hellinger.rds")

# 3) Homogeneity of dispersion (PERMDISP)
sdf <- as(sample_data(soil_ps_hellinger), "data.frame")
sdf <- sdf[labels(bray), , drop = FALSE]              # align order to distance labels

disp_mgmt <- betadisper(bray, sdf$Management)
disp_time <- betadisper(bray, sdf$Sampling.time)

permdisp_mgmt <- permutest(disp_mgmt, permutations = 999)
permdisp_time <- permutest(disp_time, permutations = 999)

permdisp_results <- data.frame(
  Factor  = c("Management", "Sampling Time"),
  F_value = c(permdisp_mgmt$tab[1, "F"], permdisp_time$tab[1, "F"]),
  p_value = c(permdisp_mgmt$tab[1, "Pr(>F)"], permdisp_time$tab[1, "Pr(>F)"])
)
permdisp_results$p_adj_fdr <- p.adjust(permdisp_results$p_value, method = "fdr")
write_xlsx(permdisp_results, "results/PERMDISP_results.xlsx")

# 4) PERMANOVA (marginal terms + interaction)
perm_terms <- adonis2(bray ~ Management + Sampling.time, data = sdf, permutations = 999, by = "term")
perm_int   <- adonis2(bray ~ Management * Sampling.time, data = sdf, permutations = 999)

permanova_results_df <- data.frame(
  Term        = rownames(perm_terms)[1:2],
  R2          = perm_terms$R2[1:2],
  F_value     = perm_terms$F[1:2],
  p_value     = perm_terms$`Pr(>F)`[1:2],
  p_value_fdr = p.adjust(perm_terms$`Pr(>F)`[1:2], method = "fdr")
)

permanova_interaction_df <- data.frame(
  Term        = rownames(perm_int),
  R2          = perm_int$R2,
  F_value     = perm_int$F,
  p_value     = perm_int$`Pr(>F)`,
  p_value_fdr = p.adjust(perm_int$`Pr(>F)`, method = "fdr")
)

write_xlsx(
  list(
    "PERMANOVA_Results"     = permanova_results_df,
    "PERMANOVA_Interaction" = permanova_interaction_df
  ),
  "results/PERMANOVA_results.xlsx"
)

# 5) PCoA plots
pcoa <- cmdscale(bray, k = 2, eig = TRUE)
coords <- as.data.frame(pcoa$points)
colnames(coords) <- c("Axis.1", "Axis.2")
coords <- cbind(coords, sdf)

# variance explained (for labels)
var_exp <- pcoa$eig / sum(pcoa$eig)
xlab <- paste0("PCoA Axis 1 (", round(var_exp[1] * 100, 2), "%)")
ylab <- paste0("PCoA Axis 2 (", round(var_exp[2] * 100, 2), "%)")

group_colors <- c("Organic" = "#9ACD32", "Conventional" = "#EE9A00")
time_colors  <- c("February" = "#8EE5EE", "May" = "#43CD80",
                  "July" = "#FFD700", "October" = "#DA70D6")
shape_values <- c("February" = 16, "May" = 15, "July" = 17, "October" = 18)

# By Management
p_mgmt <- ggplot(coords, aes(x = Axis.1, y = Axis.2, color = Management)) +
  geom_point(aes(shape = Sampling.time), size = 3) +
  stat_ellipse(aes(fill = Management), geom = "polygon", alpha = 0.3) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors, guide = guide_legend(override.aes = list(alpha = 0.3))) +
  labs(x = xlab, y = ylab) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))
ggsave("results/PCoA_BrayCurtis_Management.pdf", p_mgmt, width = 8, height = 6)

# By Sampling Time
p_time <- ggplot(coords, aes(x = Axis.1, y = Axis.2, color = Sampling.time)) +
  geom_point(aes(shape = Sampling.time), size = 3) +
  stat_ellipse(aes(fill = Sampling.time), geom = "polygon", alpha = 0.3) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = time_colors) +
  scale_fill_manual(values = time_colors, guide = guide_legend(override.aes = list(alpha = 0.3))) +
  labs(x = xlab, y = ylab) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.line = element_line(color = "black"))
ggsave("results/PCoA_BrayCurtis_SamplingTime.pdf", p_time, width = 8, height = 6)