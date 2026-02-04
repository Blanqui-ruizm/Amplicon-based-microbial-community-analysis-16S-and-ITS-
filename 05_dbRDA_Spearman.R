# Title: db-RDA of 16S or ITS community composition using soil PCs
# Description: Evaluates the effects of soil properties on microbial community
#              structure. Multicollinearity among soil variables is assessed
#              using VIF, reduced via PCA, and tested using db-RDA. Permutations
#              are blocked by sampling time. Environmental vectors are fitted
#              for visualization only.
# Author: Blanca Ruiz

## ---------------------------------------------------------
## Load packages
## ---------------------------------------------------------

library(phyloseq)
library(vegan)
library(permute)
library(ggplot2)
library(ggnewscale)
library(grid)        # for unit()
library(openxlsx)

## ---------------------------------------------------------
## Ensure output directories exist
## ---------------------------------------------------------

if (!dir.exists("results"))         dir.create("results")
if (!dir.exists("results/figures")) dir.create("results/figures")
if (!dir.exists("results/tables"))  dir.create("results/tables")

## ---------------------------------------------------------
## Load input data
## ---------------------------------------------------------

distance_matrix <- readRDS("data/distance_matrix_hellinger.rds")
soil_16S_clean  <- readRDS("data/soil_16S_clean.rds")

meta <- as.data.frame(sample_data(soil_16S_clean))

## ---------------------------------------------------------
## Define soil variables
## ---------------------------------------------------------

soil_vars <- c("pH",
               "Available_Water_Percent",
               "Available_Phosphorus",
               "C.N",
               "Organic_matter",
               "Nitrate")

# Sanity checks
stopifnot(all(soil_vars %in% colnames(meta)))
stopifnot("Sampling.time" %in% colnames(meta))

env <- meta[, soil_vars]

## ---------------------------------------------------------
## Check multicollinearity (VIF)
## ---------------------------------------------------------

env_z <- scale(env)

rda_env <- rda(env_z)
vif_vals <- vif.cca(rda_env)

vif_table <- data.frame(
  Variable = names(vif_vals),
  VIF      = round(vif_vals, 2)
)

vif_table

write.xlsx(
  vif_table,
  file = "results/tables/Supplementary_Table_S0_Soil_VIF.xlsx",
  rowNames = FALSE
)

## ---------------------------------------------------------
## PCA on soil variables
## ---------------------------------------------------------

pca_env <- rda(env_z)

PCs <- as.data.frame(scores(pca_env,
                            display = "sites",
                            choices = 1:2))
colnames(PCs) <- c("Soil_PC1", "Soil_PC2")

dat <- cbind(meta, PCs)

## ---------------------------------------------------------
## dbRDA: community ~ soil PCs
## ---------------------------------------------------------

mod_pc <- capscale(
  distance_matrix ~ Soil_PC1 + Soil_PC2,
  data = dat,
  add  = TRUE
)

## ---------------------------------------------------------
## Permutation design (blocked by sampling time)
## ---------------------------------------------------------

ctrl <- how(nperm = 9999)
setBlocks(ctrl) <- dat$Sampling.time

anova(mod_pc, permutations = ctrl)
anova(mod_pc, by = "margin", permutations = ctrl)

## ---------------------------------------------------------
## dbRDA ordination plot with envfit overlay
## ---------------------------------------------------------

# Percentage of constrained variance
eig <- mod_pc$CCA$eig
pct <- 100 * eig / sum(eig)

x_label <- paste0("CAP1 (", round(pct[1], 1), "% constrained)")
y_label <- paste0("CAP2 (", round(pct[2], 1), "% constrained)")

# Site scores
sites_df <- as.data.frame(scores(mod_pc,
                                 display = "sites",
                                 choices = 1:2))
colnames(sites_df) <- c("CAP1", "CAP2")

# Fit original soil variables (visual interpretation only)
fit <- envfit(mod_pc,
              dat[, soil_vars, drop = FALSE],
              permutations = 9999)

env_vec <- as.data.frame(scores(fit, display = "vectors"))
colnames(env_vec)[1:2] <- c("CAP1", "CAP2")
env_vec$pval <- fit$vectors$pvals[rownames(env_vec)]

env_vec$PCgroup <- ifelse(abs(env_vec$CAP1) >= abs(env_vec$CAP2),
                          "PC1", "PC2")

arrow_scale <- max(abs(sites_df[, c("CAP1", "CAP2")])) * 0.8

env_coords <- env_vec
env_coords[, c("CAP1", "CAP2")] <-
  env_coords[, c("CAP1", "CAP2")] * arrow_scale

env_coords_adj <- env_coords
env_coords_adj[, c("CAP1", "CAP2")] <-
  env_coords_adj[, c("CAP1", "CAP2")] * 1.05

abbr <- c(
  "Available_Water_Percent" = "AWP",
  "Available_Phosphorus"    = "P",
  "Organic_matter"          = "OM",
  "Nitrate"                 = "NO3",
  "C.N"                     = "C:N",
  "pH"                      = "pH"
)

rownames(env_coords)     <- abbr[rownames(env_coords)]
rownames(env_coords_adj) <- abbr[rownames(env_coords_adj)]

p <- ggplot(sites_df, aes(CAP1, CAP2)) +
  geom_point(size = 3.3) +
  labs(x = x_label, y = y_label) +
  new_scale_color() +
  geom_segment(data = env_coords,
               aes(x = 0, y = 0,
                   xend = CAP1, yend = CAP2,
                   color = PCgroup),
               arrow = arrow(length = unit(0.3, "cm")),
               linewidth = 0.9) +
  geom_text(data = env_coords_adj,
            aes(label = rownames(env_coords_adj),
                color = PCgroup),
            size = 3.2,
            fontface = "bold",
            vjust = -0.5,
            hjust = -0.1) +
  scale_color_manual(values = c("PC1" = "#53868B",
                                "PC2" = "#8B475D")) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title   = element_text(face = "bold")
  )

p

ggsave(
  filename = "results/figures/Fig2c_dbRDA_soilPCs.pdf",
  plot     = p,
  width    = 10,
  height   = 6.5,
  dpi      = 300
)

## ---------------------------------------------------------
## Supplementary tables
## ---------------------------------------------------------

eig_pca <- eigenvals(pca_env)
var_exp <- eig_pca / sum(eig_pca) * 100

var_table <- data.frame(
  PC = c("PC1", "PC2"),
  Variance_explained_percent = round(var_exp[1:2], 1)
)

loadings <- as.data.frame(scores(pca_env,
                                 display = "species",
                                 choices = 1:2))
loadings$Variable <- rownames(loadings)
loadings <- loadings[, c("Variable", "PC1", "PC2")]

anova_global <- as.data.frame(anova(mod_pc,
                                    permutations = ctrl))
anova_global$Term <- rownames(anova_global)
rownames(anova_global) <- NULL

anova_margin <- as.data.frame(anova(mod_pc,
                                    by = "margin",
                                    permutations = ctrl))
anova_margin$Predictor <- rownames(anova_margin)
rownames(anova_margin) <- NULL

write.xlsx(
  list(
    "PCA_loadings"           = loadings,
    "PCA_variance_explained" = var_table
  ),
  file = "results/tables/Supplementary_Table_S1_PCA_soil.xlsx",
  rowNames = FALSE
)

write.xlsx(
  list(
    "dbRDA_global_model"   = anova_global,
    "dbRDA_marginal_terms" = anova_margin
  ),
  file = "results/tables/Supplementary_Table_S2_dbRDA_PC1_PC2.xlsx",
  rowNames = FALSE
)



