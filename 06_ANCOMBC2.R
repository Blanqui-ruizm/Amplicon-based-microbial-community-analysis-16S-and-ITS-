
# Title: ANCOM-BC2 at Family level + combined figure (LFC + boxplots)
# Description: Runs ANCOM-BC2 at Family level with filters, exports the table,
#              and builds a two-panel figure: log2 fold-change with CI and
#              relative-abundance boxplots by group.
# Author: Blanca Ruiz


library(phyloseq)
library(mia)
library(ANCOMBC)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr)

set.seed(353)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# Input
ps_file <- "data/soil_16S_clean.rds"
ps_base <- readRDS(ps_file)
tax_table(ps_base)[is.na(tax_table(ps_base))] <- "Unclassified"

# Config
group_var <- "Management"      # "Conventional"/"Organic"
threshold_count   <- 10        # min counts per sample
threshold_samples <- 5         # min samples meeting that count
alpha_q <- 0.05
em_fast <- list(tol = 1e-3, max_iter = 25)

# Reference: Conventional
sample_data(ps_base)[[group_var]] <- relevel(
  factor(sample_data(ps_base)[[group_var]]),
  ref = "Conventional"
)

# Helpers
as_taxa_rows <- function(ps){
  M <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) M <- t(M); M
}

filter_by_count_prevalence <- function(ps, cmin, smin){
  M <- as_taxa_rows(ps)
  keep <- apply(M, 1, function(x) sum(x >= cmin) >= smin)
  prune_taxa(rownames(M)[keep], ps)
}

prune_zero_variance_all <- function(ps, group_var){
  M <- as_taxa_rows(ps); grp <- factor(sample_data(ps)[[group_var]])
  v_all <- apply(M, 1, var)
  keep_total <- is.finite(v_all) & v_all > 0
  keep_bygrp <- sapply(1:nrow(M), function(i) all(tapply(M[i, ], grp, var) > 0, na.rm = TRUE))
  keep <- rownames(M)[keep_total & keep_bygrp]
  prune_taxa(keep, ps)
}

# 1) Collapse to Family and filter
ps_family <- tax_glom(ps_base, taxrank = "Family")
ps_family <- filter_by_count_prevalence(ps_family, threshold_count, threshold_samples)
ps_family <- prune_zero_variance_all(ps_family, group_var)
stopifnot(ntaxa(ps_family) > 0)

# 2) Taxonomy and names
tx <- as.data.frame(tax_table(ps_family))
nm <- as.character(tx[["Family"]]); nm[is.na(nm) | nm == ""] <- "Unclassified"
taxa_names(ps_family) <- make.unique(nm)

ord_tax <- c("Kingdom","Phylum","Class","Order","Family","Genus")
for (rname in ord_tax)
  if (!rname %in% names(tx)) tx[[rname]] <- NA_character_
tx$Taxon <- taxa_names(ps_family)
tx$Lineage <- apply(tx[, intersect(colnames(tx), c(ord_tax,"Species")), drop=FALSE], 1, function(z){
  z <- as.character(z); z <- z[z != "" & !is.na(z)]
  if (!length(z)) "Unclassified" else paste(z, collapse = "; ")
})

# 3) ANCOM-BC2
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps_family)
fit <- ancombc2(
  data = tse, assay_name = "counts",
  fix_formula = group_var, group = group_var,
  p_adj_method = "fdr",
  prv_cut = 0, lib_cut = 0,
  struc_zero = FALSE, neg_lb = FALSE,
  alpha = alpha_q, global = FALSE,
  em_control = em_fast
)
tab <- as.data.frame(fit$res)
taxcol <- if ("taxon" %in% names(tab)) "taxon" else colnames(tab)[1]

# 4) Output table (Family)
out_family <- data.frame(Taxon = tab[[taxcol]], stringsAsFactors = FALSE)
out_family <- merge(out_family, tx[, c("Taxon", ord_tax, "Lineage")],
                    by = "Taxon", all.x = TRUE, sort = FALSE)
out_family$Rank <- "Family"
out_family <- cbind(out_family, tab[ , setdiff(names(tab), taxcol), drop = FALSE])
out_family <- out_family[, c("Rank", ord_tax, "Taxon", "Lineage",
                             setdiff(names(out_family), c("Rank", ord_tax, "Taxon","Lineage"))),
                         drop = FALSE]

write.xlsx(out_family, "results/ANCOMBC2_Family_Organic_vs_Conventional.xlsx", overwrite = TRUE)
cat("Exported: results/ANCOMBC2_Family_Organic_vs_Conventional.xlsx\n")

# 5) Figure (LFC + CI) + Relative abundance boxplots (log10)
# Significant families (FDR < 0.05 for Management contrast)
sig <- out_family %>%
  filter(q_ManagementOrganic < 0.05) %>%
  rename(
    LFC = lfc_ManagementOrganic,
    SE  = se_ManagementOrganic,
    q   = q_ManagementOrganic,
    p   = p_ManagementOrganic
  ) %>%
  mutate(
    association = ifelse(LFC > 0, "Organic", "Conventional"),
    CleanName   = Taxon
  )

if (nrow(sig) > 0) {
  # Order: positives on top, negatives below
  pos <- sig %>% filter(LFC > 0) %>% arrange(desc(LFC)) %>% pull(CleanName)
  neg <- sig %>% filter(LFC < 0) %>% arrange(LFC)        %>% pull(CleanName)
  order_taxa <- rev(c(pos, neg))
  sig$CleanName <- make.unique(as.character(sig$CleanName))
  sig$CleanName <- factor(sig$CleanName, levels = order_taxa)
  
  # Relative abundances from ps_family
  ps_rel <- transform_sample_counts(ps_family, function(x) x / sum(x))
  df_family2 <- psmelt(ps_rel) %>%
    mutate(
      FamilyClean2 = as.character(OTU),  # collapsed family name
      Management   = factor(trimws(as.character(Management)),
                            levels = c("Conventional","Organic"))
    )
  
  df_family_sig <- df_family2 %>%
    filter(FamilyClean2 %in% levels(sig$CleanName)) %>%
    mutate(FamilyClean2 = factor(FamilyClean2, levels = levels(sig$CleanName)))
  
  # Left panel: LFC with 95% CI
  sig_plot <- sig %>% mutate(ci_lower = LFC - 1.96 * SE, ci_upper = LFC + 1.96 * SE)
  p_bar_lfc <- ggplot(sig_plot, aes(x = LFC, y = CleanName)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.2, color = "gray20", linewidth = 1.1) +
    geom_point(aes(color = association), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Organic" = "#9ACD32", "Conventional" = "#EE9A00")) +
    labs(x = "Log2 Fold Change", y = NULL, color = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 12),
      panel.grid.major.y = element_line(color = "gray85")
    )
  
  # Right panel: relative abundance (log10)
  p_box_horiz <- ggplot(df_family_sig, aes(y = FamilyClean2, x = Abundance, fill = Management)) +
    geom_boxplot(outlier.size = 1) +
    scale_fill_manual(values = c("Conventional" = "#EE9A00", "Organic" = "#9ACD32")) +
    scale_x_log10() +
    labs(x = "Relative Abundance (log10)", y = NULL, fill = "Group") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
  
  final_fig <- ggpubr::ggarrange(
    p_bar_lfc, p_box_horiz,
    ncol = 2, align = "h",
    widths = c(0.8, 1.2)
  )
  # Optional small left margin
  final_fig <- ggpubr::annotate_figure(final_fig, left = ggpubr::text_grob(" ", rot = 90, size = 14))
  
  ggsave("results/ANCOMBC2_family_combined_plot.pdf", final_fig, width = 14, height = 8)
  cat("Exported figure: results/ANCOMBC2_family_combined_plot.pdf\n")
} else {
  message("No significant families at FDR < 0.05; skipping figure.")
}
