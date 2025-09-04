
# Title: Taxonomic summary tables (Phylum / Family)
# Description: Generates relative-abundance summary tables at specified
#              taxonomic ranks (Phylum, Family). Calculates means and SDs
#              by management and sampling time. Exports results to Excel.
# Author: Blanca Ruiz


library(phyloseq)
library(dplyr)
library(tidyr)
library(writexl)

# Ensure output folder exists
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# 1) Load phyloseq object (change filename if needed)
soil_ps <- readRDS("data/soil_16S_clean.rds")

# Transform to relative abundance (per-sample)
soil_rel <- transform_sample_counts(soil_ps, function(x) x / sum(x))

# 2) Helper to build summary tables for a taxonomic rank
summarize_taxa <- function(ps_obj, taxrank, output_file) {
  if (!taxrank %in% rank_names(ps_obj)) {
    warning("Rank '", taxrank, "' not found in tax_table. Skipping.")
    return(invisible(NULL))
  }
  
  rank_col <- rlang::sym(taxrank)
  
  # Aggregate counts at the given rank
  tax_table_long <- psmelt(tax_glom(ps_obj, taxrank = taxrank))
  tax_table_long[[taxrank]] <- ifelse(
    is.na(tax_table_long[[taxrank]]) | tax_table_long[[taxrank]] == "",
    "Unclassified",
    as.character(tax_table_long[[taxrank]])
  )
  
  # Summed abundances per sample (already relative per-sample)
  base_long <- tax_table_long %>%
    group_by(!!rank_col, Sample, Management, Sampling.time) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # General means by management
  general_summary <- base_long %>%
    group_by(!!rank_col, Management) %>%
    summarise(
      Mean = round(mean(Abundance) * 100, 2),
      SD   = round(sd(Abundance) * 100, 2),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = Management, values_from = c(Mean, SD), names_sep = "_")
  
  # Means by management and sampling time
  means_table <- base_long %>%
    group_by(!!rank_col, Management, Sampling.time) %>%
    summarise(Mean = round(mean(Abundance) * 100, 2), .groups = "drop") %>%
    mutate(label = paste0("Mean_", Sampling.time, "_", Management)) %>%
    select(!!rank_col, label, Mean) %>%
    pivot_wider(names_from = label, values_from = Mean)
  
  # SDs by management and sampling time
  sd_table <- base_long %>%
    group_by(!!rank_col, Management, Sampling.time) %>%
    summarise(SD = round(sd(Abundance) * 100, 2), .groups = "drop") %>%
    mutate(label = paste0("SD_", Sampling.time, "_", Management)) %>%
    select(!!rank_col, label, SD) %>%
    pivot_wider(names_from = label, values_from = SD)
  
  # Merge all parts
  final_table <- general_summary %>%
    left_join(means_table, by = taxrank) %>%
    left_join(sd_table,    by = taxrank)
  
  # Optional: order by Conventional mean if present
  if ("Mean_Conventional" %in% colnames(final_table)) {
    final_table <- final_table %>% arrange(desc(Mean_Conventional))
  }
  
  # Export to results/
  out_path <- file.path("results", output_file)
  write_xlsx(final_table, out_path)
  cat("Saved summary table:", out_path, "\n")
  
  invisible(final_table)
}

# 3) Run for Phylum and Family
summarize_taxa(soil_rel, "Phylum", "Phylum_summary.xlsx")
summarize_taxa(soil_rel, "Family", "Family_summary.xlsx")