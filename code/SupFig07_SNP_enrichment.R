#-------------------------------------------------------------------------------
# Supplementary Figure 7: Enrichment of NPD-associated SNPs within global cortical morphology GWAS signals
#-------------------------------------------------------------------------------
# Legend: Enrichment of NPD-associated SNPs within global cortical morphology GWAS signals. SNPs are
# the independent genome-wide significant SNPs of the European-ancestry analyses for ADHD (Demontis et
# al., 2023), BD (O'Connell et al., 2025), MDD (Adams et al., 2025) and SCZ (Bigdeli et al., 2026).
# Percentile ranking of NPD risk variants in the cortical GWAS; 50% expected by chance; dashed red
# line, observed median. Significance: median percentile compared with 10,000 random SNP sets matched
# on minor allele frequency (one-sided); * FDR q<0.05 across the four disorders for each metric. Red,
# CT; blue, SA.
# Abbreviations: ADHD=attention deficit hyperactivity disorder; BD=bipolar disorder; CT=cortical
# thickness; FDR=false discovery rate; GWAS=genome-wide association study; MDD=major depressive
# disorder; NPD=neurodevelopmental and psychiatric disorders; SA=surface area; SCZ=schizophrenia;
# SNP=single nucleotide polymorphism.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")

array_NPD = c("ADHD", "BD", "MDD", "SCZ")

## Load data
in_df_plot = read.csv(paste0(data_dir, "data_SupFig07_SNP_percentiles.csv"))
df_enrichment = read.csv(paste0(data_dir, "data_SupFig07_enrichment_tests.csv"))

df_enrichment = df_enrichment[which(df_enrichment[, "Role"] == "primary"), ]

in_df_plot[, "pct"] = 100 * in_df_plot[, "Percentile"]
in_df_plot[, "metric"] = factor(ifelse(in_df_plot[, "Metric"] == "meanCT", "CT", "SA"), levels = c("CT", "SA"))
in_df_plot[, "NPD"] = factor(in_df_plot[, "Disorder"], levels = array_NPD)

df_enrichment[, "metric"] = factor(ifelse(df_enrichment[, "Metric"] == "meanCT", "CT", "SA"), levels = c("CT", "SA"))
df_enrichment[, "NPD"] = factor(df_enrichment[, "Disorder"], levels = array_NPD)

## Stats Supplementary Figure 7: FDR across the four disorders within each metric
df_enrichment[, "FDR"] = NA
for (test_metric in c("CT", "SA")) {
  index_metric = which(df_enrichment[, "metric"] == test_metric)
  p_array = df_enrichment[index_metric, "P_MAFdecile_null"]
  p_array_adj = p.adjust(p_array, method = "fdr")
  df_enrichment[index_metric, "FDR"] = p_array_adj
}
df_enrichment[, "star"] = ifelse(df_enrichment[, "FDR"] < 0.05, "*", "")

df_median = aggregate(pct ~ NPD + metric, in_df_plot, median)
df_enrichment = merge(df_enrichment, df_median, by = c("NPD", "metric"), all.x = TRUE)
df_enrichment = df_enrichment[order(df_enrichment[, "metric"], df_enrichment[, "NPD"]), ]
print(df_enrichment[, c("NPD", "metric", "N_in_Grasby", "Median_percentile", "pct", "P_MAFdecile_null", "FDR", "star")])

## Plot Supplementary Figure 7
col_values_metric_SNP = c("CT" = "#FF3B30", "SA" = "#2340FF")

p_SupFig07 = ggplot(in_df_plot, aes(x = NPD, y = pct)) +
  theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.75, base_fontface = "plain") +
  geom_violin(aes(fill = metric), colour = NA, scale = "width", width = 0.9, alpha = 0.4) +
  geom_point(position = position_jitter(width = 0.08, seed = 1), size = 0.6, alpha = 0.25) +
  geom_boxplot(aes(fill = metric), width = 0.2, outlier.shape = NA, colour = "black") +
  geom_hline(yintercept = 50, linetype = "dotted", linewidth = 0.8) +
  geom_errorbar(data = df_median, aes(ymin = pct, ymax = pct), width = 1, linetype = "dashed", colour = "red", linewidth = 0.9) +
  geom_text(data = df_enrichment, aes(x = NPD, y = 108, label = star), size = 7) +
  facet_grid(metric ~ .) +
  scale_fill_manual(values = col_values_metric_SNP) +
  scale_y_continuous(breaks = c(0, 50, 100), limits = c(0, 112)) +
  labs(x = NULL, y = "SNP Rank Percentile (%)",
       caption = "* FDR < 0.05 (dashed red line = observed median, dotted black = null)") +
  theme(legend.position = "none", strip.text.y = element_text(face = "bold", angle = 270),
        plot.caption = element_text(hjust = 0.5, size = input_base_font_size + 1))

save_panel(p_SupFig07, "SupFig07", 9, 7)

## Save source data
save_source_data(in_df_plot[, c("SNP", "NPD", "metric", "Percentile")], "SupFig07_SNPs")
df_source_SupFig07 = df_enrichment[, c("NPD", "metric", "N_in_Grasby", "Median_percentile", "P_MAFdecile_null", "FDR", "N_null")]
colnames(df_source_SupFig07) = c("NPD", "metric", "n_SNPs", "median_percentile", "pval_perm_MAF_matched", "FDR", "n_perm")
save_source_data(df_source_SupFig07, "SupFig07")
