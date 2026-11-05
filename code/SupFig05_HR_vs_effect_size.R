#-------------------------------------------------------------------------------
# Supplementary Figure 5: Psychiatric risk of CNVs versus their effects on global cortical metrics
#-------------------------------------------------------------------------------
# Legend: Hazard ratios for any neurodevelopmental and psychiatric disorder (NPD) from the Danish
# iPSYCH2015 dataset plotted against absolute effect sizes (Cohen's d) of the same CNVs (n = 9) on
# A) total surface area (SA) and B) mean cortical thickness (CT). Circles indicate deletions and
# triangles duplications. The dashed grey line shows the linear fit. Spearman correlation
# coefficients (rho) and p-values are shown on each panel.
# Abbreviations: CNV=copy number variant; CT=cortical thickness; Del=deletion; Dup=duplication;
# NPD=neurodevelopmental and psychiatric disorders; SA=surface area.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(patchwork)

## Load data
in_df_plot = read.csv(paste0(data_dir, "data_SupFig05_HR_abs_es_CNVs.csv"))
in_df_plot[, "Type_Del_Dup"] = ifelse(grepl("Del$", in_df_plot[, "CNV"]), "Del", "Dup")

## Stats Supplementary Figure 5A and 5B: Spearman correlations
array_metric = c("abs_es_TotalSA", "abs_es_MeanCT")
array_rho = c()
p_array = c()
for (test_metric in array_metric) {
  res_cor = cor.test(in_df_plot[, "HR_iPSYCH2015_ANY"], in_df_plot[, test_metric], method = "spearman", exact = FALSE)
  array_rho = c(array_rho, unname(res_cor$estimate))
  p_array = c(p_array, res_cor$p.value)
}

df_cor_SupFig05 = data.frame(metric = c("totalSA", "meanCT"),
                             n_CNVs = nrow(in_df_plot),
                             spearman_rho = array_rho,
                             pval = p_array)
print(df_cor_SupFig05)

## Plot Supplementary Figure 5A and 5B
fPlot_HR_vs_es = function(test_metric, in_col, in_ylab, in_rho, in_pval) {
  ggplot(in_df_plot, aes(x = HR_iPSYCH2015_ANY, y = .data[[test_metric]])) +
    theme_prism(base_size = input_base_font_size, base_line_size = 0.75) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey60", linetype = "dashed", linewidth = 0.6) +
    geom_point(aes(shape = Type_Del_Dup), colour = in_col, size = 3) +
    geom_text_repel(aes(label = CNV), size = 3.5, colour = "grey20", seed = 1, max.overlaps = 20) +
    scale_shape_manual(values = c("Del" = 16, "Dup" = 17)) +
    xlab("Hazard ratio, any NPD (iPSYCH2015)") + ylab(in_ylab) +
    labs(subtitle = sprintf("Spearman rho = %.2f, p = %s", in_rho, formatC(in_pval, digits = 2, format = "g"))) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0), plot.subtitle = element_text(size = input_base_font_size, hjust = 0))
}

p_SupFig05A = fPlot_HR_vs_es("abs_es_TotalSA", col_SA, "Effect size SA (absolute Cohen's d)", array_rho[1], p_array[1]) + ggtitle("A")
p_SupFig05B = fPlot_HR_vs_es("abs_es_MeanCT", col_CT, "Effect size CT (absolute Cohen's d)", array_rho[2], p_array[2]) + ggtitle("B")

save_panel(p_SupFig05A, "SupFig05A", 5.5, 5)
save_panel(p_SupFig05B, "SupFig05B", 5.5, 5)

p_SupFig05 = p_SupFig05A + p_SupFig05B
save_panel(p_SupFig05, "SupFig05", 11, 5)

## Save source data
df_source_SupFig05A = data.frame(CNV = in_df_plot[, "CNV"],
                                 type = in_df_plot[, "Type_Del_Dup"],
                                 HR_iPSYCH2015_any_NPD = in_df_plot[, "HR_iPSYCH2015_ANY"],
                                 abs_es_TotalSA = in_df_plot[, "abs_es_TotalSA"],
                                 n_CNVs = nrow(in_df_plot),
                                 spearman_rho = array_rho[1],
                                 pval = p_array[1])
df_source_SupFig05B = data.frame(CNV = in_df_plot[, "CNV"],
                                 type = in_df_plot[, "Type_Del_Dup"],
                                 HR_iPSYCH2015_any_NPD = in_df_plot[, "HR_iPSYCH2015_ANY"],
                                 abs_es_MeanCT = in_df_plot[, "abs_es_MeanCT"],
                                 n_CNVs = nrow(in_df_plot),
                                 spearman_rho = array_rho[2],
                                 pval = p_array[2])
save_source_data(df_source_SupFig05A, "SupFig05A")
save_source_data(df_source_SupFig05B, "SupFig05B")
