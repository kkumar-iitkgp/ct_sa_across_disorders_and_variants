#-------------------------------------------------------------------------------
# Supplementary Figure 6: Genetic correlations between global cortical metrics and NPDs from Grasby et al. 2020
#-------------------------------------------------------------------------------
# Legend: Genetic correlations (rg) between neurodevelopmental and psychiatric disorders (NPDs) and
# global cortical metrics: mean cortical thickness (CT) and total surface area (SA). Correlations and
# standard errors were derived using summary statistics from Grasby et al. 10. * indicates FDR
# significant (q<0.1). Note: We used the lenient FDR q < 0.1 threshold, as these values come from
# Grasby et al. (2020), where FDR was applied across multiple psychiatric disorders and MRI metrics.
# Abbreviations: ADHD=attention deficit hyperactivity disorder; BD=bipolar disorder; CT=cortical
# thickness; FDR=false discovery rate; MDD=major depressive disorder; NPD=neurodevelopmental and
# psychiatric disorders; SA=surface area; SCZ=schizophrenia.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")

## Load data: genetic correlations reported in Grasby et al. 2020
df_rg_grasby = data.frame(metric = c(rep("CT", 4), rep("SA", 4)),
                          disorder = c("ADHD", "BD", "MDD", "SCZ", "ADHD", "BD", "MDD", "SCZ"),
                          rg = c(-0.0238, 0.0402, 0.0036, 0.0081,
                                 -0.1698, 0.0662, -0.0861, -0.0234),
                          rg_se = c(0.0462, 0.044, 0.0332, 0.031,
                                    0.039, 0.0375, 0.03, 0.0294),
                          rg_pval = c(0.607, 0.3602, 0.9126, 0.7941,
                                      0.000013347, 0.0778, 0.0041, 0.4258))

## Stats Supplementary Figure 6: significance at q < 0.1
in_df_plot = df_rg_grasby
in_df_plot[, "metric"] = factor(in_df_plot[, "metric"], levels = c("CT", "SA"))
in_df_plot[, "sig_label"] = ifelse(in_df_plot[, "rg_pval"] < 0.1, "*", NA)
in_df_plot[, "star_y"] = ifelse(in_df_plot[, "rg"] >= 0, in_df_plot[, "rg"] + in_df_plot[, "rg_se"] + 0.0001, in_df_plot[, "rg"] - in_df_plot[, "rg_se"] - 0.005)
in_df_plot[, "vjust_val"] = ifelse(in_df_plot[, "rg"] >= 0, 0.5, 1)
print(in_df_plot[, c("metric", "disorder", "rg", "rg_se", "rg_pval", "sig_label")])

## Plot Supplementary Figure 6
p_SupFig06 = ggplot(in_df_plot, aes(x = disorder, y = rg, fill = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_col(color = "black", alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = rg - rg_se, ymax = rg + rg_se), width = 0.2) +
  geom_text(aes(y = star_y, label = sig_label, vjust = vjust_val), size = 8, na.rm = TRUE) +
  facet_grid(metric ~ .) +
  scale_fill_manual(values = c("CT" = "red", "SA" = "blue")) +
  theme_prism(base_size = 18, base_line_size = 0.75, base_rect_size = 0.75, base_fontface = "plain") +
  labs(x = NULL, y = "Genetic correlation (rg)", caption = "* FDR-q < 0.1") +
  theme(legend.position = "none")

save_panel(p_SupFig06, "SupFig06", 8, 6)

## Save source data
save_source_data(in_df_plot[, c("metric", "disorder", "rg", "rg_se", "rg_pval", "sig_label")], "SupFig06")
