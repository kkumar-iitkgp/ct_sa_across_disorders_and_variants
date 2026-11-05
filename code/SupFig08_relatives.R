#-------------------------------------------------------------------------------
# Supplementary Figure 8: Cortical alterations in psychiatric disorders versus 1st-degree relatives
#-------------------------------------------------------------------------------
# Legend: A) Bar plots comparing mean cortical thickness (CT) and total surface area (SA) effect sizes
# for individuals with bipolar disorder (BD) or schizophrenia (SCZ) versus their unaffected
# first-degree relatives. Data are derived from the ENIGMA Relatives study11. Blue bars indicate the
# disorder group, and orange bars indicate first-degree relatives. # indicates significant
# case-control differences as reported in the original study. B) Violin plots comparing absolute
# effect sizes of relative sub-groups against the disorder group (BD / SCZ). The dashed horizontal
# lines indicate the effect size of the disorder group. Labeled points represent specific relative
# categories. * indicates that the median effect size of the relative sub-groups is significantly
# lower than the disorder effect size (FDR q<0.05; one-sample Wilcoxon signed-rank test). The y-axis
# in panel B represents absolute Cohen's d plotted on a log10 scale.
# Abbreviations: BD=bipolar disorder; CT=cortical thickness; DZ=dizygotic twin; MZ=monozygotic twin;
# SA=surface area; SCZ=schizophrenia.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(patchwork)

col_value_FDRs = c("Disorder" = "#56B4E9", "1st degree relatives" = "#E69F00")
col_values_Relatives = c("BD" = "#4daf4a", "SCZ" = "#4d4d4d")

## Load data: effect sizes reported in the ENIGMA Relatives study (de Zwarte et al.)
# Panel A: disorder and all first-degree relatives
df_es_fig_A = data.frame(disorder = rep(c("BD", "BD", "SCZ", "SCZ"), 2),
                         group_name = rep(c("1st degree relatives", "Disorder"), 4),
                         metric = c(rep("Mean CT", 4), rep("Total SA", 4)),
                         es = c(-0.01, -0.35, -0.13, -0.62,
                                0.15, 0.08, -0.01, -0.21),
                         sig_label = c(NA, "#", NA, "#", NA, NA, NA, NA))

# Panel B: disorder and relative sub-groups
array_group_name = c("Disorder", "MZ", "DZ", "Offspring", "Sibling", "Parent")
array_es_TotalSA = c(0.08, 0.04, 0.05, 0.15, 0.18, 0.22,
                     -0.21, -0.2, -0.04, -0.27, 0.1, 0.03)
array_es_MeanCT = c(-0.35, 0.12, -0.03, -0.07, -0.02, 0.13,
                    -0.62, -0.29, -0.52, -0.07, -0.11, 0.01)

df_Global_ES_Relatives = data.frame(Column = rep(array_group_name, 2),
                                    Type = factor(c(rep("BD", 6), rep("SCZ", 6)), levels = c("BD", "SCZ")),
                                    es_TotalSA = array_es_TotalSA,
                                    es_MeanCT = array_es_MeanCT)
df_Global_ES_Relatives[, "abs_es_MeanCT"] = abs(df_Global_ES_Relatives[, "es_MeanCT"])
df_Global_ES_Relatives[, "abs_es_TotalSA"] = abs(df_Global_ES_Relatives[, "es_TotalSA"])

violin_data = df_Global_ES_Relatives[which(df_Global_ES_Relatives$Column != "Disorder"), ]
disorder_data = df_Global_ES_Relatives[which(df_Global_ES_Relatives$Column == "Disorder"), ]

## Stats Supplementary Figure 8B: relative sub-groups vs disorder (one-sample Wilcoxon), FDR per metric
df_stats_fig_B = data.frame()
for (test_metric in c("abs_es_MeanCT", "abs_es_TotalSA")) {
  p_array = c()
  array_median = c()
  for (input_type in c("BD", "SCZ")) {
    w1 = wilcox.test(x = violin_data[which(violin_data[, "Type"] == input_type), test_metric],
                     mu = disorder_data[which(disorder_data[, "Type"] == input_type), test_metric],
                     alternative = "less")
    p_array = c(p_array, w1$p.value)
    array_median = c(array_median, median(violin_data[which(violin_data[, "Type"] == input_type), test_metric]))
  }
  p_array_adj = p.adjust(p_array, method = "fdr")
  df_stats_fig_B = rbind(df_stats_fig_B,
                         data.frame(metric = test_metric,
                                    Type = c("BD", "SCZ"),
                                    disorder_abs_es = disorder_data[, test_metric],
                                    median_relatives_abs_es = array_median,
                                    n_relative_groups = 5,
                                    pval = p_array,
                                    pval_adj = p_array_adj))
}
df_stats_fig_B[, "sig_label"] = ifelse(df_stats_fig_B[, "pval_adj"] < 0.05, "*", "")
print(df_stats_fig_B)

## Plot Supplementary Figure 8A: bar plots
df_es_fig_A[, "group_name"] = factor(df_es_fig_A[, "group_name"], levels = c("Disorder", "1st degree relatives"))

p_SupFig08A = ggplot(df_es_fig_A, aes(x = disorder, y = es, fill = group_name)) +
  theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.75) +
  geom_bar(stat = "identity", width = 0.4, position = position_dodge(width = 0.5)) +
  geom_text(aes(y = es - 0.05, label = sig_label), position = position_dodge(width = 0.5), size = 7, na.rm = TRUE) +
  facet_wrap(~ metric) +
  scale_fill_manual(values = col_value_FDRs) +
  scale_y_continuous(limits = c(-0.72, 0.2), breaks = c(-0.5, -0.25, 0)) +
  xlab(NULL) + ylab("Effect size") +
  ggtitle("A. Effect sizes: disorder & all 1st degree relatives") +
  theme(legend.position = "bottom", strip.text = element_text(size = input_base_font_size + 2),
        plot.title = element_text(size = input_base_font_size + 2, face = "plain", hjust = 0))

## Plot Supplementary Figure 8B: violins of absolute effect sizes
fPlot_relatives_violin = function(test_metric, in_title) {
  df_stats = df_stats_fig_B[which(df_stats_fig_B[, "metric"] == test_metric), ]
  df_stats[, "x"] = c(1, 2) + 0.3
  df_stats = df_stats[which(df_stats[, "sig_label"] == "*"), ]
  ggplot(violin_data, aes(x = Type, y = .data[[test_metric]], color = Type)) +
    theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.75) +
    geom_violin(draw_quantiles = c(0.5), linetype = "dashed", trim = TRUE) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = Column), size = 5, min.segment.length = 0.3, max.overlaps = 30, seed = 1) +
    geom_segment(data = disorder_data, aes(x = as.numeric(Type) - 0.35, xend = as.numeric(Type) + 0.35,
                                           y = .data[[test_metric]], yend = .data[[test_metric]]), linetype = "dashed", linewidth = 0.4) +
    geom_point(data = disorder_data, size = 2) +
    geom_text_repel(data = disorder_data, aes(label = Column), size = 5.5, colour = "black", fontface = "bold",
                    nudge_y = 0.1, min.segment.length = 1, seed = 1) +
    geom_segment(data = df_stats, aes(x = x, xend = x, y = disorder_abs_es, yend = median_relatives_abs_es), colour = "black", inherit.aes = FALSE) +
    geom_segment(data = df_stats, aes(x = x - 0.04, xend = x, y = disorder_abs_es, yend = disorder_abs_es), colour = "black", inherit.aes = FALSE) +
    geom_segment(data = df_stats, aes(x = x - 0.04, xend = x, y = median_relatives_abs_es, yend = median_relatives_abs_es), colour = "black", inherit.aes = FALSE) +
    geom_text(data = df_stats, aes(x = x + 0.08, y = sqrt(disorder_abs_es * median_relatives_abs_es), label = sig_label),
              colour = "black", size = 8, inherit.aes = FALSE) +
    scale_color_manual(values = col_values_Relatives) +
    scale_y_continuous(trans = "log10", limits = c(0.009, 1.01), breaks = c(seq(0.01, 0.09, 0.01), seq(0.1, 1, 0.1)),
                       labels = c(0.01, rep("", 8), 0.1, rep("", 8), 1)) +
    xlab(NULL) + ylab("Effect size (absolute Cohen's d)") + ggtitle(in_title) +
    theme(axis.line.x = element_blank(), legend.position = "none",
          plot.title = element_text(size = input_base_font_size + 2, face = "plain"))
}

p_SupFig08B_CT = fPlot_relatives_violin("abs_es_MeanCT", "Mean CT")
p_SupFig08B_SA = fPlot_relatives_violin("abs_es_TotalSA", "Total SA") + ylab(NULL)
p_SupFig08B = (p_SupFig08B_CT | p_SupFig08B_SA) +
  plot_annotation(title = "B. Effect sizes: disorder & 1st degree relative sub-groups",
                  theme = theme(plot.title = element_text(size = input_base_font_size + 2)))

save_panel(p_SupFig08A, "SupFig08A", 6.5, 6.5)
save_panel(p_SupFig08B, "SupFig08B", 8, 6.5)

p_SupFig08 = wrap_elements(p_SupFig08A) | wrap_elements(p_SupFig08B)
p_SupFig08 = p_SupFig08 + plot_layout(widths = c(0.45, 0.55))
save_panel(p_SupFig08, "SupFig08", 15, 7)

## Save source data
save_source_data(df_es_fig_A, "SupFig08A")
df_source_SupFig08B = merge(df_Global_ES_Relatives,
                            reshape(df_stats_fig_B[, c("metric", "Type", "pval", "pval_adj")], idvar = "Type", timevar = "metric", direction = "wide"),
                            by = "Type")
save_source_data(df_source_SupFig08B, "SupFig08B")
