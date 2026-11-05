#-------------------------------------------------------------------------------
# Supplementary Figure 4: Deletion-Duplication-SCA stratifications
#-------------------------------------------------------------------------------
# Legend: Violin plots comparing the absolute effect sizes for A) mean CT and B) total SA; for
# deletions (Del), duplications (Dup), sex chromosome aneuploidies (SCA); and C) CT/SA effect size
# ratios.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)

## Load data
load(paste0(data_dir, "data_fig2_fig3.RData"))

## CNVs and aneuploidies, split into deletions, duplications and sex chromosome aneuploidies (SCA)
in_df_plot = data_fig2_fig3[which(data_fig2_fig3$Type == "CNV"), ]
in_df_plot[in_df_plot$Column %in% c("XXX", "XXY", "XYY", "XXYY", "X (TS)"), "CNV_NPD_Type"] = "SCA"
in_df_plot$CNV_NPD_Type = factor(in_df_plot$CNV_NPD_Type, levels = c("Del", "Dup", "SCA"))

## Absolute effect sizes
in_df_plot$abs_es_MeanCT = abs(in_df_plot$es_MeanCT)
in_df_plot$abs_es_TotalSA = abs(in_df_plot$es_TotalSA)

print(table(in_df_plot$CNV_NPD_Type))


## Stats Supplementary Figure 4A-C: pairwise Wilcoxon rank-sum tests between Del, Dup and SCA
array_pairs = list(c("Del", "Dup"), c("Del", "SCA"), c("Dup", "SCA"))
array_metrics = c("abs_es_MeanCT", "abs_es_TotalSA", "ratio_CT_SA")

df_pairwise_stats = data.frame()
for (test_metric in array_metrics) {
  p_array = sapply(array_pairs, function(pair) {
    wilcox.test(in_df_plot[in_df_plot$CNV_NPD_Type == pair[1], test_metric],
                in_df_plot[in_df_plot$CNV_NPD_Type == pair[2], test_metric])$p.value
  })
  p_array_adj = p.adjust(p_array, method = "fdr")
  df_pairwise_stats = rbind(df_pairwise_stats,
                            data.frame(metric = test_metric,
                                       comparison = sapply(array_pairs, paste, collapse = "-"),
                                       pval = p_array,
                                       pval_adj = p_array_adj))
}
print(df_pairwise_stats)


## Stats Supplementary Figure 4C: shift of the median ratio from 1
test_metric = "ratio_CT_SA"
array_types = c("Del", "Dup", "SCA")
p_array = sapply(array_types, function(t) {
  wilcox.test(in_df_plot[in_df_plot$CNV_NPD_Type == t, test_metric], mu = 1, alternative = "less")$p.value
})
p_array_adj = p.adjust(p_array, method = "fdr")

df_median_es_ratio_supfig4C = data.frame(type = array_types,
                                         n = as.vector(table(in_df_plot$CNV_NPD_Type)[array_types]),
                                         median_ratio = sapply(array_types, function(t) median(in_df_plot[in_df_plot$CNV_NPD_Type == t, test_metric])),
                                         pval = p_array,
                                         pval_adj = p_array_adj,
                                         row.names = NULL)
print(df_median_es_ratio_supfig4C)


## Plot Supplementary Figure 4A-C
col_values_Del_Dup_SCA = c("Del" = "#e41a1c", "Dup" = "#377eb8", "SCA" = "#4d4d4d")
input_label_font_size = 4

violin_plot = function(y_var, y_lab, y_lim, y_breaks, y_labs, plot_title) {
  ggplot(in_df_plot, aes(x = CNV_NPD_Type, y = .data[[y_var]], color = CNV_NPD_Type)) +
    theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.5) +
    geom_violin(draw_quantiles = 0.5, linetype = "dashed", trim = TRUE) +
    scale_color_manual(values = col_values_Del_Dup_SCA) +
    geom_point(size = 1.5, position = position_jitter(width = 0.3, seed = 1)) +
    geom_text_repel(aes(label = Column), size = input_label_font_size, min.segment.length = 0.3, max.overlaps = 20, seed = 1) +
    scale_y_continuous(trans = "log10", limits = y_lim, breaks = y_breaks, labels = y_labs) +
    xlab(NULL) + ylab(y_lab) +
    ggtitle(plot_title) +
    theme(axis.line.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "plain", hjust = 0),
          plot.margin = margin(5, 5, 5, 5))
}

y_lim_es = c(0.05, 2.65)
y_breaks_es = c(seq(0.05, 0.09, 0.01), seq(0.1, 0.9, 0.1), 1, 2)
y_labs_es = c(rep("", 5), 0.1, rep("", 8), 1, "")

p_supfig4A = violin_plot("abs_es_MeanCT", "Effect size (abs)", y_lim_es, y_breaks_es, y_labs_es, "A. Cortical thickness")
p_supfig4B = violin_plot("abs_es_TotalSA", NULL, y_lim_es, y_breaks_es, y_labs_es, "B. Surface area") +
  theme(axis.text.y = element_blank())
p_supfig4C = violin_plot("ratio_CT_SA", "Cohen's d CT / Cohen's d SA ratio", c(0.045, 3.5),
                         c(seq(0.05, 0.09, 0.01), seq(0.1, 0.9, 0.1), seq(1, 3, 1)),
                         c(0.05, rep("", 4), 0.1, rep("", 8), 1, rep("", 2)), "C. Effect size ratio") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black")

p_supfig4 = plot_grid(p_supfig4A, p_supfig4B, p_supfig4C, nrow = 1, align = "h")
print(p_supfig4)

save_panel(p_supfig4A, "SupFig04A", 6, 5)
save_panel(p_supfig4B, "SupFig04B", 6, 5)
save_panel(p_supfig4C, "SupFig04C", 6, 5)
save_panel(p_supfig4, "SupFig04", 18, 5.5)


## Save source data
source_data_panel = function(y_var, test_metric) {
  out = data.frame(condition = in_df_plot$Column, group = in_df_plot$CNV_NPD_Type, value = in_df_plot[, y_var])
  names(out)[3] = y_var
  out$n_group = as.vector(table(in_df_plot$CNV_NPD_Type)[as.character(out$group)])
  df_stats = df_pairwise_stats[df_pairwise_stats$metric == test_metric, ]
  for (i in seq_len(nrow(df_stats))) {
    out[, paste0("pval_", df_stats$comparison[i])] = df_stats$pval[i]
    out[, paste0("pval_adj_", df_stats$comparison[i])] = df_stats$pval_adj[i]
  }
  out
}

save_source_data(source_data_panel("abs_es_MeanCT", "abs_es_MeanCT"), "SupFig04A")
save_source_data(source_data_panel("abs_es_TotalSA", "abs_es_TotalSA"), "SupFig04B")
df_supfig4C = source_data_panel("ratio_CT_SA", "ratio_CT_SA")
df_supfig4C = cbind(df_supfig4C, df_median_es_ratio_supfig4C[match(df_supfig4C$group, df_median_es_ratio_supfig4C$type),
                                                              c("median_ratio", "pval", "pval_adj")])
names(df_supfig4C)[(ncol(df_supfig4C) - 2):ncol(df_supfig4C)] = c("median_ratio_group", "pval_ratio_vs_1", "pval_adj_ratio_vs_1")
save_source_data(df_supfig4C, "SupFig04C")
