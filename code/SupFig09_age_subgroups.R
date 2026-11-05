#-------------------------------------------------------------------------------
# Supplementary Figure 9: Adult-young-pediatric age sub-groups
#-------------------------------------------------------------------------------
# Legend: Violin plots comparing the absolute effect sizes for mean CT, total SA, and CT/SA effect
# size ratios for A) NPDs (top row), B) NPDs-medicated (middle row), and C) NPDs-unmedicated (bottom
# row); Adult-Adolescence-Pediatric abbreviations: -y=adolescence/young; -p=pediatric.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)

## Load data
load(paste0(data_dir, "data_fig2_fig3.RData"))

## NPDs, NPDs-medicated and NPDs-unmedicated, by age group
in_df_plot = data_fig2_fig3[which(data_fig2_fig3$Type != "CNV"), ]
in_df_plot$age_group = factor(in_df_plot$age_group, levels = c("Adult", "Young", "Ped"))

## Absolute effect sizes
in_df_plot$abs_es_MeanCT = abs(in_df_plot$es_MeanCT)
in_df_plot$abs_es_TotalSA = abs(in_df_plot$es_TotalSA)

print(table(in_df_plot$Type, in_df_plot$age_group))


## Stats Supplementary Figure 9A-C: pairwise Wilcoxon rank-sum tests between age groups
array_types = c("NPD", "NPD-med", "NPD-umed")
array_pairs = list(c("Adult", "Young"), c("Adult", "Ped"), c("Young", "Ped"))
array_metrics = c("abs_es_MeanCT", "abs_es_TotalSA", "ratio_CT_SA")

df_pairwise_stats = data.frame()
for (input_type in array_types) {
  in_df_type = in_df_plot[in_df_plot$Type == input_type, ]
  for (test_metric in array_metrics) {
    p_array = sapply(array_pairs, function(pair) {
      wilcox.test(in_df_type[in_df_type$age_group == pair[1], test_metric],
                  in_df_type[in_df_type$age_group == pair[2], test_metric])$p.value
    })
    p_array_adj = p.adjust(p_array, method = "fdr")
    df_pairwise_stats = rbind(df_pairwise_stats,
                              data.frame(type = input_type,
                                         metric = test_metric,
                                         comparison = sapply(array_pairs, paste, collapse = "-"),
                                         pval = p_array,
                                         pval_adj = p_array_adj))
  }
}
print(df_pairwise_stats)


## Plot Supplementary Figure 9A-C
col_values_age_group = c("Adult" = "#1b9e77", "Young" = "#d95f02", "Ped" = "#7570b3")
input_label_font_size = 4

violin_plot = function(in_df, y_var, y_lab, y_lim, y_breaks, y_labs) {
  ggplot(in_df, aes(x = age_group, y = .data[[y_var]], color = age_group)) +
    theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.5) +
    geom_violin(draw_quantiles = 0.5, linetype = "dashed", trim = TRUE) +
    scale_color_manual(values = col_values_age_group) +
    geom_point(size = 1.5, position = position_jitter(width = 0.3, seed = 1)) +
    geom_text_repel(aes(label = Column), size = input_label_font_size, min.segment.length = 0.3, max.overlaps = 20, seed = 1) +
    scale_y_continuous(trans = "log10", limits = y_lim, breaks = y_breaks, labels = y_labs) +
    xlab(NULL) + ylab(y_lab) +
    theme(axis.line.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(5, 5, 5, 5))
}

y_lim_es = c(0.00295, 0.85)
y_breaks_es = c(seq(0.003, 0.009, 0.001), seq(0.01, 0.09, 0.01), seq(0.1, 0.8, 0.1))
y_labs_es = c(rep("", 7), 0.01, rep("", 8), 0.1, rep("", 3), 0.5, rep("", 3))

y_lim_ratio = c(0.045, 41.5)
y_breaks_ratio = c(seq(0.05, 0.09, 0.01), seq(0.1, 0.9, 0.1), seq(1, 9, 1), 10, 20, 30, 40)
y_labs_ratio = c(0.05, rep("", 4), 0.1, rep("", 8), 1, rep("", 8), 10, rep("", 3))

array_row_titles = c("A. NPDs", "B. NPDs-medicated", "C. NPDs-unmedicated")
array_panel_names = c("SupFig09A", "SupFig09B", "SupFig09C")
list_rows = list()
for (i in seq_along(array_types)) {
  in_df_type = in_df_plot[in_df_plot$Type == array_types[i], ]
  p_CT = violin_plot(in_df_type, "abs_es_MeanCT", "Effect size (abs)", y_lim_es, y_breaks_es, y_labs_es)
  p_SA = violin_plot(in_df_type, "abs_es_TotalSA", NULL, y_lim_es, y_breaks_es, y_labs_es) +
    theme(axis.text.y = element_blank())
  p_ratio = violin_plot(in_df_type, "ratio_CT_SA", "Cohen's d CT / Cohen's d SA ratio", y_lim_ratio, y_breaks_ratio, y_labs_ratio) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black")
  p_row = plot_grid(p_CT, p_SA, p_ratio, nrow = 1, align = "h")
  p_row_title = ggdraw() + draw_label(array_row_titles[i], fontface = "bold", size = 16, x = 0.005, hjust = 0)
  list_rows[[i]] = plot_grid(p_row_title, p_row, ncol = 1, rel_heights = c(0.07, 1))
  save_panel(list_rows[[i]], array_panel_names[i], 18, 6)
}

p_column_titles = plot_grid(ggdraw() + draw_label("Cortical thickness", size = 16),
                            ggdraw() + draw_label("Surface area", size = 16),
                            ggdraw() + draw_label("Effect size ratio", size = 16), nrow = 1)
p_supfig9 = plot_grid(p_column_titles, list_rows[[1]], list_rows[[2]], list_rows[[3]], ncol = 1,
                      rel_heights = c(0.12, 1, 1, 1))
print(p_supfig9)

save_panel(p_supfig9, "SupFig09", 18, 18)


## Save source data
for (i in seq_along(array_types)) {
  in_df_type = in_df_plot[in_df_plot$Type == array_types[i], ]
  df_source = data.frame(condition = in_df_type$Column, type = in_df_type$Type, age_group = in_df_type$age_group,
                         abs_es_MeanCT = in_df_type$abs_es_MeanCT, abs_es_TotalSA = in_df_type$abs_es_TotalSA,
                         ratio_CT_SA = in_df_type$ratio_CT_SA)
  df_source$n_age_group = as.vector(table(in_df_type$age_group)[as.character(df_source$age_group)])
  df_stats = df_pairwise_stats[df_pairwise_stats$type == array_types[i], ]
  for (j in seq_len(nrow(df_stats))) {
    df_source[, paste0("pval_", df_stats$metric[j], "_", df_stats$comparison[j])] = df_stats$pval[j]
    df_source[, paste0("pval_adj_", df_stats$metric[j], "_", df_stats$comparison[j])] = df_stats$pval_adj[j]
  }
  save_source_data(df_source, array_panel_names[i])
}
