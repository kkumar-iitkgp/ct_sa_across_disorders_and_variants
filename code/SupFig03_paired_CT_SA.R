#-------------------------------------------------------------------------------
# Supplementary Figure 3: Paired comparison of CT and SA effect sizes
#-------------------------------------------------------------------------------
# Legend: A-B) Paired boxplots comparing absolute effect sizes of mean cortical thickness (CT) and
# total surface area (SA) for A) rare genetic variants (CNVs, paired Wilcoxon FDR q=2E-2); and B)
# neurodevelopmental and psychiatric disorders (NPDs, paired Wilcoxon FDR q=1.6E-1). C-D) Comparison
# of absolute CT and SA effect sizes for NPD diagnostic sub-groups stratified by medication status:
# C) Medicated (paired Wilcoxon FDR q=1.6E-2); and D) Un-medicated (paired Wilcoxon FDR q=1). Case-
# control differences were calculated after adjusting for age, sex, and site. Individual points
# represent specific CNVs or disorders, with gray lines connecting paired values to visualize the
# within-disorder relationship between metrics. Boxplots represent the median and interquartile
# range (IQR). Red indicates CT and blue indicates SA. The y-axis represents absolute Cohen's d
# values plotted on a log10 scale.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(ggpubr)
library(cowplot)
library(gridExtra)
library(grid)

## Load data
load(paste0(data_dir, "data_fig2_fig3.RData"))

## CNVs and adult NPDs, plus conduct disorder (CD)
in_df_fig2 = data_fig2_fig3[which((data_fig2_fig3$Type %in% c("CNV", "NPD")) & (data_fig2_fig3$age_group %in% c("CNV", "Adult"))), ]
in_df_fig2 = rbind(in_df_fig2, data_fig2_fig3[which(data_fig2_fig3$Column == "CD"), ])

## Adult NPD sub-groups with and without medication
in_df_fig3 = data_fig2_fig3[which((data_fig2_fig3$Type %in% c("NPD-med", "NPD-umed")) & (data_fig2_fig3$age_group == "Adult")), ]

in_df_plot = rbind(in_df_fig2, in_df_fig3)
in_df_plot$abs_es_MeanCT = abs(in_df_plot$es_MeanCT)
in_df_plot$abs_es_TotalSA = abs(in_df_plot$es_TotalSA)

## Long format for paired plots
in_df_plot_long = data.frame(Column = rep(in_df_plot$Column, 2),
                             Type = rep(in_df_plot$Type, 2),
                             metric = rep(c("CT", "SA"), each = nrow(in_df_plot)),
                             ES = c(in_df_plot$abs_es_MeanCT, in_df_plot$abs_es_TotalSA))


## Stats Supplementary Figure 3: paired Wilcoxon signed-rank test, CT vs SA
array_types = c("CNV", "NPD", "NPD-med", "NPD-umed")
list_tests = lapply(array_types, function(t) {
  df = in_df_plot[in_df_plot$Type == t, ]
  wilcox.test(df$abs_es_MeanCT, df$abs_es_TotalSA, paired = TRUE)
})
p_array = sapply(list_tests, function(w) w$p.value)
p_array_adj = p.adjust(p_array, method = "fdr")

df_paired_CT_SA_supfig3 = data.frame(type = array_types,
                                     n_pairs = as.vector(table(in_df_plot$Type)[array_types]),
                                     statistic = sapply(list_tests, function(w) w$statistic),
                                     pval = p_array,
                                     pval_adj = p_array_adj)
print(df_paired_CT_SA_supfig3)


## Plot Supplementary Figure 3A-D: paired boxplots
names(col_values_metric) = c("CT", "SA")

paired_boxplot = function(df, plot_title, q, show_ylab) {
  df$metric = factor(df$metric, levels = c("CT", "SA"))
  p = ggpaired(df, x = "metric", y = "ES", color = "metric", shape = "metric",
               line.color = "gray", line.size = 0.6, palette = "jco") +
    theme_prism(base_size = 14, base_line_size = 0.5) +
    scale_color_manual(values = col_values_metric) +
    scale_y_continuous(trans = "log10", limits = c(0.002, 3),
                       breaks = c(seq(0.002, 0.009, 0.001), seq(0.01, 0.09, 0.01), seq(0.1, 0.9, 0.1), 1, 2, 3),
                       labels = c(rep("", 8), 0.01, rep("", 8), 0.1, rep("", 8), 1, "", ""),
                       oob = scales::oob_keep) +
    coord_cartesian(clip = "off") +
    theme(axis.line.x = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 15))) +
    xlab(NULL) + ylab(if (show_ylab) "Effect size (abs)" else NULL) +
    ggtitle(plot_title)
  if (q < 0.05) {
    p = p +
      annotate("segment", x = 1, xend = 2, y = 4.3, yend = 4.3, linewidth = 0.6) +
      annotate("segment", x = c(1, 2), xend = c(1, 2), y = 4.3, yend = 3.7, linewidth = 0.6) +
      annotate("text", x = 1.5, y = 4.9, label = "*", size = 9)
  }
  p
}

array_titles = c("CNV" = "Rare genetic\nvariants (CNVs)", "NPD" = "NPD\ndiagnosis",
                 "NPD-med" = "NPD\nmedicated", "NPD-umed" = "NPD\nun-medicated")
list_panels = lapply(seq_along(array_types), function(i) {
  paired_boxplot(in_df_plot_long[in_df_plot_long$Type == array_types[i], ], array_titles[i],
                 p_array_adj[i], show_ylab = i %in% c(1, 3))
})

p_legend = ggplot(in_df_plot_long, aes(x = metric, y = ES, color = metric)) +
  geom_point() + scale_color_manual(values = col_values_metric) +
  theme_prism(base_size = 18) +
  theme(legend.position = "bottom", legend.title = element_blank())

title_1 = textGrob("NPDs and associated rare genetic variants (CNVs)", gp = gpar(fontsize = 14, fontface = "bold"))
title_2 = textGrob("Influence of non-genetic factors (Medication status)", gp = gpar(fontsize = 14, fontface = "bold"))
row_1 = plot_grid(list_panels[[1]], list_panels[[2]], labels = c("A.", "B."), ncol = 2, align = "h")
row_2 = plot_grid(list_panels[[3]], list_panels[[4]], labels = c("C.", "D."), ncol = 2, align = "h")
p_supfig3 = arrangeGrob(title_1, row_1, title_2, row_2, get_legend(p_legend), ncol = 1,
                        heights = c(0.1, 1, 0.1, 1, 0.1))

save_panel(list_panels[[1]], "SupFig03A", 4, 4.5)
save_panel(list_panels[[2]], "SupFig03B", 4, 4.5)
save_panel(list_panels[[3]], "SupFig03C", 4, 4.5)
save_panel(list_panels[[4]], "SupFig03D", 4, 4.5)
save_panel(p_supfig3, "SupFig03", 8, 10)


## Save source data
array_panels = c("SupFig03A", "SupFig03B", "SupFig03C", "SupFig03D")
for (i in seq_along(array_types)) {
  df = in_df_plot[in_df_plot$Type == array_types[i], c("Column", "Type", "abs_es_MeanCT", "abs_es_TotalSA")]
  df = cbind(df, df_paired_CT_SA_supfig3[i, -1], row.names = NULL)
  save_source_data(df, array_panels[i])
}
