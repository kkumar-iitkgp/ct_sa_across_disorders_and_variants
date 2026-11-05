#-------------------------------------------------------------------------------
# Supplementary Figure 2: Concordance between CT and SA effect sizes across NPDs, genetic variants,
# and medication status
#-------------------------------------------------------------------------------
# Legend: A-B) Scatter plots comparing mean cortical thickness (CT) and total surface area (SA)
# effect sizes for A) rare genetic variants (CNVs); and B) neurodevelopmental and psychiatric
# disorders (NPDs). C-D) Comparison of CT and SA effect sizes for NPD diagnostic sub-groups split by
# medication status: C) Medicated; and D) Un-medicated. Each point represents the signed Cohen's d
# value, where positive and negative values indicate increased and decreased metrics in cases
# compared to controls, respectively. The blue solid line represents the linear regression fit. The
# gray dot-dash line (y=x) indicates the line of perfect concordance where CT and SA effect sizes
# are equal. Dotted lines at x=0 and y=0 delineate the quadrants of directionality. Marginal density
# plots (top and right) visualize the distribution of effect sizes for SA and CT, respectively.
# Arrows within the marginal plots indicate the direction and magnitude of the mean effect size
# shift from zero. Case-control differences were adjusted for age, sex, and site. Axes represent
# signed effect sizes.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
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


## Stats Supplementary Figure 2: linear fit (points within the axis limits) and mean signed effect sizes
array_types = c("CNV", "NPD", "NPD-med", "NPD-umed")
list_limits = list("CNV" = c(-2, 1.5), "NPD" = c(-0.6, 0.3), "NPD-med" = c(-0.8, 0.5), "NPD-umed" = c(-0.3, 0.2))
in_df_plot = rbind(in_df_fig2, in_df_fig3)
in_df_plot$within_axis_limits = NA
for (t in array_types) {
  idx = which(in_df_plot$Type == t)
  in_df_plot$within_axis_limits[idx] = (in_df_plot$es_TotalSA[idx] >= list_limits[[t]][1]) & (in_df_plot$es_TotalSA[idx] <= list_limits[[t]][2]) &
    (in_df_plot$es_MeanCT[idx] >= list_limits[[t]][1]) & (in_df_plot$es_MeanCT[idx] <= list_limits[[t]][2])
}

df_stats_supfig2 = do.call(rbind, lapply(array_types, function(t) {
  df = in_df_plot[in_df_plot$Type == t, ]
  fit = lm(es_MeanCT ~ es_TotalSA, data = df[df$within_axis_limits, ])
  data.frame(type = t, n = nrow(df), n_outside_axis_limits = sum(!df$within_axis_limits),
             lm_intercept = coef(fit)[1], lm_slope = coef(fit)[2],
             mean_es_TotalSA = mean(df$es_TotalSA), mean_es_MeanCT = mean(df$es_MeanCT), row.names = NULL)
}))
print(df_stats_supfig2)


## Plot Supplementary Figure 2A-D: scatter plot with marginal densities
scatter_panel = function(df, plot_title, pt_color, axis_lims) {

  start_point = axis_lims[1]
  label_pos = start_point + (axis_lims[2] - start_point) * 0.1

  p_main = ggplot(df, aes(x = es_TotalSA, y = es_MeanCT)) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotdash", color = "gray50", linewidth = 0.8) +
    annotate("text", x = label_pos, y = label_pos, label = "y=x", angle = 45, vjust = -0.5, color = "gray50", size = 3.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "blue", linewidth = 1) +
    geom_point(color = pt_color, size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = Column), size = 3, max.overlaps = 20, box.padding = 0.3) +
    scale_x_continuous(limits = axis_lims) +
    scale_y_continuous(limits = axis_lims) +
    theme_prism(base_size = input_base_font_size, base_line_size = 0.5) +
    labs(title = plot_title, x = "SA Effect Size", y = "CT Effect Size") +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"), legend.position = "none")

  mean_sa = mean(df$es_TotalSA)
  mean_ct = mean(df$es_MeanCT)
  arrow_y_sa = max(density(df$es_TotalSA)$y) * 0.5
  arrow_y_ct = max(density(df$es_MeanCT)$y) * 0.5

  p_top = ggplot(df, aes(x = es_TotalSA)) +
    geom_density(fill = pt_color, alpha = 0.4, color = NA) +
    geom_vline(xintercept = mean_sa, linetype = "dashed", alpha = 0.5) +
    annotate("segment", x = 0, xend = mean_sa, y = arrow_y_sa, yend = arrow_y_sa,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.8) +
    geom_rug(color = pt_color, alpha = 0.7) +
    scale_x_continuous(limits = axis_lims) +
    theme_void() +
    theme(plot.margin = margin(b = -5))

  p_right = ggplot(df, aes(x = es_MeanCT)) +
    geom_density(fill = pt_color, alpha = 0.4, color = NA) +
    geom_vline(xintercept = mean_ct, linetype = "dashed", alpha = 0.5) +
    annotate("segment", x = 0, xend = mean_ct, y = arrow_y_ct, yend = arrow_y_ct,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 0.8) +
    geom_rug(color = pt_color, alpha = 0.7) +
    scale_x_continuous(limits = axis_lims) +
    coord_flip() +
    theme_void() +
    theme(plot.margin = margin(l = -5))

  p_combined = insert_xaxis_grob(p_main, p_top, position = "top", height = unit(0.2, "null"))
  insert_yaxis_grob(p_combined, p_right, position = "right", width = unit(0.2, "null"))
}

array_titles = c("CNV" = "A. Rare genetic variants (CNVs)", "NPD" = "B. NPD diagnosis",
                 "NPD-med" = "C. NPD medicated", "NPD-umed" = "D. NPD un-medicated")
list_panels = lapply(array_types, function(t) {
  df = in_df_plot[in_df_plot$Type == t, ]
  scatter_panel(df, array_titles[t], col_values_CNV_NPD[t], list_limits[[t]])
})
names(list_panels) = array_types

title_1 = textGrob("NPDs and associated rare genetic variants (CNVs)", gp = gpar(fontsize = 14, fontface = "bold"))
title_2 = textGrob("Influence of non-genetic factors (Medication status)", gp = gpar(fontsize = 14, fontface = "bold"))
p_supfig2 = arrangeGrob(title_1, arrangeGrob(list_panels[["CNV"]], list_panels[["NPD"]], ncol = 2),
                        title_2, arrangeGrob(list_panels[["NPD-med"]], list_panels[["NPD-umed"]], ncol = 2),
                        ncol = 1, heights = c(0.1, 1, 0.1, 1))

save_panel(list_panels[["CNV"]], "SupFig02A", 6, 5.5)
save_panel(list_panels[["NPD"]], "SupFig02B", 6, 5.5)
save_panel(list_panels[["NPD-med"]], "SupFig02C", 6, 5.5)
save_panel(list_panels[["NPD-umed"]], "SupFig02D", 6, 5.5)
save_panel(p_supfig2, "SupFig02", 12, 12)


## Save source data
array_panels = c("CNV" = "SupFig02A", "NPD" = "SupFig02B", "NPD-med" = "SupFig02C", "NPD-umed" = "SupFig02D")
for (t in array_types) {
  df = in_df_plot[in_df_plot$Type == t, c("Column", "Type", "es_TotalSA", "es_MeanCT", "within_axis_limits")]
  df = cbind(df, df_stats_supfig2[df_stats_supfig2$type == t, -1], row.names = NULL)
  save_source_data(df, array_panels[t])
}
