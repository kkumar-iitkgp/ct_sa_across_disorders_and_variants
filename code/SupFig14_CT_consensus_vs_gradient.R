#-------------------------------------------------------------------------------
# Supplementary Figure 14: Correlations between cortical gradient and CT consensus maps
#-------------------------------------------------------------------------------
# Legend: Correlation between cortical gradient and cortical thickness (CT) consensus maps for A) mean
# absolute effect sizes, B) % significance, and C) standard deviation of effect sizes. *: spin
# permutation p-value < 0.05, which preserves the spatial structure of the cortex while generating null
# distributions to account for spatial autocorrelation.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)
library(grid)

hemi = "CT"
hemi_label = "Thickness (CT)"
fig_name = "SupFig14"

## Load data
## Consensus maps and cortical gradient (-1 x PC1 of AHBA gene expression, abagen)
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))
gradient_map = array_cortical_gradient

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

## Standard deviation of effect sizes per region (stored as "variance" in the data file)
list_consensus = list(meanAbsES = df_concensus_profiles_meanAbsES,
                      significance = df_concensus_profiles_significance,
                      SD = df_concensus_profiles_variance)

## Stats: Pearson correlation with the gradient and spin permutation p-value
df_stats = data.frame()
for (measure in names(list_consensus)) {
  df_consensus = list_consensus[[measure]]
  for (map_name in grep(paste0("^", hemi, "_"), colnames(df_consensus), value = TRUE)) {
    df_stats = rbind(df_stats, data.frame(measure = measure,
                                          map = map_name,
                                          r = cor(df_consensus[, map_name], gradient_map),
                                          p_spin = perm.sphere.p(df_consensus[, map_name], gradient_map, df_perm_ids_LH)))
  }
}
df_stats$label = paste0("r = ", round(df_stats$r, 2), ifelse(df_stats$p_spin < 0.05, "*", ""))
print(df_stats)

## Plot functions
meanES_colors = colorRampPalette(c("#fcfbfd", "#efedf5", "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#54278f", "#3f007d"))
sig_colors = colorRampPalette(c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d"))
sd_colors = colorRampPalette(c("#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"))

## Left hemisphere brain map (lateral and medial views)
fBrainPlot_LH = function(array_vals, color_fn, cbar_min, cbar_max) {
  array_vals = pmin(pmax(array_vals, cbar_min), cbar_max)
  df_brain = data.frame(region = names(dk$palette), value = c(array_vals[1:3], NA, array_vals[4:34]))
  p = ggseg(df_brain, atlas = "dk", position = "dispersed", hemisphere = "left",
            mapping = aes(fill = value), colour = "grey60", size = 0.3) +
    scale_fill_gradientn(limits = c(cbar_min, cbar_max), colours = color_fn(100), na.value = "grey90") +
    theme_void() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
  return(ggplotGrob(p))
}

## Spread label positions so that neighbouring brain maps do not overlap
fDodge = function(x, min_gap = 0.6, x_lim = 1) {
  ord = order(x)
  x_disp = x[ord]
  for (iter in 1:50) {
    for (i in seq_along(x_disp)[-1]) {
      if (x_disp[i] - x_disp[i - 1] < min_gap) {
        x_mid = (x_disp[i] + x_disp[i - 1]) / 2
        x_disp[i - 1] = x_mid - min_gap / 2
        x_disp[i] = x_mid + min_gap / 2
      }
    }
    x_disp = pmin(pmax(x_disp, -x_lim), x_lim)
  }
  x_out = x
  x_out[ord] = x_disp
  return(x_out)
}

## Background: sensorimotor (blue) to association (yellow)
background_colours = c(colorRampPalette(c("#c6d4ef", "#a8bee6", "#c3d1ee", "#f7f8fc"))(130),
                       colorRampPalette(c("#fffcf2", "#fdeec0", "#fde293", "#fde9b5"))(130))

## Gradient axis with consensus maps above the axis
fGradientAxisPanel = function(in_df_items, list_grobs, y_bottom, y_top) {
  p = ggplot() + theme_void() +
    annotation_raster(matrix(background_colours, nrow = 1), xmin = -1.3, xmax = 1.3,
                      ymin = y_bottom, ymax = y_top, interpolate = TRUE) +
    annotate("rect", xmin = -0.012, xmax = 0.012, ymin = y_bottom, ymax = y_top, fill = "#a6a6a6") +
    annotate("segment", x = -1.12, xend = 1.12, y = 0, yend = 0, linewidth = 1.2, colour = "grey30",
             arrow = arrow(length = unit(0.12, "inches"), type = "closed", ends = "both")) +
    annotate("segment", x = c(-1, -0.5, 0, 0.5, 1), xend = c(-1, -0.5, 0, 0.5, 1),
             y = c(-0.15, -0.1, -0.15, -0.1, -0.15), yend = c(0.15, 0.1, 0.15, 0.1, 0.15),
             colour = c("black", "grey50", "black", "grey50", "black")) +
    annotate("text", x = c(-1, 0, 1), y = -0.45, label = c("-1", "0", "1"), size = 4, colour = "grey30") +
    annotate("text", x = c(-1.1, 1.1), y = 0.45, label = c("sensorimotor", "association"),
             size = 4.5, fontface = "italic", colour = "grey25")
  for (i in seq_len(nrow(in_df_items))) {
    x_true = in_df_items$r[i]
    x_disp = in_df_items$x_display[i]
    p = p +
      annotate("segment", x = c(x_true, x_true), xend = c(x_true, x_disp), y = c(0.15, 0.5),
               yend = c(0.5, 0.85), linewidth = 0.4, colour = "grey50") +
      annotate("point", x = x_true, y = 0, size = 3, shape = 21, fill = "white", colour = "grey30", stroke = 0.8) +
      annotate("text", x = x_disp, y = 1.0, label = in_df_items$label[i], size = 3.8, fontface = "bold") +
      annotate("text", x = x_disp, y = 1.6, label = in_df_items$label_text[i], size = 3.5,
               colour = "grey20", lineheight = 0.85) +
      annotation_custom(list_grobs[[i]], xmin = x_disp - 0.28, xmax = x_disp + 0.28, ymin = 2.3, ymax = 4.6)
  }
  p = p + coord_cartesian(xlim = c(-1.3, 1.3), ylim = c(y_bottom, y_top), expand = FALSE, clip = "off") +
    theme(plot.title = element_text(face = "bold", size = 15, hjust = 0, margin = margin(b = 8)),
          plot.margin = margin(4, 4, 4, 4, "mm"))
  return(p)
}

## Colour bar legends
fColorbar = function(color_fn, cbar_min, cbar_max, title_text) {
  df_dummy = data.frame(x = 1:100, val = seq(cbar_min, cbar_max, length.out = 100))
  p = ggplot(df_dummy, aes(x = x, y = 1, fill = val)) +
    geom_tile() +
    scale_fill_gradientn(colours = color_fn(100), limits = c(cbar_min, cbar_max),
                         breaks = round(seq(cbar_min, cbar_max, length.out = 3), 2)) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11)) +
    guides(fill = guide_colourbar(title = title_text, title.position = "top", title.hjust = 0.5,
                                  barwidth = 10, barheight = 0.8))
  return(get_legend(p))
}

## Colour ranges shared by CT and SA (SNP maps excluded: their values are on a different scale)
cols_no_SNPs = c("CT_CNVs", "SA_CNVs", "CT_NPDs", "SA_NPDs")
list_cbar = list(meanAbsES = c(0, ceiling(max(df_concensus_profiles_meanAbsES[, cols_no_SNPs]) * 100) / 100),
                 significance = c(0, ceiling(max(df_concensus_profiles_significance) / 5) * 5),
                 SD = c(0, ceiling(max(df_concensus_profiles_variance[, cols_no_SNPs]) * 10) / 10))
list_colors = list(meanAbsES = meanES_colors, significance = sig_colors, SD = sd_colors)
list_titles = c(meanAbsES = "A. Mean absolute effect sizes", significance = "B. % Significance maps",
                SD = "C. Standard deviation maps")
label_text = c(CNVs = "Rare genetic\nvariants (CNVs)", NPDs = "NPD\ndiagnosis", SNPs = "Common genetic\nvariants (SNPs)")

## Plot Supplementary Figure 14 A to C
list_panels = list()
for (measure in names(list_consensus)) {
  df_consensus = list_consensus[[measure]]
  in_df_items = df_stats[df_stats$measure == measure, ]
  in_df_items$label_text = label_text[sub(paste0("^", hemi, "_"), "", in_df_items$map)]
  in_df_items$x_display = fDodge(in_df_items$r)
  cbar = list_cbar[[measure]]
  list_grobs = list()
  for (i in seq_len(nrow(in_df_items))) {
    map_vals = df_consensus[, in_df_items$map[i]]
    ## SNP maps (mean abs Beta and SD of Beta) rescaled to the CNV/NPD range for display only
    if (grepl("SNPs", in_df_items$map[i]) && measure != "significance") {
      vals_ref = unlist(df_consensus[, setdiff(grep(paste0("^", hemi, "_"), colnames(df_consensus), value = TRUE), in_df_items$map[i])])
      map_vals = min(vals_ref) + (map_vals - min(map_vals)) / (max(map_vals) - min(map_vals)) * (max(vals_ref) - min(vals_ref))
    }
    list_grobs[[i]] = fBrainPlot_LH(map_vals, list_colors[[measure]], cbar[1], cbar[2])
  }
  list_panels[[measure]] = fGradientAxisPanel(in_df_items, list_grobs, y_bottom = -0.9, y_top = 4.8) +
    ggtitle(list_titles[measure])
}

save_panel(list_panels[["meanAbsES"]], paste0(fig_name, "A"), 12, 3.2)
save_panel(list_panels[["significance"]], paste0(fig_name, "B"), 12, 3.2)
save_panel(list_panels[["SD"]], paste0(fig_name, "C"), 12, 3.2)

## Plot Supplementary Figure 14: composite
title_grob = textGrob(paste0("Consensus maps — Correlation with cortical gradient: ", hemi_label),
                      gp = gpar(fontsize = 18, fontface = "bold"))
axis_title = textGrob("Correlation with cortical gradient (r)", gp = gpar(fontsize = 16, fontface = "bold.italic", col = "grey15"))
legend_row = plot_grid(fColorbar(meanES_colors, list_cbar$meanAbsES[1], list_cbar$meanAbsES[2], "Mean abs. effect size"),
                       fColorbar(sig_colors, list_cbar$significance[1], list_cbar$significance[2], "% Significant"),
                       textGrob("* p-spin < 0.05", gp = gpar(fontsize = 13, fontface = "italic")),
                       fColorbar(sd_colors, list_cbar$SD[1], list_cbar$SD[2], "SD of effect sizes"),
                       ncol = 4, rel_widths = c(1, 1, 0.6, 1))
p_composite = plot_grid(title_grob, list_panels[["meanAbsES"]], list_panels[["significance"]], list_panels[["SD"]],
                        axis_title, legend_row, ncol = 1, rel_heights = c(0.06, 1, 1, 1, 0.1, 0.18))
save_panel(p_composite, fig_name, 12, 12)

## Save source data
for (measure in names(list_consensus)) {
  panel_letter = c(meanAbsES = "A", significance = "B", SD = "C")[measure]
  df_consensus = list_consensus[[measure]]
  map_names = df_stats$map[df_stats$measure == measure]
  df_regional = data.frame(region = rep(rois_desikan, length(map_names)),
                           map = rep(map_names, each = 34),
                           value = unlist(df_consensus[, map_names], use.names = FALSE),
                           cortical_gradient = rep(gradient_map, length(map_names)))
  df_regional = merge(df_regional, df_stats[df_stats$measure == measure, c("map", "r", "p_spin", "label")], by = "map", sort = FALSE)
  save_source_data(df_regional, paste0(fig_name, panel_letter))
}
