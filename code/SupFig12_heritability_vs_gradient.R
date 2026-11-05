#-------------------------------------------------------------------------------
# Supplementary Figure 12: Correlations between cortical gradient and regional twin and SNP heritability estimates
#-------------------------------------------------------------------------------
# Legend: Regional profiles of twin and SNP heritability estimates for cortical thickness and surface
# area across 34 Desikan cortical regions from Grasby et al. 10, plotted against cortical gradient
# (sensorimotor-to-association regions). Each point represents heritability estimates. The surface area
# estimates are adjusted (adj) for Total SA. X-axis: cortical regions ordered according to the cortical
# gradient from sensorimotor to association regions. Dotted line: correlation with the cortical
# gradient. Each panel displays Pearson correlation and *:spin-permutation significant, p-spin < 0.05.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)

## Load data
## Cortical gradient (-1 x PC1 of AHBA gene expression, abagen)
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))
gradient_map = array_cortical_gradient

## Regional twin and SNP heritability estimates (Grasby et al. 2020)
load(paste0(data_dir, "df_Grasby_h2_regional_RData_file_with_95CI_apr2025.RData"))

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

df_h2_info = data.frame(h2_column = c("CT_wo_adj_Twin_h2", "SA_adj_Twin_h2", "CT_wo_adj_SNP_h2", "SA_adj_SNP_h2"),
                        y_label = c("CT twin h2", "SA twin h2 (adj)", "CT snp h2", "SA snp h2 (adj)"),
                        metric = c("CT", "SA", "CT", "SA"),
                        ylim_min = c(0.28, 0.17, 0.06, 0.06),
                        ylim_max = c(0.72, 0.77, 0.28, 0.33))

## Stats: Pearson correlation with the gradient and spin permutation p-value
df_h2_info$r = NA
df_h2_info$p_spin = NA
for (i in 1:nrow(df_h2_info)) {
  h2_vals = df_Grasby_h2_regional[, df_h2_info$h2_column[i]]
  df_h2_info$r[i] = cor(h2_vals, gradient_map)
  df_h2_info$p_spin[i] = perm.sphere.p(h2_vals, gradient_map, df_perm_ids_LH)
}
df_h2_info$label = paste0("r= ", round(df_h2_info$r, 2), ifelse(df_h2_info$p_spin < 0.05, "*", ""))
print(df_h2_info[, c("h2_column", "r", "p_spin", "label")])

## Plot: regions ordered by the gradient (sensorimotor to association)
input_order = rois_desikan[order(gradient_map)]

fPoint_plot_h2 = function(i) {
  in_df_plot = data.frame(ROI = factor(rois_desikan, levels = input_order),
                          h2 = df_Grasby_h2_regional[, df_h2_info$h2_column[i]])
  in_df_plot$ROI_numeric = as.numeric(in_df_plot$ROI)
  p = ggplot(in_df_plot, aes(x = ROI, y = h2)) +
    theme_prism(base_size = 18, base_line_size = 0.75) +
    geom_point(size = 4, colour = col_values_metric[match(df_h2_info$metric[i], c("CT", "SA"))]) +
    geom_smooth(aes(x = ROI_numeric, y = h2), method = "lm", formula = y ~ x, colour = "black", linetype = "dashed", se = FALSE) +
    annotate("text", x = 34.5, y = df_h2_info$ylim_max[i], label = df_h2_info$label[i], hjust = 1, vjust = 1,
             size = 6.5, fontface = "bold") +
    ylim(c(df_h2_info$ylim_min[i], df_h2_info$ylim_max[i])) +
    labs(x = NULL, y = df_h2_info$y_label[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(p)
}

list_panels = lapply(1:nrow(df_h2_info), fPoint_plot_h2)
for (i in 1:3) list_panels[[i]] = list_panels[[i]] + theme(axis.text.x = element_blank())

p_composite = plot_grid(plotlist = list_panels, ncol = 1, align = "v", rel_heights = c(0.21, 0.21, 0.21, 0.37))
save_panel(p_composite, "SupFig12", 15, 13)

## Save source data
df_source = data.frame()
for (i in 1:nrow(df_h2_info)) {
  df_source = rbind(df_source, data.frame(h2_measure = df_h2_info$y_label[i],
                                          region = rois_desikan,
                                          gradient_rank = match(rois_desikan, input_order),
                                          cortical_gradient = gradient_map,
                                          h2 = df_Grasby_h2_regional[, df_h2_info$h2_column[i]],
                                          r = df_h2_info$r[i],
                                          p_spin = df_h2_info$p_spin[i]))
}
save_source_data(df_source, "SupFig12")
