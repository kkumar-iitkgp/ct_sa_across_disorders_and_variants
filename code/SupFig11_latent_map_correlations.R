#-------------------------------------------------------------------------------
# Supplementary Figure 11: Correlations between consensus maps
#-------------------------------------------------------------------------------
# Legend: Correlation between pairs of consensus maps for A) mean absolute effect sizes, B) %
# significance, C) standard deviation of effect sizes, and D) principal component 1 (PC1). *: spin
# permutation p-value < 0.05, which preserves the spatial structure of the cortex while generating null
# distributions to account for spatial autocorrelation.
# Abbreviations. CNV=copy number variants; Cor.=correlation; CT: cortical thickness; SA=surface area;
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(grid)

## Load data
## Consensus maps: mean absolute effect size, % significance, standard deviation (stored as
## "variance") and PC1, for SNPs SA, CNVs CT, CNVs SA, NPDs CT and NPDs SA
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

list_consensus = list(A = df_concensus_profiles_meanAbsES,
                      B = df_concensus_profiles_significance,
                      C = df_concensus_profiles_variance,
                      D = df_concensus_profiles_latentDim_PC1)
list_titles = c(A = "A. Mean absolute effect size per ROI",
                B = "B. % significance per ROI",
                C = "C. Standard deviation of effect sizes per ROI",
                D = "D. Latent dimension (PC1)")
array_map_names = c("SNPs SA", "CNVs CT", "CNVs SA", "NPDs CT", "NPDs SA")

## Stats: pair-wise Pearson correlation and spin permutation p-value (lower triangle)
list_stats = list()
for (panel in names(list_consensus)) {
  df_maps = list_consensus[[panel]]
  df_stats = data.frame()
  for (loop_col in 1:(ncol(df_maps) - 1)) {
    for (loop_row in (loop_col + 1):ncol(df_maps)) {
      df_stats = rbind(df_stats, data.frame(map_row = array_map_names[loop_row],
                                            map_col = array_map_names[loop_col],
                                            r = cor(df_maps[, loop_row], df_maps[, loop_col]),
                                            p_spin = perm.sphere.p(df_maps[, loop_row], df_maps[, loop_col], df_perm_ids_LH)))
    }
  }
  df_stats$star = ifelse(df_stats$p_spin < 0.05, "*", "")
  list_stats[[panel]] = df_stats
  print(list_titles[panel])
  print(df_stats)
}

## Plot: lower-triangle correlation heatmap
fCorHeatmap_lower = function(df_stats) {
  mat_cor = diag(length(array_map_names))
  mat_pspin = matrix(1, length(array_map_names), length(array_map_names))
  rownames(mat_cor) = colnames(mat_cor) = array_map_names
  for (i in seq_len(nrow(df_stats))) {
    mat_cor[df_stats$map_row[i], df_stats$map_col[i]] = df_stats$r[i]
    mat_pspin[match(df_stats$map_row[i], array_map_names), match(df_stats$map_col[i], array_map_names)] = df_stats$p_spin[i]
  }
  col_fun = colorRamp2(c(-1, 0, 1), hcl_palette = "RdYlBu", reverse = TRUE)
  ht = Heatmap(mat_cor, name = "Cor", rect_gp = gpar(type = "none"), col = col_fun,
               cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left",
               column_names_rot = 60, row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 16),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (i >= j) {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = fill, fill = fill))
                   if (i > j) {
                     is_sig = mat_pspin[i, j] < 0.05
                     grid.text(paste0(round(mat_cor[i, j], 2), "\n", ifelse(is_sig, "*", " ")), x, y,
                               gp = gpar(fontsize = ifelse(is_sig, 17, 16), fontface = ifelse(is_sig, "bold", "plain")))
                   }
                 }
               },
               heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1), labels = c("-1", "0", "1"),
                                           legend_height = unit(4, "cm"), grid_width = unit(0.6, "cm")))
  return(grid.grabExpr(draw(ht)))
}

list_panels = list()
for (panel in names(list_consensus)) {
  p_title = ggdraw() + draw_label(list_titles[panel], x = 0.02, hjust = 0, size = 15)
  list_panels[[panel]] = plot_grid(p_title, fCorHeatmap_lower(list_stats[[panel]]), ncol = 1, rel_heights = c(0.08, 1))
  save_panel(list_panels[[panel]], paste0("SupFig11", panel), 6, 5.5)
}

## Plot Supplementary Figure 11: composite
p_composite = plot_grid(plotlist = list_panels, ncol = 2)
save_panel(p_composite, "SupFig11", 12, 11)

## Save source data
for (panel in names(list_consensus)) {
  save_source_data(list_stats[[panel]], paste0("SupFig11", panel))
}
