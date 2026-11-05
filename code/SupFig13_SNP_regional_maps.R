#-------------------------------------------------------------------------------
# Supplementary Figure 13: Regional cortical alterations across NPD-associated common variants (SNPs)
#-------------------------------------------------------------------------------
# Legend: A) Surface area (SA) cortical GWAS 10 regression estimate maps (Z = Beta/SE) for 20
# NPD-associated genetic variants (SNPs). FDR (q<0.05) significant regions are shown in black
# boundaries. B-C) Principal component analysis across SNP regional profiles. B) PC1 variable loadings
# for SA. C) The variance explained by the first two principal components (PC1 and PC2).
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)
library(grid)

## Load data
## Regional SA estimates (Beta, Z = Beta/SE, FDR-adjusted p) of the 20 SNPs, Grasby et al. 2020 GWAS;
## columns are regions in alphabetical order
load(paste0(data_dir, "CrossDisorder_stacked_FDRsnps_df_beta_se_z_pval_SNP_x_ROI_surfavg_ENIGMA3_withGlobal_5feb2025.RData"))

## PC1 loadings of the SNP SA PCA (as submitted)
load(paste0(data_dir, "res_df_summnary_nsig_results_gradient_regional_SA_SNPs_FDR_5feb2025.RData"))

snp_names = rownames(col_loadings_SA_SNPs)
rownames(df_z_cross_dis_SA) = snp_names
df_pval_FDR = as.data.frame(df_pval_cross_dis_SA_FDR)
dimnames(df_pval_FDR) = list(snp_names, colnames(df_z_cross_dis_SA))

## Regions in ggseg order
df_z_SA = df_z_cross_dis_SA[, rois_desikan]
df_pval_FDR = df_pval_FDR[, rois_desikan]

## Stats Supplementary Figure 13B and 13C: PCA across the 20 SNP regional Beta profiles
## (34 regions x 20 SNPs, variables scaled)
res_pca = prcomp(t(df_beta_cross_dis_SA), scale. = TRUE)
array_var_exp = 100 * res_pca$sdev^2 / sum(res_pca$sdev^2)
df_var_exp = data.frame(PC = c("PC1", "PC2"), var_exp = array_var_exp[1:2], var_exp_rounded = round(array_var_exp[1:2]))
print(df_var_exp)

df_loadings = data.frame(SNP = snp_names, PC1_loading = col_loadings_SA_SNPs[, 1])
print(df_loadings)

## Plot Supplementary Figure 13A: Z maps with FDR-significant regions outlined
coolwarm = colorRampPalette(c("#3B4CC0", "#6B8DF0", "#A9C5FC", "#D6DFE8", "#F7F7F7",
                              "#F5CBA7", "#F0946C", "#D7191C", "#B40426"))
z_cbar = ceiling(min(max(abs(df_z_SA)), 5))

fBrainPlot_Z_FDR = function(array_z, array_pval) {
  array_z = pmin(pmax(array_z, -z_cbar), z_cbar)
  array_pval = c(array_pval[1:3], 1, array_pval[4:34])
  df_brain = data.frame(region = names(dk$palette),
                        value = c(array_z[1:3], NA, array_z[4:34]),
                        sig_col = ifelse(array_pval > 0.05, "grey60", "black"),
                        sig_thick = ifelse(array_pval > 0.05, 0.05, pmin(-0.2 * log10(array_pval), 0.3)))
  p = ggseg(df_brain, atlas = "dk", position = "dispersed", hemisphere = "left",
            mapping = aes(fill = value, colour = I(sig_col), size = I(sig_thick))) +
    scale_fill_gradientn(limits = c(-z_cbar, z_cbar), colours = coolwarm(100)) +
    theme_brain(text.family = "sans", text.size = 14) +
    theme(legend.position = "bottom", legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14)) +
    guides(fill = guide_colourbar(title = "Z (Beta/SE)", title.position = "top", title.hjust = 0.5,
                                  barwidth = 12, barheight = 1))
  return(p)
}

list_maps = list()
for (i in seq_along(snp_names)) {
  p_map = fBrainPlot_Z_FDR(as.numeric(df_z_SA[i, ]), as.numeric(df_pval_FDR[i, ]))
  if (i == 1) legend_A = get_legend(p_map)
  list_maps[[i]] = p_map + ggtitle(snp_names[i]) +
    theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(0, 0, 0, 0)),
          plot.margin = margin(0, 0, 0, 0))
}
p_SupFig13A = plot_grid(plotlist = list_maps, ncol = 1)
p_SupFig13A = plot_grid(ggdraw() + draw_label("A. SA effect size maps for NPD SNPs", fontface = "bold", size = 14, x = 0.02, hjust = 0),
                        p_SupFig13A, legend_A, ncol = 1, rel_heights = c(0.02, 1, 0.05))
save_panel(p_SupFig13A, "SupFig13A", 5, 18)

## Plot Supplementary Figure 13B: PC1 loadings
pca_loadings_colors = colorRampPalette(c("#5e3c99", "#b2abd2", "#f7f7f7", "#fdb863", "#e66101"))
max_abs_load = max(abs(df_loadings$PC1_loading))
in_df_plot = df_loadings
in_df_plot$SNP = factor(in_df_plot$SNP, levels = rev(snp_names))

p_SupFig13B = ggplot(in_df_plot, aes(x = PC1_loading, y = SNP, fill = PC1_loading)) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
  scale_fill_gradientn(colours = pca_loadings_colors(100), limits = c(-max_abs_load, max_abs_load),
                       breaks = c(-round(max_abs_load, 1), 0, round(max_abs_load, 1))) +
  labs(x = "Loadings", y = NULL, title = "B. SA PC1 loadings") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.3), axis.text.y = element_text(size = 11, face = "bold"),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom", legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)) +
  guides(fill = guide_colourbar(title = "Loadings", title.position = "top", title.hjust = 0.5,
                                barwidth = 12, barheight = 1))
save_panel(p_SupFig13B, "SupFig13B", 7, 14)

## Plot Supplementary Figure 13C: variance explained by PC1 and PC2
in_df_plot = df_var_exp
in_df_plot$PC = factor(in_df_plot$PC, levels = c("PC2", "PC1"))

p_SupFig13C = ggplot(in_df_plot, aes(x = var_exp_rounded, y = PC)) +
  geom_col(fill = "#2B3A67", width = 0.4) +
  geom_text(aes(label = paste0(var_exp_rounded, "%")), hjust = -0.1, size = 5) +
  scale_x_continuous(limits = c(0, max(in_df_plot$var_exp_rounded) * 1.3), expand = c(0, 0)) +
  labs(x = "% Variance explained", y = NULL, title = "C. SA % var. explained") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14))
save_panel(p_SupFig13C, "SupFig13C", 7, 3.5)

## Plot Supplementary Figure 13: composite
p_right = plot_grid(p_SupFig13B, p_SupFig13C, ncol = 1, rel_heights = c(0.8, 0.2))
p_composite = plot_grid(p_SupFig13A, p_right, ncol = 2, rel_widths = c(1, 1.1))
save_panel(p_composite, "SupFig13", 12, 18)

## Save source data
df_source_A = data.frame(SNP = rep(snp_names, each = 34),
                         region = rep(rois_desikan, length(snp_names)),
                         Z = as.vector(t(as.matrix(df_z_SA))),
                         p_FDR = as.vector(t(as.matrix(df_pval_FDR))))
df_source_A$FDR_significant = df_source_A$p_FDR < 0.05
save_source_data(df_source_A, "SupFig13A")
save_source_data(df_loadings, "SupFig13B")
save_source_data(df_var_exp, "SupFig13C")
