#-------------------------------------------------------------------------------
# Supplementary Figure 16: Correlations between latent dimensions (PC1) and other cortical reference maps
#-------------------------------------------------------------------------------
# Legend: Correlation (Pearson's r) between the latent dimension (PC1) maps of rare genetic variants
# (CNVs), NPD diagnosis (NPDs) and common genetic variants (SNPs) for cortical thickness (CT) and
# surface area (SA), and two reference maps from neuromaps: A) the principal functional connectivity
# gradient (Margulies et al. 2016) and B) cortical evolutionary expansion (Hill et al. 2010), across 34
# left-hemisphere Desikan cortical regions. Filled circles and *: spin-permutation significant, p-spin
# < 0.05.
# Abbreviations: CNV=copy number variant; CT=cortical thickness; NPD=neurodevelopmental and psychiatric
# disorders; PC=principal component; SA=surface area; SNP=single nucleotide polymorphism.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)

## Load data
## PC1 latent profiles
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))
df_PC1 = df_concensus_profiles_latentDim_PC1

## Reference maps (neuromaps, Desikan atlas); left hemisphere rows in ggseg order
df_neuromaps = read.csv(paste0(data_dir, "df_neuromaps_20maps_Desikan_68rois.csv"), row.names = 1)
df_neuromaps_info = read.csv(paste0(data_dir, "df_all_neuromaps_info_20maps.csv"), row.names = 1)
df_neuromaps_LH = df_neuromaps[paste0("lh_", rois_desikan), ]

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

array_ref_maps = c("fcgradient01", "evoexp")
array_PC1_maps = c("CT_CNVs", "SA_CNVs", "CT_NPDs", "SA_NPDs", "SA_SNPs")

## Stats: Pearson correlation and spin permutation p-value
df_stats = data.frame()
for (ref_map in array_ref_maps) {
  for (map_name in array_PC1_maps) {
    df_stats = rbind(df_stats, data.frame(reference_map = ref_map,
                                          reference_label = df_neuromaps_info$map_plot_label[df_neuromaps_info$map_colname == ref_map],
                                          map = map_name,
                                          r = cor(df_PC1[, map_name], df_neuromaps_LH[, ref_map]),
                                          p_spin = perm.sphere.p(df_PC1[, map_name], df_neuromaps_LH[, ref_map], df_perm_ids_LH)))
  }
}
df_stats$sig_status = ifelse(df_stats$p_spin < 0.05, "Sig", "ns")
df_stats$label = paste0("r = ", round(df_stats$r, 2), ifelse(df_stats$p_spin < 0.05, "*", ""))
print(df_stats)

## Plot Supplementary Figure 16A and 16B
map_labels = c(CT_CNVs = "PC1 CNVs CT", SA_CNVs = "PC1 CNVs SA", CT_NPDs = "PC1 NPDs CT",
               SA_NPDs = "PC1 NPDs SA", SA_SNPs = "PC1 SNPs SA")
df_stats$map_label = factor(map_labels[df_stats$map], levels = rev(map_labels))
df_stats$metric = substr(df_stats$map, 1, 2)

fRefMapPanel = function(ref_map, title_text) {
  in_df_plot = df_stats[df_stats$reference_map == ref_map, ]
  p = ggplot(in_df_plot, aes(x = r, y = map_label)) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", colour = "grey80") +
    geom_segment(aes(x = 0, xend = r, yend = map_label, colour = metric), linewidth = 0.8) +
    geom_point(aes(colour = metric, fill = interaction(metric, sig_status)), shape = 21, size = 5, stroke = 1) +
    geom_text(aes(label = label, x = r + ifelse(r > 0, 0.08, -0.08), hjust = ifelse(r > 0, 0, 1)), size = 5) +
    scale_colour_manual(values = c(CT = col_CT, SA = col_SA)) +
    scale_fill_manual(values = c(CT.Sig = col_CT, SA.Sig = col_SA, CT.ns = "white", SA.ns = "white")) +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    labs(x = "Correlation (r)", y = NULL, title = title_text) +
    theme_prism(base_size = 16, base_line_size = 0.75) +
    theme(legend.position = "none")
  return(p)
}

p_SupFig16A = fRefMapPanel("fcgradient01", "A. Functional connectivity gradient")
p_SupFig16B = fRefMapPanel("evoexp", "B. Evolutionary expansion")
save_panel(p_SupFig16A, "SupFig16A", 7, 4.5)
save_panel(p_SupFig16B, "SupFig16B", 7, 4.5)

## Plot Supplementary Figure 16: composite
p_composite = plot_grid(p_SupFig16A, p_SupFig16B + theme(axis.text.y = element_blank()), ncol = 2, rel_widths = c(1.25, 1))
p_composite = plot_grid(p_composite, ggdraw() + draw_label("* p-spin < 0.05", fontface = "italic", size = 13),
                        ncol = 1, rel_heights = c(1, 0.07))
save_panel(p_composite, "SupFig16", 13, 5)

## Save source data
df_source = df_stats[, c("reference_map", "reference_label", "map", "r", "p_spin")]
save_source_data(df_source[df_source$reference_map == "fcgradient01", ], "SupFig16A")
save_source_data(df_source[df_source$reference_map == "evoexp", ], "SupFig16B")
