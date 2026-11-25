#-------------------------------------------------------------------------------
# Supplementary Figure 17: Correlation between latent dimensions and other cortical reference maps
#-------------------------------------------------------------------------------
# Legend: Pearson correlations between the first principal component (PC1) maps of regional effect sizes
# (columns: SA for NPD-associated SNPs; CT and SA for CNVs; CT and SA for NPDs) and 20 cortical reference
# maps (rows). The reference maps are organizational, developmental and structural maps (allometric
# scaling, intersubject functional variability, functional gradient, evolutionary and developmental
# expansion, cortical thickness, T1w/T2w ratio, PC1 of gene expression, PC1 of NeuroSynth meta-analytic
# terms); metabolic and perfusion maps (glucose and oxygen metabolism, cerebral blood volume and flow);
# and MEG power maps (delta to high-gamma bands). All maps were obtained from neuromaps and parcellated
# into 34 left-hemisphere Desikan-Killiany regions. Colour and values indicate Pearson r. *: p-spin<0.05
# (one-sided, 10,000 spin permutations, not corrected for multiple comparisons). Correlations and
# p-values are provided in Supplementary Table 10.
# Abbreviations: CNV=copy number variant; CT=cortical thickness; NPD=neurodevelopmental and psychiatric
# disorders; PC=principal component; SA=surface area; SNP=single nucleotide polymorphism.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")

## Load data
## PC1 latent profiles
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))
df_PC1 = df_concensus_profiles_latentDim_PC1

## Reference maps (neuromaps, Desikan atlas); left hemisphere rows in ggseg order
df_neuromaps = read.csv(paste0(data_dir, "df_neuromaps_20maps_Desikan_68rois.csv"), row.names = 1)
df_neuromaps_info = read.csv(paste0(data_dir, "df_all_neuromaps_info_20maps.csv"), row.names = 1)
df_neuromaps_LH = df_neuromaps[paste0("lh_", rois_desikan), ]

## PC1 of gene expression is the cortical gradient of the main analyses (abagen x -1: sensorimotor to
## association)
df_neuromaps_LH[, "abagen"] = -df_neuromaps_LH[, "abagen"]
df_neuromaps_info$map_annotation = df_neuromaps_info$map_colname
df_neuromaps_info$map_annotation[df_neuromaps_info$map_colname == "abagen"] = "abagen (x -1; main gradient)"

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

array_ref_maps = rev(df_neuromaps_info$map_colname)
array_PC1_maps = c("SA_SNPs", "CT_CNVs", "SA_CNVs", "CT_NPDs", "SA_NPDs")

## Stats: Pearson correlation and spin permutation p-value (20 reference maps x 5 PC1 maps)
df_stats = data.frame()
for (ref_map in array_ref_maps) {
  for (map_name in array_PC1_maps) {
    df_stats = rbind(df_stats, data.frame(reference_map = df_neuromaps_info$map_annotation[df_neuromaps_info$map_colname == ref_map],
                                          reference_label = df_neuromaps_info$map_plot_label[df_neuromaps_info$map_colname == ref_map],
                                          map = map_name,
                                          r = cor(df_PC1[, map_name], df_neuromaps_LH[, ref_map]),
                                          p_spin = perm.sphere.p(df_PC1[, map_name], df_neuromaps_LH[, ref_map], df_perm_ids_LH)))
  }
}
print(df_stats)

## Plot Supplementary Figure 17: heatmap of correlations
in_df_plot = df_stats
in_df_plot$reference_label = factor(in_df_plot$reference_label, levels = rev(unique(in_df_plot$reference_label)))
in_df_plot$map = factor(in_df_plot$map, levels = array_PC1_maps)
in_df_plot$label = paste0(sprintf("%.2f", in_df_plot$r), ifelse(in_df_plot$p_spin < 0.05, "*", ""))

p_SupFig17 = ggplot(in_df_plot, aes(x = map, y = reference_label, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "#3B7BBF", mid = "white", high = "#E8392B", midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(x = "PC1 map", y = "Reference map (neuromaps, Desikan LH 34)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        legend.key.height = unit(1.6, "cm"))
save_panel(p_SupFig17, "SupFig17", 7.5, 8.5)

## Save source data
save_source_data(df_stats, "SupFig17")
