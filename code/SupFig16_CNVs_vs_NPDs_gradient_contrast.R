#-------------------------------------------------------------------------------
# Supplementary Figure 16: Difference between genetic variants and NPDs in their correlation with the cortical gradient
#-------------------------------------------------------------------------------
# Legend: Correlation (Pearson's r) between the cortical gradient and regional maps of A) rare genetic
# variants (CNVs) and NPDs for cortical thickness (CT), B) CNVs and NPDs for surface area (SA), and
# C) NPD-associated common variants (SNPs) and NPDs for SA. Rows show the latent dimension (PC1) and the
# three consensus maps (mean absolute effect size, standard deviation, % significance). The line shows
# the difference between the two correlations (delta r). p-spin: spin permutation p-value of this
# difference (one-sided): the cortical gradient was rotated 10,000 times and the same rotation was
# applied to both correlations. Values, and the results with the updated MDD lead SNPs, are given in
# Supplementary Table 9.
# Abbreviations: CNV=copy number variant; CT=cortical thickness; NPD=neurodevelopmental and
# psychiatric disorders; PC: principal component; SA=surface area; SNP: single nucleotide polymorphism.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")

## Load data
## Consensus maps: mean absolute effect size, % significance, standard deviation (stored as
## "variance") and PC1, for SNPs SA, CNVs CT, CNVs SA, NPDs CT and NPDs SA; cortical gradient
load(paste0(data_dir, "data_fig5_concensus_profiles_meanAbsES_variance_significance_latentDimPC1.RData"))

## Pre-computed spin rotations (10,000) for the 34 left-hemisphere Desikan regions
load(paste0(data_dir, "df_perm_ids_DesikanLH34_nIterNull_10000.RData"))

list_consensus = list("PC1" = df_concensus_profiles_latentDim_PC1,
                      "Mean abs. effect size" = df_concensus_profiles_meanAbsES,
                      "Standard deviation" = df_concensus_profiles_variance,
                      "% significance" = df_concensus_profiles_significance)

df_pairs = data.frame(panel = c("A", "B", "C"),
                      comparison = c("CT: CNVs vs NPDs", "SA: CNVs vs NPDs", "SA: SNPs vs NPDs"),
                      map_1 = c("CT_CNVs", "SA_CNVs", "SA_SNPs"),
                      map_2 = c("CT_NPDs", "SA_NPDs", "SA_NPDs"))

## Stats Supplementary Figure 16: difference between two correlations with the cortical gradient,
## spin permutation p-value (one-sided, map 1 more sensorimotor than map 2; same rotation for both)
perm.sphere.p.diff = function(x1, x2, y, perm.id) {
  diff.emp = cor(x1, y) - cor(x2, y)
  diff.null = apply(perm.id, 2, function(i) cor(x1, y[i]) - cor(x2, y[i]))
  return(mean(diff.null <= diff.emp))
}

df_stats = data.frame()
for (map_type in names(list_consensus)) {
  for (i in seq_len(nrow(df_pairs))) {
    array_map_1 = list_consensus[[map_type]][, df_pairs$map_1[i]]
    array_map_2 = list_consensus[[map_type]][, df_pairs$map_2[i]]
    df_stats = rbind(df_stats, data.frame(panel = df_pairs$panel[i], comparison = df_pairs$comparison[i],
                                          consensus_measure = map_type,
                                          map_1 = df_pairs$map_1[i], map_2 = df_pairs$map_2[i],
                                          r_map_1 = cor(array_map_1, array_cortical_gradient),
                                          r_map_2 = cor(array_map_2, array_cortical_gradient),
                                          delta_r = cor(array_map_1, array_cortical_gradient) - cor(array_map_2, array_cortical_gradient),
                                          p_spin = perm.sphere.p.diff(array_map_1, array_map_2, array_cortical_gradient, df_perm_ids_LH)))
  }
}
df_stats = df_stats[order(df_stats$panel), ]
print(df_stats)

## Plot Supplementary Figure 16
col_values_group = c("CNVs" = "#E7298A", "SNPs" = "#377eb8", "NPDs" = "#1B9E77")

in_df_label = df_stats
in_df_label$facet = paste0(in_df_label$panel, ". ", in_df_label$comparison)
in_df_label$consensus_measure = factor(in_df_label$consensus_measure, levels = rev(names(list_consensus)))
in_df_label$label = sprintf("\u0394r = %.2f\np-spin = %s", in_df_label$delta_r, formatC(in_df_label$p_spin, digits = 2, format = "g"))

in_df_plot = rbind(data.frame(in_df_label[, c("facet", "consensus_measure")], group = sub("^.*_", "", in_df_label$map_1), r = in_df_label$r_map_1),
                   data.frame(in_df_label[, c("facet", "consensus_measure")], group = sub("^.*_", "", in_df_label$map_2), r = in_df_label$r_map_2))
in_df_plot$group = factor(in_df_plot$group, levels = names(col_values_group))

p_SupFig16 = ggplot(in_df_plot, aes(x = r, y = consensus_measure)) +
  theme_prism(base_size = input_base_font_size, base_line_size = 0.75) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = in_df_label, aes(x = r_map_1, xend = r_map_2, y = consensus_measure, yend = consensus_measure),
               colour = "grey40", linewidth = 0.8) +
  geom_point(aes(colour = group), size = 4) +
  geom_text(data = in_df_label, aes(x = 1.02, y = consensus_measure, label = label), hjust = 0, size = 3.5, lineheight = 0.9) +
  scale_colour_manual(values = col_values_group, name = NULL) +
  scale_x_continuous(limits = c(-1, 1.8), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  facet_wrap(~ facet, nrow = 1) +
  xlab("Correlation with cortical gradient (r)\n\u2190 sensorimotor      association \u2192") + ylab(NULL) +
  theme(legend.position = "top", strip.text = element_text(face = "bold", size = input_base_font_size, hjust = 0))
save_panel(p_SupFig16, "SupFig16", 13, 4.6)

## Save source data
save_source_data(df_stats, "SupFig16")
