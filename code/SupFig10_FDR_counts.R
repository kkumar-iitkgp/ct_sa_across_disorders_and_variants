#-------------------------------------------------------------------------------
# Supplementary Figure 10: Number of FDR significant regions per cortical map
#-------------------------------------------------------------------------------
# Legend: Number of FDR significant cortical regions (out of 34) per CNV/NPD for cortical thickness
# (up, red) and surface area (down, blue) for A) CNVs; B) NPDs; C) NPD medication; and D)
# NPD-unmedicated.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(cowplot)

## Load data
load(paste0(data_dir, "df_pval_combined_and_sep_CT_CNVs_NPDs_med_umed_5feb2025.RData"))
load(paste0(data_dir, "df_pval_combined_and_sep_SA_CNVs_NPDs_med_umed_5feb2025.RData"))


## Stats Supplementary Figure 10: number of FDR significant regions per condition
## CNVs: regional p-values as provided in the input
## NPDs, NPD-med, NPD-umed: FDR across all regions and conditions, CT and SA together, within each block
fdr_block = function(df_pval_CT, df_pval_SA) {
  p_stacked = cbind(df_pval_CT, df_pval_SA)
  p_stacked_adj = matrix(p.adjust(as.vector(as.matrix(p_stacked)), method = "fdr"), ncol = ncol(p_stacked))
  list(CT = p_stacked_adj[, 1:ncol(df_pval_CT)], SA = p_stacked_adj[, (ncol(df_pval_CT) + 1):ncol(p_stacked)])
}

list_pval_adj = list("CNV" = list(CT = df_pval_CT_CNVs, SA = df_pval_SA_CNVs),
                     "NPD" = fdr_block(df_pval_CT_NPDs_all, df_pval_SA_NPDs_all),
                     "NPD-med" = fdr_block(df_pval_CT_NPDs_med, df_pval_SA_NPDs_med),
                     "NPD-umed" = fdr_block(df_pval_CT_NPDs_unmed, df_pval_SA_NPDs_unmed))
list_names = list("CNV" = colnames(df_pval_CT_CNVs), "NPD" = colnames(df_pval_CT_NPDs_all),
                  "NPD-med" = colnames(df_pval_CT_NPDs_med), "NPD-umed" = colnames(df_pval_CT_NPDs_unmed))

df_nsig_fdr = data.frame()
for (input_type in names(list_pval_adj)) {
  df_nsig_fdr = rbind(df_nsig_fdr,
                      data.frame(condition = list_names[[input_type]],
                                 type = input_type,
                                 n_regions = nrow(list_pval_adj[[input_type]]$CT),
                                 nsig_fdr_CT = colSums(as.matrix(list_pval_adj[[input_type]]$CT) < 0.05),
                                 nsig_fdr_SA = colSums(as.matrix(list_pval_adj[[input_type]]$SA) < 0.05),
                                 row.names = NULL))
}
print(df_nsig_fdr)


## Plot Supplementary Figure 10A-D
barplot_nsig = function(df, input_type) {
  df_subset = df[df$type == input_type, ]
  ggplot(df_subset, aes(x = condition)) +
    geom_bar(aes(y = nsig_fdr_CT), stat = "identity", fill = "salmon", color = "black") +
    geom_bar(aes(y = -nsig_fdr_SA), stat = "identity", fill = "skyblue", color = "black") +
    geom_text(aes(y = nsig_fdr_CT, label = nsig_fdr_CT), vjust = -1, size = 3) +
    geom_text(aes(y = -nsig_fdr_SA, label = nsig_fdr_SA), vjust = 1.5, size = 3) +
    scale_y_continuous(limits = c(-35, 35), breaks = seq(-30, 30, by = 10), labels = function(x) abs(x)) +
    labs(title = paste("FDR significant Counts -", input_type), x = NULL, y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

array_types = c("CNV", "NPD", "NPD-med", "NPD-umed")
array_panel_names = c("SupFig10A", "SupFig10B", "SupFig10C", "SupFig10D")
list_plots = lapply(array_types, function(t) barplot_nsig(df_nsig_fdr, t))

for (i in seq_along(array_types)) {
  save_panel(list_plots[[i]], array_panel_names[i], 6, 4.5)
}

p_supfig10 = plot_grid(plotlist = list_plots, nrow = 2, ncol = 2, labels = c("A.", "B.", "C.", "D."))
print(p_supfig10)

save_panel(p_supfig10, "SupFig10", 12, 9)


## Save source data
for (i in seq_along(array_types)) {
  save_source_data(df_nsig_fdr[df_nsig_fdr$type == array_types[i], ], array_panel_names[i])
}
