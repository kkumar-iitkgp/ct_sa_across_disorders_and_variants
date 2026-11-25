#-------------------------------------------------------------------------------
# Supplementary Figure 18: Cortical alterations across NPD-associated common variants and polygenic risk scores
#-------------------------------------------------------------------------------
# Legend: A) Percentage of negative (-ve) and positive (+ve) regression estimates (from cortical
# structure GWAS for mean CT and total SA) for the risk alleles of BD and SCZ genome-wide associated
# SNPs, after allele alignment (SA: 234 BD and 259 SCZ SNPs; CT: 225 and 249, excluding SNPs with a zero
# estimate; two-sided binomial test against 50%). B-C) CT and SA effect sizes for BD and SCZ polygenic
# risk scores (PRS) for B) global, and C) regional measures. Notably, none of the PRS effect sizes were
# FDR significant (panels B and C). Effect sizes for PRS were derived using a linear regression in UK
# Biobank participants, adjusting for age, sex, site, and ancestry principal components. Analyses were
# restricted to 31,413 UK Biobank participants of European ancestry who did not carry any of the
# recurrent CNVs examined in this study. Regional SA effect sizes were adjusted for total SA.
# Abbreviations: BD=bipolar disorder; CT=cortical thickness; PRS=polygenic risk scores; SA: surface
# area; SCZ=schizophrenia.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")
library(patchwork)

col_values_sign = c("Positive (+)" = "#E8890C", "Negative (-)" = "#2E8B57")
names(col_values_metric) = c("CT", "SA")

## Load data
df_sign_split = read.csv(paste0(data_dir, "data_SupFig18A_sign_split.csv"))
df_PRS = read.csv(paste0(data_dir, "data_PRS_BD_SCZ_global_regional.csv"))

## Stats Supplementary Figure 18A: two-sided binomial test of the proportion of negative estimates against 50%
df_sign_split[, "Binomial_p_two_sided"] = sapply(seq_len(nrow(df_sign_split)), function(i)
  signif(binom.test(df_sign_split[i, "N_negative"], df_sign_split[i, "N"], p = 0.5)$p.value, 2))
print(df_sign_split)

## Stats Supplementary Figure 18B and 18C: FDR across all PRS estimates (2 PRS x 2 metrics x 35 measures)
df_PRS[, "pval_adj"] = p.adjust(df_PRS[, "PRS_lm_pval"], method = "fdr")
df_PRS[, "FDR_significant"] = df_PRS[, "pval_adj"] < 0.05
print(df_PRS[which(df_PRS[, "Desikan_roi"] == "global"), c("Disorder_PRS", "metric", "PRS_lm_est", "PRS_lm_pval", "pval_adj")])
print(paste0("Number of FDR significant PRS estimates: ", sum(df_PRS[, "FDR_significant"])))

## Plot Supplementary Figure 18A: sign split of allele-aligned cortical GWAS estimates
df_sign_split[, "panel"] = paste0(df_sign_split[, "Disorder"], " - ", df_sign_split[, "Metric"], "\nn = ", df_sign_split[, "N"], " SNPs")
df_sign_split[, "panel"] = factor(df_sign_split[, "panel"], levels = df_sign_split[order(df_sign_split[, "Disorder"], df_sign_split[, "Metric"]), "panel"])

in_df_plot = rbind(data.frame(df_sign_split, sign = "Negative (-)", pct = df_sign_split[, "Pct_negative"]),
                   data.frame(df_sign_split, sign = "Positive (+)", pct = 100 - df_sign_split[, "Pct_negative"]))
in_df_plot[, "sign"] = factor(in_df_plot[, "sign"], levels = c("Positive (+)", "Negative (-)"))

p_SupFig18A = ggplot(in_df_plot, aes(x = 1, y = pct, fill = sign)) +
  geom_col(width = 1, colour = "white", na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4.5, na.rm = TRUE) +
  coord_polar(theta = "y", start = 0, direction = 1) +
  facet_wrap(~ panel, ncol = 2, drop = FALSE) +
  scale_x_continuous(limits = c(0.5, 1.5)) +
  scale_fill_manual(values = col_values_sign, name = NULL) +
  ggtitle("A. % of -ve and +ve beta estimates for BD/SCZ associated SNPs") +
  theme_void(base_size = input_base_font_size + 1) +
  theme(legend.position = "top", strip.text = element_text(face = "bold", margin = margin(b = 4)),
        plot.title = element_text(size = input_base_font_size + 1))

## Plot Supplementary Figure 18B: PRS effect sizes on global metrics
df_PRS_global = df_PRS[which(df_PRS[, "Desikan_roi"] == "global"), ]

p_SupFig18B = ggplot(df_PRS_global, aes(x = metric, y = PRS_lm_est, fill = metric)) +
  theme_prism(base_size = input_base_font_size + 4, base_line_size = 0.75) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(~ Disorder_PRS) +
  scale_fill_manual(values = col_values_metric) +
  scale_y_continuous(limits = c(-0.02, 0.02)) +
  xlab(NULL) + ylab("Effect size") +
  ggtitle("B. PRS BD/SCZ effect sizes on global metrics") +
  theme(legend.position = "none", strip.text = element_text(face = "plain"),
        plot.title = element_text(size = input_base_font_size + 1, face = "plain", hjust = 0))

## Plot Supplementary Figure 18C: PRS effect sizes on regional metrics (left hemisphere)
df_PRS_regional = df_PRS[which(df_PRS[, "Desikan_roi"] != "global"), ]
df_PRS_regional[, "region"] = names(dk$palette)[c(1:3, 5:35)][match(df_PRS_regional[, "Desikan_roi"], rois_desikan)]
df_PRS_regional[, "sig_col"] = ifelse(df_PRS_regional[, "FDR_significant"], "black", "grey")

p_SupFig18C = ggplot(dplyr::group_by(df_PRS_regional, Disorder_PRS, metric)) +
  geom_brain(atlas = dk, hemi = "left",
             aes(fill = PRS_lm_est, colour = I(sig_col))) +
  facet_grid(Disorder_PRS ~ metric, switch = "y") +
  scale_fill_gradientn(limits = c(-0.02, 0.02), colours = pals::coolwarm(100), na.value = "grey50") +
  guides(fill = guide_colourbar(title = "Effect size", title.position = "top", title.hjust = 0.5, barwidth = 12, barheight = 1)) +
  ggtitle("C. PRS BD/SCZ effect sizes on regional metrics") +
  theme_void(base_size = input_base_font_size + 1) +
  theme(legend.position = "bottom", strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(size = input_base_font_size + 1))

save_panel(p_SupFig18A, "SupFig18A", 6, 6.5)
save_panel(p_SupFig18B, "SupFig18B", 6, 6.5)
save_panel(p_SupFig18C, "SupFig18C", 12, 5.5)

p_SupFig18 = (p_SupFig18A | p_SupFig18B) / p_SupFig18C + plot_layout(heights = c(1, 0.85))
save_panel(p_SupFig18, "SupFig18", 12, 12.5)

## Save source data
save_source_data(df_sign_split[, c("Disorder", "Metric", "N", "N_negative", "Pct_negative", "Binomial_p_two_sided")], "SupFig18A")
save_source_data(df_PRS_global[, c("Disorder_PRS", "metric", "PRS_lm_est", "PRS_lm_pval", "pval_adj", "FDR_significant")], "SupFig18B")
save_source_data(df_PRS_regional[, c("Desikan_roi", "Disorder_PRS", "metric", "PRS_lm_est", "PRS_lm_pval", "pval_adj", "FDR_significant")], "SupFig18C")
