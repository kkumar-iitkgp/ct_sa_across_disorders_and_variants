#-------------------------------------------------------------------------------
# Supplementary Table 7: Enrichment of NPD-associated SNPs within global cortical morphology GWAS signals
#-------------------------------------------------------------------------------
# Legend: Median percentile of the independent genome-wide significant SNPs of ADHD (Demontis et al.,
# 2023), BD (O'Connell et al., 2025), MDD (Adams et al., 2025) and SCZ (Bigdeli et al., 2026) in the
# mean CT and total SA GWAS, compared with 10,000 random SNP sets matched on minor allele frequency and
# with 10,000 random SNP sets (one-sided permutation p-values). FDR q across the four disorders for
# each metric. Sensitivity analyses: SCZ, PGC3 European (Trubetskoy et al., 2022); MDD, without UK
# Biobank (P<5e-8); FDR q with the sensitivity set in place of the primary set of the same disorder.
# Abbreviations: ADHD=attention deficit hyperactivity disorder; BD=bipolar disorder; CT=cortical
# thickness; FDR=false discovery rate; GWAS=genome-wide association study; MDD=major depressive
# disorder; NPD=neurodevelopmental and psychiatric disorders; SA=surface area; SCZ=schizophrenia;
# SNP=single nucleotide polymorphism.
#
# Author: Kuldeep Kumar

source("Rscript_utility_functions.R")

tables_dir = "../tables/"
dir.create(tables_dir, showWarnings = FALSE)

## Load data
df_enrichment = read.csv(paste0(data_dir, "data_SupFig07_enrichment_tests.csv"))

## Function: FDR across the four disorders within each metric
fFDR_within_metric = function(df) {
  for (test_metric in unique(df[, "Metric"])) {
    index_metric = which(df[, "Metric"] == test_metric)
    df[index_metric, "FDR_MAFdecile_null"] = p.adjust(df[index_metric, "P_MAFdecile_null"], method = "fdr")
    df[index_metric, "FDR_random_null"] = p.adjust(df[index_metric, "P_random_null"], method = "fdr")
  }
  return(df)
}

## Stats Supplementary Table 7: primary sets
df_primary = fFDR_within_metric(df_enrichment[which(df_enrichment[, "Role"] == "primary"), ])

## Stats Supplementary Table 7: sensitivity sets, in place of the primary set of the same disorder
df_sensitivity = data.frame()
for (i in which(df_enrichment[, "Role"] == "sensitivity")) {
  df_swap = rbind(df_enrichment[which(df_enrichment[, "Role"] == "primary" & df_enrichment[, "Disorder"] != df_enrichment[i, "Disorder"]), ],
                  df_enrichment[i, ])
  df_swap = fFDR_within_metric(df_swap)
  df_sensitivity = rbind(df_sensitivity, df_swap[which(df_swap[, "Role"] == "sensitivity" & df_swap[, "Metric"] == df_enrichment[i, "Metric"]), ])
}

df_SupTable07 = rbind(df_primary, df_sensitivity)
df_SupTable07 = df_SupTable07[, c("Disorder", "Role", "Lead_set", "Metric", "N_lead_table", "N_in_Grasby", "Median_percentile",
                                  "P_MAFdecile_null", "FDR_MAFdecile_null", "P_random_null", "FDR_random_null",
                                  "Null_MAFdecile_median", "N_null", "Seed")]
print(df_SupTable07[, c("Disorder", "Role", "Metric", "N_in_Grasby", "Median_percentile", "P_MAFdecile_null", "FDR_MAFdecile_null")])

write.csv(df_SupTable07, paste0(tables_dir, "SupTable07_enrichment.csv"), row.names = FALSE)
