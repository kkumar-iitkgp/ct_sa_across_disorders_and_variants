#-------------------------------------------------------------------------------
# Source Data: combine the per-panel CSV files into one XLSX file
#-------------------------------------------------------------------------------
# Author: Kuldeep Kumar

library(openxlsx)

source_data_dir = "../source_data/"

## Figure descriptions (Index sheet)
df_fig_description = data.frame(
  figure = c("Fig2", "Fig3", "Fig4", "Fig5",
             "SupFig02", "SupFig03", "SupFig04", "SupFig05", "SupFig06", "SupFig07", "SupFig08", "SupFig09",
             "SupFig10", "SupFig11", "SupFig12", "SupFig13", "SupFig14", "SupFig15", "SupFig16", "SupFig17"),
  description = c(
    "Global CT and SA effect sizes: CNVs vs NPDs, medicated vs unmedicated NPDs",
    "Regional effect sizes and latent dimension (PC1) of CNVs",
    "Regional effect sizes and latent dimension (PC1) of NPDs",
    "Correlation between the cortical gradient and latent dimensions and heritability maps",
    "Signed CT vs SA effect sizes",
    "Paired comparison of absolute CT and SA effect sizes",
    "Deletions, duplications and sex chromosome aneuploidies",
    "Psychiatric risk (iPSYCH2015 hazard ratio) vs absolute effect sizes of CNVs",
    "Genetic correlations between NPDs and global cortical metrics (Grasby et al. 2020)",
    "Enrichment of NPD-associated SNPs within global cortical GWAS signals",
    "NPDs vs first-degree relatives",
    "Adult, young and pediatric age sub-groups",
    "Number of FDR significant regions per cortical map",
    "Correlations between consensus maps",
    "Cortical gradient vs regional twin and SNP heritability",
    "Regional SA associations of NPD-associated SNPs",
    "Cortical gradient vs CT consensus maps",
    "Cortical gradient vs SA consensus maps",
    "Latent dimensions (PC1) vs other reference maps of cortical organization",
    "Signs of SNP effects on CT and SA, and PRS effect sizes"))

## List panel files in figure order
panel_files = list.files(source_data_dir, pattern = "^(Fig|SupFig)[0-9].*\\.csv$")
panel_names = sub("\\.csv$", "", panel_files)
fig_names = sub("^((Sup)?Fig[0-9]+).*$", "\\1", panel_names)
panel_order = order(match(fig_names, df_fig_description$figure), panel_names)
panel_names = panel_names[panel_order]
fig_names = fig_names[panel_order]

## Index sheet
df_index = data.frame(sheet = panel_names,
                      figure = fig_names,
                      panel = sub("^(Sup)?Fig[0-9]+_?", "", panel_names),
                      description = df_fig_description$description[match(fig_names, df_fig_description$figure)])
df_index$panel[df_index$panel == ""] = "all"
print(df_index)

## Write workbook
wb = createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial")
addWorksheet(wb, "Index")
writeData(wb, "Index", df_index)
for (i in seq_along(panel_names)) {
  addWorksheet(wb, panel_names[i])
  writeData(wb, panel_names[i], read.csv(paste0(source_data_dir, panel_names[i], ".csv"), check.names = FALSE))
}
saveWorkbook(wb, paste0(source_data_dir, "SourceData_CommMed_R1.xlsx"), overwrite = TRUE)
