#-------------------------------------------------------------------------------
# Run all figure code: main figures (Rmd), supplementary figures (R), Supplementary Table 7, Source Data
#-------------------------------------------------------------------------------
# Usage (from the repository folder): Rscript run_all.R
# Author: Kuldeep Kumar

setwd("code")
pdf(NULL)
dir.create("../figures", showWarnings = FALSE)
dir.create("../source_data", showWarnings = FALSE)

## Main figures
for (fig_file in c("Fig2.Rmd", "Fig3.Rmd", "Fig4.Rmd", "Fig5.Rmd")) {
  rmarkdown::render(fig_file, output_options = list(html_preview = FALSE), envir = new.env())
}

## Supplementary figures
for (fig_file in sort(list.files(".", pattern = "^SupFig[0-9]{2}_.*\\.R$"))) {
  source(fig_file, local = new.env())
}

## Supplementary Table 7
source("SupTable07_enrichment_FDR.R", local = new.env())

## Source Data
source("make_source_data.R", local = new.env())
