# Copy number variants reveal divergent genetic and diagnostic cortical signatures across psychiatric disorders

Code and summary-level data to reproduce the figures of:

Kumar K. et al. *Copy number variants reveal divergent genetic and diagnostic cortical signatures across psychiatric disorders.* Communications Medicine (under review).

## Figures
| Figure | Code |
|---|---|
| Fig. 2: Global cortical differences across genetic variants and psychiatric diagnoses | `code/Fig2.Rmd` ([output](code/Fig2.md)) |
| Fig. 3: Regional cortical alterations across rare genetic variants | `code/Fig3.Rmd` ([output](code/Fig3.md)) |
| Fig. 4: Regional cortical alterations across neurodevelopmental and psychiatric disorders | `code/Fig4.Rmd` ([output](code/Fig4.md)) |
| Fig. 5: Correlation between the cortical gradient and latent dimensions of regional cortical differences | `code/Fig5.Rmd` ([output](code/Fig5.md)) |
| Supplementary Figures 2 to 17 | `code/SupFig02_*.R` to `code/SupFig17_*.R` |
| Supplementary Table 7: enrichment of NPD-associated SNPs in cortical GWAS signals | `code/SupTable07_enrichment_FDR.R` |

Figure 1 and Supplementary Figure 1 are overview illustrations.

## Run
```
Rscript run_all.R
```
This knits the main figure notebooks (requires pandoc, included with RStudio), runs the supplementary figure scripts, saves each panel to `figures/` (PNG and SVG) and the plotted values to `source_data/` (one CSV per panel, combined in `SourceData_CommMed_R1.xlsx`).

Shared functions (spin permutation test, colours, plot settings) are in `code/Rscript_utility_functions.R`.

## Data
`data/` holds summary-level inputs only: effect sizes and p-values of CNVs and psychiatric disorders, regional maps, principal component loadings, percentile ranks of disorder-associated SNPs in the cortical GWAS and enrichment test results, counts of SNP effect signs, published heritability and genetic correlation estimates, PRS model estimates, reference maps of cortical organization and the precomputed spin permutations (Desikan atlas, left hemisphere, 10,000 rotations). No individual-level data or GWAS summary statistics are included.

Only derived, summary-level data are shared here. Individual-level data from UK Biobank and the clinical CNV cohorts, and GWAS summary statistics from the PGC and ENIGMA, must be obtained from the original sources and used following their data access guidelines.


## Data sources

### Cortical GWAS (mean CT, total SA and regional SA)
- Grasby, K. L. et al. The genetic architecture of the human cerebral cortex. *Science* 367, eaay6690 (2020).

### Disorder-associated SNPs (independent genome-wide significant SNPs)
- ADHD: Demontis, D. et al. Genome-wide analyses of ADHD identify 27 risk loci, refine the genetic architecture and implicate several cognitive domains. *Nat. Genet.* 55, 198–208 (2023).
- BD: O'Connell, K. S. et al. Genomics yields biological and phenotypic insights into bipolar disorder. *Nature* 639, 968–975 (2025).
- MDD: Adams, M. J. et al. Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies. *Cell* 188, 640–652 (2025). The largest MDD GWAS from the PGC (European-ancestry COJO SNPs).
- MDD (regional SNP maps): Meng, X. et al. Multi-ancestry genome-wide association study of major depression aids locus discovery, fine mapping, gene prioritization and causal inference. *Nat. Genet.* 56, 222–233 (2024).
- SCZ: Bigdeli, T. B. et al. Biological insights into schizophrenia from ancestrally diverse populations. *Nature* 651, 404–413 (2026). The largest schizophrenia meta-analysis, combining the Million Veteran Program, All of Us and PGC data (European-ancestry COJO SNPs).
- SCZ (PGC3 European): Trubetskoy, V. et al. Mapping genomic loci implicates genes and synaptic biology in schizophrenia. *Nature* 604, 502–508 (2022).

### Psychiatric disorders (ENIGMA)
- Thompson, P. M. et al. ENIGMA and global neuroscience: A decade of large-scale studies of the brain in health and disease across more than 40 countries. *Transl. Psychiatry* 10, 100 (2020).
- ADHD: Hoogman, M. et al. *Am. J. Psychiatry* 176, 531–542 (2019).
- ASD: van Rooij, D. et al. *Am. J. Psychiatry* 175, 359–369 (2018).
- ADHD, ASD and OCD: Boedhoe, P. S. W. et al. *Am. J. Psychiatry* 177, 834–843 (2020).
- BD: Hibar, D. P. et al. *Mol. Psychiatry* 23, 932–942 (2018).
- CD: Gao, Y. et al. *Lancet Psychiatry* 11, 620–632 (2024).
- CHR-PS: ENIGMA Clinical High Risk for Psychosis Working Group. *JAMA Psychiatry* 78, 753–766 (2021).
- MDD: Schmaal, L. et al. *Mol. Psychiatry* 22, 900–909 (2017).
- OCD: Boedhoe, P. S. W. et al. *Am. J. Psychiatry* 175, 453–462 (2018).
- SCZ: van Erp, T. G. M. et al. *Biol. Psychiatry* 84, 644–654 (2018).
- First-degree relatives: de Zwarte, S. M. C. et al. *Biol. Psychiatry* 86, 545–556 (2019).
- CNVs: Sønderby, I. E. et al. *Hum. Brain Mapp.* 43, 300–328 (2022).

### UK Biobank polygenic risk scores
- Standard PRS for BD and SCZ: Thompson, D. J. et al. A systematic evaluation of the performance and properties of the UK Biobank Polygenic Risk Score (PRS) Release. *PLoS One* 19, e0307270 (2024).
- Sudlow, C. et al. UK Biobank: an open access resource for identifying the causes of a wide range of complex diseases of middle and old age. *PLoS Med.* 12, e1001779 (2015).


## Requirements
R (4.2) with ggplot2, ggprism, ggrepel, ggseg, cowplot, patchwork, rmarkdown, knitr, svglite and openxlsx.

## Interactive results
An interactive website to explore the main results of the paper:
https://kkumar-iitkgp.github.io/cortical-signatures-CNVs-NPDs/

## Pre-print
Cortical differences across psychiatric disorders and associated common and rare genetic variants. doi: https://doi.org/10.1101/2025.04.16.25325971
