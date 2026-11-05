#-------------------------------------------------------------------------------
# Utility functions for global and regional analyses (Figures 2-5, Supplementary Figures 2-17)
#-------------------------------------------------------------------------------
# Author: Kuldeep Kumar

library(ggplot2)
library(ggprism)
library(ggrepel)

## Paths (relative to the code/ folder)
data_dir = "../data/"
fig_dir = "../figures/"
source_data_dir = "../source_data/"

## Plot settings
input_base_font_size = 12
col_CT = "#e41a1c"
col_SA = "#377eb8"
col_values_metric = c(col_CT, col_SA)
col_values_CNV_NPD = c("CNV" = "#E7298A", "NPD" = "#1B9E77", "NPD-med" = "#D95F02", "NPD-umed" = "#7570B3")

## Desikan 34 regions (left hemisphere, ggseg order)
library(ggseg)
rois_desikan = gsub(" ", "", names(dk$palette)[c(1:3, 5:35)])


# # A .------ Spin Perm test --------------------
# Modified from REF: # Frantisek Vasa, fv247@cam.ac.uk, June 2017 - July 2018
perm.sphere.p = function(x,y,perm.id,corr.type='pearson') {

  nroi = dim(perm.id)[1]  # number of regions
  nperm = dim(perm.id)[2] # number of permutations

  rho.emp = cor(x,y,method=corr.type)  # empirical correlation

  # permutation of measures
  x.perm = y.perm = array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] = x[perm.id[i,r]]
      y.perm[i,r] = y[perm.id[i,r]]
    }
  }

  # correlation to unpermuted measures
  rho.null.xy = rho.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    rho.null.xy[r] = cor(x.perm[,r],y,method=corr.type)
    rho.null.yx[r] = cor(y.perm[,r],x,method=corr.type)
  }

  # p-value definition depends on the sign of the empirical correlation
  if (rho.emp>0) {
    p.perm.xy = sum(rho.null.xy>rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx>rho.emp)/nperm
  } else {
    p.perm.xy = sum(rho.null.xy<rho.emp)/nperm
    p.perm.yx = sum(rho.null.yx<rho.emp)/nperm
  }

  # return average p-value
  pval_avg <- (p.perm.xy+p.perm.yx)/2

  # p-value of 0 set to 1/nperm
  if(pval_avg == 0){
    pval_avg = 1/nperm
  }

  return(pval_avg)

}

## Save a panel as PNG (300 dpi) and SVG
save_panel = function(p, file_stem, width, height) {
  ggsave(paste0(fig_dir, file_stem, ".png"), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(fig_dir, file_stem, ".svg"), p, width = width, height = height, bg = "white")
}

## Save source data for one panel
save_source_data = function(df, panel_name) {
  write.csv(df, paste0(source_data_dir, panel_name, ".csv"), row.names = FALSE, na = "NA")
}
