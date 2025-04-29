


#-------------------------------------------------------------------------------
# Utility functions for analysis and plotting results
#------------------------------------------------------------------------------- 
# SSDN
# Author: Kuldeep Kumar
# 

#-------------------------------------------------------------------------------
# Working directories
#------------------------------------------------------------------------------- 

## Set project Dir
project_dir = "temp_project"

# #----- directory structrue ----
code_dir <- paste0(project_dir,"/code")
data_dir <- paste0(project_dir,"/data")

plots_dir <- paste0(project_dir,"/plots_fig")
dir.create(plots_dir, showWarnings = FALSE)

# #-----  Set Working directory (Code directory)
in_wd <- paste0(project_dir,"/code")  
setwd(in_wd)


corr.type='pearson'
padj_method = 'fdr'

## Abagen gradient: But rotated to make it sensorimotor to association
input_gradient_profile = array_cortical_gradient

#-------------------------------------------------------------------------------
# PACKAGES
#-------------------------------------------------------------------------------
library(nlme)
library(xtable)
library(dplyr)
library(ggpubr)
library(RNOmni)
library(lme4)
library(ggplot2)
library("readxl")
library(ggseg)
library(RColorBrewer)
library(pals)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(corrplot)
library(matrixStats)
library(ggprism)
library(ggrepel)
library(ggExtra)



ggseg_desikan_rois = names(dk$palette)[c(1:3,5:35)]
ggseg_desikan_rois = gsub(x=ggseg_desikan_rois,pattern=" ",replacement = "")


#------------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#---- 0. Base function to adjust for p-value in matrix format-------
fPval_adj_in_mat <- function(in_pval_mat,padj_method='fdr'){
  # REF: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
  
  # p.adjust.methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
  in_pval_mat_FDR <- matrix(p.adjust(as.vector(as.matrix(in_pval_mat)), method=padj_method),ncol=ncol(in_pval_mat))
  
  return(in_pval_mat_FDR)
}




#---------------------------------------------------
## 1. load Null permutations
load(paste0("df_perm_ids_DesikanLH34_nIterNull_10000.RData"))   # LH: 34 ROIs


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
  
  # check if p-value is 0
  if(pval_avg == 0){
    pval_avg = 1/nperm
  }
  
  return(pval_avg)
  
}



### Function to compute cor
fCompute_cor_spin_p_mat_and_array = function(df_map_mat,array_2){
  
  corr.type='pearson'
  
  array_cor_vals = c()
  array_pspin_vals = c()
  for( loop_a in c(1:ncol(df_map_mat))){
    array_cor_vals = c(array_cor_vals,cor(df_map_mat[c(1:34),loop_a],array_2))
    array_pspin_vals = c(array_pspin_vals,perm.sphere.p(df_map_mat[c(1:34),loop_a],array_2,perm.id = df_perm_ids_LH,corr.type=corr.type))
  }
  
  df_cor_pspin = data.frame( map_names = colnames(df_map_mat),
                             cor = array_cor_vals,
                             pspin = array_pspin_vals)
  
  return(df_cor_pspin)
}



### Function to compute cor
fCompute_cor_spin_p_mat_and_mat2 = function(df_map_mat,df_map_mat2){
  
  corr.type='pearson'
  
  array_map_name1 = c()
  array_map_name2 = c()
  array_cor_vals = c()
  array_pspin_vals = c()
  for( loop_b in c(1:ncol(df_map_mat2))){
    array_2 = df_map_mat2[,loop_b]
    for( loop_a in c(1:ncol(df_map_mat))){
      array_map_name2 = c(array_map_name2,colnames(df_map_mat2)[loop_b])
      array_map_name1 = c(array_map_name1,colnames(df_map_mat)[loop_a])
      array_cor_vals = c(array_cor_vals,cor(df_map_mat[c(1:34),loop_a],array_2))
      array_pspin_vals = c(array_pspin_vals,perm.sphere.p(df_map_mat[c(1:34),loop_a],array_2,perm.id = df_perm_ids_LH,corr.type=corr.type))
    }
  }
  
  df_cor_pspin_mat1_v_mat2 = data.frame( map_name1 = array_map_name1,
                             map_name2 = array_map_name2,
                             cor = array_cor_vals,
                             pspin = array_pspin_vals)
  
  return(df_cor_pspin_mat1_v_mat2)
}




#---------------------------------------------------
## function to get mean-abs; var; %sig and PC1

fCor_with_gradient_in_df_es_pval = function(df_beta_cross_dis,df_pval_cross_dis,input_gradient_profile){
# 1. PCA
res.pca <- prcomp(df_beta_cross_dis, scale = TRUE,retx = TRUE)
temp_row <- get_pca_ind(res.pca)
row_loadings = temp_row$coord
temp_PC1 = row_loadings[,1]


#-------------- 2. Stats
temp_roi_names = rownames(df_pval_cross_dis)

## Transpose to make them 34 columns
df_beta_cross_dis = t(df_beta_cross_dis)
df_pval_cross_dis = t(df_pval_cross_dis)

p_val_threshold <- 0.05

df_summnary_nsig <- data.frame(
  roi = ggseg_desikan_rois,
  abs_meanEs = apply(df_beta_cross_dis, 2, function(x) mean(abs(x), na.rm = TRUE)),
  var = apply(df_beta_cross_dis, 2, function(x) sd(x, na.rm = TRUE)),
  nsig = apply(df_pval_cross_dis, 2, function(x) sum(x < p_val_threshold, na.rm = TRUE))
)

# Calculate percentages outside the loop for efficiency
df_summnary_nsig[,"nsig_percent"] <- 100 * df_summnary_nsig[,"nsig"] / nrow(df_pval_cross_dis)
if( cor(temp_PC1,df_summnary_nsig[,"var"]) < 0 ){
  df_summnary_nsig[,"PC1"] = -1*temp_PC1   # align the variance and PC1 direction
} else {
  df_summnary_nsig[,"PC1"] = temp_PC1
}
#df_summnary_nsig[,"PC1"] = temp_PC1
#df_summnary_nsig[,"gradient"] = df_neuromaps_abagen
temp_cor_w_gradient = cor(df_summnary_nsig[,c("abs_meanEs","var","nsig_percent","PC1")],input_gradient_profile)

return(temp_cor_w_gradient)

}



#---------------------------------------------------
## Corr + p-spin function to get mean-abs; var; %sig and PC1

fCor_pspin_with_gradient_in_df_es_pval = function(df_beta_cross_dis,df_pval_cross_dis,input_gradient_profile){
  # 1. PCA
  res.pca <- prcomp(df_beta_cross_dis, scale = TRUE,retx = TRUE)
  temp_row <- get_pca_ind(res.pca)
  row_loadings = temp_row$coord
  temp_PC1 = row_loadings[,1]
  
  
  #-------------- 2. Stats
  temp_roi_names = rownames(df_pval_cross_dis)
  
  ## Transpose to make them 34 columns
  df_beta_cross_dis = t(df_beta_cross_dis)
  df_pval_cross_dis = t(df_pval_cross_dis)
  
  p_val_threshold <- 0.05
  
  df_summnary_nsig <- data.frame(
    roi = ggseg_desikan_rois,
    abs_meanEs = apply(df_beta_cross_dis, 2, function(x) mean(abs(x), na.rm = TRUE)),
    var = apply(df_beta_cross_dis, 2, function(x) sd(x, na.rm = TRUE)),
    nsig = apply(df_pval_cross_dis, 2, function(x) sum(x < p_val_threshold, na.rm = TRUE))
  )
  
  # Calculate percentages outside the loop for efficiency
  df_summnary_nsig[,"nsig_percent"] <- 100 * df_summnary_nsig[,"nsig"] / nrow(df_pval_cross_dis)
  if( cor(temp_PC1,df_summnary_nsig[,"var"]) < 0 ){
    df_summnary_nsig[,"PC1"] = -1*temp_PC1   # align the variance and PC1 direction
  } else {
    df_summnary_nsig[,"PC1"] = temp_PC1
  }
  #df_summnary_nsig[,"PC1"] = temp_PC1
  #df_summnary_nsig[,"gradient"] = df_neuromaps_abagen
  #temp_cor_w_gradient = cor(df_summnary_nsig[,c("abs_meanEs","var","nsig_percent","PC1")],input_gradient_profile)
  
  ## cor and pspin
  df_map_mat = df_summnary_nsig[,c("abs_meanEs","var","nsig_percent","PC1")]
  array_2 = input_gradient_profile
  df_cor_pspin_w_gradient = fCompute_cor_spin_p_mat_and_array(df_map_mat,array_2)
  
  return(list(df_cor_pspin_w_gradient,df_summnary_nsig))
  
}


#------------------------------------------------------------------------------
# 2. Bootstrap function

fBootstrap_in_df_es_pval = function(df_es,df_pval,n_boot){
  
# Number of bootstrap iterations
#n_boot <- 1000

# Store results of each bootstrap iteration
results_cor <- matrix(NA, nrow = n_boot, ncol = 4)  # 4 columns for the output of fCor_with_gradient_in_df_es_pval
colnames(results_cor) <- c("abs_meanEs","var","nsig_percent","PC1")

### 1. CNVs
# Bootstrap loop
for (i in 1:n_boot) {
  # 1. Randomly select columns (with replacement)
  sampled_cols <- sample(1:ncol(df_es), size = ncol(df_es), replace = TRUE) # Sample 20 columns (adjust as needed)
  
  # 2. Extract the sampled columns from both dataframes
  es_subset <- df_es[, sampled_cols]
  pval_subset <- df_pval[, sampled_cols]
  
  # 3. Apply the analysis function to the subsets.
  analysis_output <- fCor_with_gradient_in_df_es_pval(es_subset, pval_subset, input_gradient_profile)
  results_cor[i, ] <- analysis_output
}

  return(results_cor)
}

#------------------------------------------------------------------------------
# 3. Leave 1 out

leave_one_out_analysis_in_df_es_pval = function(df_es,df_pval){
  
  n_cols <- ncol(df_es)  # Number of columns (ROIs)
  results_loo <- matrix(NA, nrow = n_cols, ncol = 4) # Store results for each LOO iteration
  colnames(results_loo) <- c("abs_meanEs","var","nsig_percent","PC1") # Name the columns
  
  for (i in 1:n_cols) {
    # 1. Leave one column (ROI) out
    loo_cols <- setdiff(1:n_cols, i)  # Indices of columns to *keep*
    
    # 2. Extract the subsets (leaving one out)
    es_subset <- df_es[, loo_cols]
    pval_subset <- df_pval[, loo_cols]
    
    # 3. Apply the analysis function to the subsets.
    analysis_output <- fCor_with_gradient_in_df_es_pval(es_subset, pval_subset, input_gradient_profile)
    results_loo[i, ] <- analysis_output
  }
  
  return(results_loo)
}


#------------------------------------------------------------------------------

#---- 3Function to make single Brain maps Variance or meanES
fBrainPlot_ggseg_Cortex_LH_RH_pos_MeanES_values <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  
  p_stacked = ggplot(data=someData) +
    geom_brain(
      atlas = dk, 
      #show.legend = FALSE,
      #hemi = 'left',
      position = position_brain(side ~ hemi ),
      #position = position_brain(hemi ~ side),
      aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
      colour = "grey",
      size = 0.3
    ) +
    theme_brain(text.family = "sans",text.size = 14)+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d'))(n = 90)),1)   # Use your custom colour palette
  
  
  #print(p_stacked)
  
  ## Add legend on bottom
  p_stacked_legend_bottom = p_stacked +  theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_bottom)
  
  # legend on right
  p_stacked_legend_right = p_stacked +  theme(legend.position = "right") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  #barwidth = 0.8,  # Extend bar length
                                  #barheight =5,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_right)
  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
  
  return(list_plot_and_legend)
}



#------------------------------------------------------------------------------

#---------------------

#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LH_RH_pos_SigMap_values <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  
  p_stacked = ggplot(data=someData) +
    geom_brain(
      atlas = dk, 
      #show.legend = FALSE,
      #hemi = 'left',
      position = position_brain(side ~ hemi ),
      #position = position_brain(hemi ~ side),
      aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
      colour = "grey",
      size = 0.3
    ) +
    theme_brain(text.family = "sans",text.size = 14)+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#fee5d9','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#99000d'))(n = 140)),1)   # Use your custom colour palette
  
  
  #print(p_stacked)
  
  ## Add legend on bottom
  p_stacked_legend_bottom = p_stacked +  theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_bottom)
  
  # legend on right
  p_stacked_legend_right = p_stacked +  theme(legend.position = "right") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  #barwidth = 0.8,  # Extend bar length
                                  #barheight =5,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_right)
  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
  
  return(list_plot_and_legend)
}


#------------------------------------------------------------------------------

#---- 3Function to make single Brain maps Variance or meanES
fBrainPlot_ggseg_Cortex_LH_RH_pos_VarMap_values <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  
  p_stacked = ggplot(data=someData) +
    geom_brain(
      atlas = dk, 
      #show.legend = FALSE,
      #hemi = 'left',
      position = position_brain(side ~ hemi ),
      #position = position_brain(hemi ~ side),
      aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
      colour = "grey",
      size = 0.3
    ) +
    theme_brain(text.family = "sans",text.size = 14)+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b'))(n = 90)),1)   # Use your custom colour palette
  
  
  
  #print(p_stacked)
  
  ## Add legend on bottom
  p_stacked_legend_bottom = p_stacked +  theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_bottom)
  
  # legend on right
  p_stacked_legend_right = p_stacked +  theme(legend.position = "right") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  #barwidth = 0.8,  # Extend bar length
                                  #barheight =5,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_right)
  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
  
  return(list_plot_and_legend)
}

#------------------------------------------------------------------------------

#---------------------

#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LH_RH_input_PC <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  p_left <-  ggseg(someData,
                   atlas = in_atlas_name,
                   position="dispersed",
                   #hemisphere = in_hem,
                   mapping =aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
                   colour = "grey",
                   size = 0.3,
                   adapt_scales = TRUE) +
    theme_brain(text.family = "sans",text.size = 14) + # "monospace") +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#5e3c99','#b2abd2','#f7f7f7','#fdb863','#e66101'))(n = 100)),0.9) +   # Use your custom colour palette
    #scale_fill_viridis(limits = c(cbar_min,cbar_max),option="plasma") +
    #scale_fill_gradientn(limits = c(cbar_min,cbar_max),
    #                     colours = colorRampPalette(c("#018571","#80CDC1","#F5F5F5","#DFC27D","#A6611A"))(n = 100)) +  #brewer.pal(5,"BrBG")
    #colours = colorRampPalette(c("#2B83BA","#ABDDA4","#FFFFBF","#FDAE61","#D7191C"))(n = 100)) +   # brewer.pal(5, "Spectral")
    #scale_fill_gradientn(colours = colorRampPalette(c("#DF8F44FF","#79AF97FF","#6A6599FF"))(n = 3)) +  # Use your custom colour palette
    #colours = colorRampPalette(c("#3C5488FF","white","#DC0000FF"))(n = 300)) +   # Nature color pallete based codes for blue and red
    #labs(fill=in_legtitle) + 
    theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = in_legtitle,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  # guides(fill = guide_colourbar(title = in_legtitle,
  #                               title.position = "top", #"left", # Legend title below bar
  #                               barwidth = 12, #7,  # Extend bar length
  #                               barheight = 1, #0.7,
  #                               # title.hjust =1,
  #                               # title.vjust = 0.9)) +
  #                               title.hjust =0.5,
  #                               title.vjust = 0)) +
  #   theme(legend.margin=margin(t = 0, unit='cm'))
  
  
  p_full = p_left
  
  # extract legend
  p_legend <- p_left
  cowplot_legend <- cowplot::get_legend(p_legend)
  
  # Remove legend and any other labels + remove some space around top and bottom 
  p_left <- p_left + theme(legend.position = "none")
  p_single_plot <- p_left + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_single_plot,cowplot_legend,p_full)
  
  return(list_plot_and_legend)
}





#------------------------------------------------------------------------------


#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LH_RH_Neuromaps_Gradient <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  p_left <-  ggseg(someData,
                   atlas = in_atlas_name,
                   position="dispersed",
                   #hemisphere = in_hem,
                   mapping =aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
                   colour = "grey",
                   size = 0.3,
                   adapt_scales = TRUE) +
    theme_brain(text.family = "sans",text.size = 14) + # "monospace") +
    #scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
    #                     colours = alpha(colorRampPalette(c('#5e3c99','#b2abd2','#f7f7f7','#fdb863','#e66101'))(n = 100)),0.9) +   # Use your custom colour palette
    #scale_fill_viridis(limits = c(cbar_min,cbar_max),option="plasma") +
    scale_fill_viridis(limits = c(cbar_min,cbar_max),option="magma") +
    #scale_fill_gradientn(limits = c(cbar_min,cbar_max),
    #                     colours = colorRampPalette(c("#018571","#80CDC1","#F5F5F5","#DFC27D","#A6611A"))(n = 100)) +  #brewer.pal(5,"BrBG")
    #colours = colorRampPalette(c("#2B83BA","#ABDDA4","#FFFFBF","#FDAE61","#D7191C"))(n = 100)) +   # brewer.pal(5, "Spectral")
    #scale_fill_gradientn(colours = colorRampPalette(c("#DF8F44FF","#79AF97FF","#6A6599FF"))(n = 3)) +  # Use your custom colour palette
    #colours = colorRampPalette(c("#3C5488FF","white","#DC0000FF"))(n = 300)) +   # Nature color pallete based codes for blue and red
    #labs(fill=in_legtitle) + 
    theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = in_legtitle,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  # guides(fill = guide_colourbar(title = in_legtitle,
  #                               title.position = "top", #"left", # Legend title below bar
  #                               barwidth = 12, #7,  # Extend bar length
  #                               barheight = 1, #0.7,
  #                               # title.hjust =1,
  #                               # title.vjust = 0.9)) +
  #                               title.hjust =0.5,
  #                               title.vjust = 0)) +
  #   theme(legend.margin=margin(t = 0, unit='cm'))
  
  
  p_full = p_left
  
  # extract legend
  p_legend <- p_left
  cowplot_legend <- cowplot::get_legend(p_legend)
  
  # Remove legend and any other labels + remove some space around top and bottom 
  p_left <- p_left + theme(legend.position = "none")
  p_single_plot <- p_left + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_single_plot,p_full)
  
  return(list_plot_and_legend)
}


#------------------------------------------------------------------------------

#---------------------

#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LH_RH_pos_neg_values <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  
  p_stacked = ggplot(data=someData) +
    geom_brain(
      atlas = dk, 
      #show.legend = FALSE,
      #hemi = 'left',
      position = position_brain(side ~ hemi ),
      #position = position_brain(hemi ~ side),
      aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
      colour = "grey",
      size = 0.3
    ) +
    theme_brain(text.family = "sans",text.size = 14)+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#5e3c99','#b2abd2','#f7f7f7','#fdb863','#e66101'))(n = 100)),1) 
  
  #print(p_stacked)
  
  ## Add legend on bottom
  p_stacked_legend_bottom = p_stacked +  theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_bottom)
  
  # legend on right
  p_stacked_legend_right = p_stacked +  theme(legend.position = "right") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  #barwidth = 0.8,  # Extend bar length
                                  #barheight =5,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_right)
  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
  
  return(list_plot_and_legend)
}

#-----------------------------------------------------------------------------
#---------------------
# REF: https://cran.r-project.org/web/packages/ggseg/vignettes/ggseg.html

#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LH_RH_pos_neg_values_1x4 <- function(array_es,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = dk$palette
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df <- as.data.frame(cbind(array_es1))
  colnames(temp_df) <- c("ES")
  
  ## Create columns for color and thickness based on the pvalue
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  # temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey", "black"))    # significnace color
  # temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.1", temp_df_es_pval$pval_log10))
  # 
  # Create data tibble with "region" column
  someData <- tibble(
    region = c(ggseg_region_names,ggseg_region_names),
    hemi  = c(rep("left",length(ggseg_region_names)),rep("right",length(ggseg_region_names))),
    effect = temp_df[,"ES"] #,
    # pval = temp_df_es_pval[,"pval"],
    # sig_col = temp_df_es_pval[,"sig_col"],
    # sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  
  p_stacked = ggplot(data=someData) +
    geom_brain(
      atlas = dk, 
      #show.legend = FALSE,
      hemi = 'left',
      #position = position_brain(side ~ hemi ),
      #position = position_brain(hemi ~ side),
      aes(fill = effect,size=I(0.3)), #,colour = I(sig_col),size = I(sig_thick)),
      colour = "grey",
      size = 0.3
    ) +
    theme_brain(text.family = "sans",text.size = 14)+
    theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),breaks=c(cbar_min,cbar_max),
                         colours = alpha(colorRampPalette(c('#5e3c99','#b2abd2','#f7f7f7','#fdb863','#e66101'))(n = 100)),1) 
  
  #print(p_stacked)
  
  ## Add legend on bottom
  p_stacked_legend_bottom = p_stacked +  theme(legend.position = "bottom") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  barwidth = 7,  # Extend bar length
                                  barheight =0.7,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_bottom)
  
  # legend on right
  p_stacked_legend_right = p_stacked +  theme(legend.position = "right") +  
    guides(fill = guide_colourbar(title = NULL,
                                  title.position = "top", # Legend title below bar
                                  #barwidth = 0.8,  # Extend bar length
                                  #barheight =5,
                                  title.hjust =0.5,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  
  #print(p_stacked_legend_right)
  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
  
  return(list_plot_and_legend)
}

#-----------------------------------------------------------------------------


#---- 2. Function to make single Brain maps for Cohen's D using ggseg-------
fBrainPlot_ggseg_Cortex_LHonly_input_ES_Pval <- function(array_es,array_pval,in_legtitle,in_hem,cbar_min,cbar_max){
  
  in_atlas_name = "dk"
  
  # limit within the max values of cbar
  array_es[array_es < cbar_min] = cbar_min
  array_es[array_es > cbar_max] = cbar_max
  
  pallete_ggseg = brain_pal(in_atlas_name)
  ggseg_region_names <- names(pallete_ggseg)
  #Note: ggseg dk has >"corpus callosum"   at index 4 (so add NA accordingly)
  if( in_hem == "left"){
    array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)])
    array_pval1 = c(array_pval[c(1:3)],1,array_pval[c(4:34)])
  } else if( in_hem == "right"){
    array_es1 = c(array_es[c(35:37)],NA,array_es[c(38:68)])
    array_pval1 = c(array_pval[c(35:37)],1,array_pval[c(38:68)])
  } else if( in_hem == "left_right"){
    array_es1 = c(array_es[c(1:3)],NA,array_es[c(4:34)],array_es[c(35:37)],NA,array_es[c(38:68)])
    array_pval1 = c(array_pval[c(1:3)],1,array_pval[c(4:34)],array_pval[c(35:37)],1,array_pval[c(38:68)])
  }
  
  
  # p-value cut-off for ROI boundary color
  cutoff_Pvalue = 0.05
  
  
  # Max value limit: keep estimates within cbar limit (# limit estimate values to max and min range for comparison)
  
  # create df with estimate and p-value
  temp_df_es_pval <- as.data.frame(cbind(array_es1,array_pval1))
  colnames(temp_df_es_pval) <- c("ES","pval")
  
  ## Create columns for color and thickness based on the pvalue
  temp_df_es_pval$pval_log10 <- -0.2*log10(array_pval1) # multiply by 0.05 => 4 will be 0.2
  temp_df_es_pval[temp_df_es_pval$pval_log10 > 0.3,"pval_log10"] <- 0.3 #0.15   # Max line width: any log10 p-value larger than 0.15 should be set at 0.15 (else lines are thick)
  colnames(temp_df_es_pval) <- c("ES","pval","pval_log10")   # column names for df
  
  #temp_df_es_pval$sig_col <- as.factor(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "white", "black"))   # significance color
  temp_df_es_pval$sig_col <- ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "grey60", "black")    # significnace color
  temp_df_es_pval$sig_thick <- as.numeric(ifelse(temp_df_es_pval$pval>cutoff_Pvalue, "0.05", temp_df_es_pval$pval_log10))
  
  # Create data tibble with "region" column
  someData <- tibble(
    region = ggseg_region_names, 
    effect = temp_df_es_pval[,"ES"],
    pval = temp_df_es_pval[,"pval"],
    sig_col = temp_df_es_pval[,"sig_col"],
    sig_thick = temp_df_es_pval[,"sig_thick"]
  )
  
  p_left <-  ggseg(someData,
                   atlas = in_atlas_name,
                   position="dispersed",
                   hemisphere = in_hem,
                   mapping =aes(fill = effect,colour = I(sig_col),size = I(sig_thick)),
                   #colour = "grey",
                   #size = 0.3,
                   adapt_scales = TRUE) +
    theme_brain(text.family = "sans",text.size = 14) + # "monospace") +
    scale_fill_gradientn(limits = c(cbar_min,cbar_max),
                         colours = coolwarm(100))+
    #colours = colorRampPalette(c("#2B83BA","#ABDDA4","#FFFFBF","#FDAE61","#D7191C"))(n = 100)) +   # brewer.pal(5, "Spectral")
    #scale_fill_gradientn(colours = colorRampPalette(c("#DF8F44FF","#79AF97FF","#6A6599FF"))(n = 3)) +  # Use your custom colour palette
    #colours = colorRampPalette(c("#3C5488FF","white","#DC0000FF"))(n = 300)) +   # Nature color pallete based codes for blue and red
    #labs(fill=in_legtitle) + 
    #theme(legend.position = "right") +  
    theme(legend.position = "bottom") +  
    # guides(fill = guide_colourbar(title = in_legtitle,
    #                               title.position = "left", # Legend title below bar
    #                               barwidth = 7,  # Extend bar length
    #                               barheight =0.7,
    #                               title.hjust =1,
    #                               title.vjust = 0.9)) +
    # theme(legend.margin=margin(t = 0, unit='cm'))
    # # Horiz: Bottom legend with Leg title on Top
    guides(fill = guide_colourbar(title = in_legtitle,
                                  title.position = "top", #"left", # Legend title below bar
                                  barwidth = 12, #7,  # Extend bar length
                                  barheight = 1, #0.7,
                                  title.hjust =1,
                                  title.vjust = 0.9)) +
    theme(legend.margin=margin(t = 0, unit='cm'))
  # # Vertical: right legend with Leg title on Top
  # guides(fill = guide_colourbar(title = in_legtitle,
  #                               title.position = "top", #"left", # Legend title below bar
  #                               barwidth = 1, #12, #7,  # Extend bar length
  #                               barheight = 8, #1, #0.7,
  #                               # title.hjust =1,
  #                               # title.vjust = 0.9)) +
  #                               title.hjust =0,
  #                               title.vjust = 0)) +
  #   theme(legend.margin=margin(t = 0, unit='cm'))
  
  # extract legend
  p_legend <- p_left
  
  # Remove legend and any other labels + remove some space around top and bottom 
  p_left <- p_left + theme(legend.position = "none")
  p_single_plot <- p_left + theme(axis.title = element_blank(),axis.text.x = element_blank()) + theme(plot.margin=grid::unit(c(-20,0,-20,0), "mm"))  
  
  # return single plot and legend
  list_plot_and_legend <- list(p_single_plot,p_legend)
  
  return(list_plot_and_legend)
}



#---------------------- PCA ------
## NOTE: We align the PCs to Variance explained per ROI profiles
fCompute_PCA_indf = function(in_df,plot_suffix,plots_dir){
  
  # 1. PCA
  res.pca <- prcomp(in_df, scale = TRUE,retx = TRUE)
  
  # 2. Plot eigen values and variance explained
  res_eig = get_eig(res.pca)
  
  p1 = fviz_eig(res.pca, col.var="blue",addlabels = TRUE)
  plot_file_name = paste0(plots_dir,"/fviz_eig_varexp_",plot_suffix,".pdf")
  pdf(plot_file_name,width=6, height=4)
  print(p1)
  dev.off()
  
  # 3. Column loadings
  temp_var = get_pca_var(res.pca)

  #------------------ Variable loadings ---------------------------------
  
  data_var_loadings <- temp_var$coord[,c(1,2)]
  
  colnames(data_var_loadings) <- c("PC1","PC2")
  
  #-------------------------------------------
  # 4. Row loadings
  temp_row <- get_pca_ind(res.pca)
  
  # return loadings coord
  col_loadings = temp_var$coord
  row_loadings = temp_row$coord
  
  # # 6. make brain-maps for PC1, PC2, and PC3
  cbar_min = -3
  cbar_max = 3
  
  for( loop_pc in c(1:2)){
    
    in_legtitle = paste0("PC",loop_pc)
    array_es = scale(row_loadings[,loop_pc]) # scale PCs
    
    # limit within the max values of cbar
    array_es[array_es < cbar_min] = cbar_min
    array_es[array_es > cbar_max] = cbar_max
    
    
    in_width = 6
    in_height =4
    in_hem = "left"
    ##list_plot_and_legend <- list(p_stacked_legend_right,p_stacked_legend_bottom)
    list_plots = fBrainPlot_ggseg_Cortex_LH_RH_pos_neg_values(array_es,in_legtitle,in_hem,cbar_min,cbar_max)
    
    pdf_filename = paste0(plots_dir,"/brianmap_PC_",loop_pc,"_",in_hem,"_",plot_suffix,"_legendRight.pdf")
    pdf(pdf_filename,width = in_width,height=in_height)
    print(list_plots[[1]])
    dev.off()
    
    pdf_filename = paste0(plots_dir,"/brianmap_PC_",loop_pc,"_",in_hem,"_",plot_suffix,"_legendBottom.pdf")
    pdf(pdf_filename,width = in_width,height=in_height)
    print(list_plots[[2]])
    dev.off()
    
    
  }
  
  
  # 7. return the col and row loadings
  res_list = list(col_loadings,row_loadings)
  
  return(res_list)
}


#-----------------------------------------------------------------------------





