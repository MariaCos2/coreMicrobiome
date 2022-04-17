# library(dplyr)
# library(ggpubr)
# library(plotly)
# library(gridExtra)
# source("functions/cores.R")

# There must be Sample_ID and SampleType in map file.

##' @title mean_occ_diff_plot
##' 
##' @description Plots a heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun the minimum abundance of the loop,
##' @param max_abun the maximum abundance of the loop,
##' @param abun_step the abundance step of the loop,
##' @param mini_thre the minimum threshold of the loop,
##' @param max_thre the maximum threshold of the loop,
##' @param thre_step the threshold step of the loop,
##' @param height the height of the figure,
##' @param weight the width of the figure,
##' @param ...
##' @details Plots a heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
##' @return A heatmap showing the mean occurrence frequency difference between core and non-core otu under different minimal abundance and threshold.
##' @examples
##'  mean_occ_diff_plot(otu, sample, min_abun = 0, max_abun = 400, abun_step = 100, min_thre = 0, max_thre = 0.1, thre_step = 0.02)
##' @export

mean_occ_diff_plot <- function(otu, sample, min_abun, max_abun, abun_step, min_thre, max_thre, thre_step){ # , height, weight
  
  otu = read.csv(otu, sep = ",", header = T, row.names = 1)
  sample = read.csv(sample, sep = ",", header = T, row.names = 1)
  p <- dis_occ_loop(otu, sample, min_abun = 0, max_abun = 400, abun_step = 100, min_thre = 0, max_thre = 0.1, thre_step = 0.02)
  p
}

