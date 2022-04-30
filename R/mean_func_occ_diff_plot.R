# There must be Sample_ID and SampleType in map file.

##' @title mean_func_occ_diff_plot
##'
##' @description Plots a heatmap showing the mean occurrence frequency difference between core and non-core functional otu under different minimal tpm and thresholds.
##' @param functional_profile a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_tpm the minimum TPM of the loop,
##' @param max_tpm the maximum TPM of the loop,
##' @param abun_step the TPM step of the loop,
##' @param mini_thre the minimum threshold of the loop,
##' @param max_thre the maximum threshold of the loop,
##' @param thre_step the threshold step of the loop,
##' @param ...
##' @details Plots a heatmap showing the mean occurrence frequency difference between core and non-core functional otu under different minimal tpm and thresholds.
##' @return A heatmap showing the mean occurrence frequency difference between core and non-core functional otu under different minimal tpm and thresholds.
##' @examples
##'  mean_func_occ_diff_plot(functional_profile, sample, min_tpm = 0, max_tpm = 4000, tpm_step = 100, min_thre = 0, max_thre = 0.1, thre_step = 0.02)
##' @export

mean_func_occ_diff_plot <- function(functional_profile, sample, min_tpm, max_tpm, tpm_step, min_thre, max_thre, thre_step){ # , height, weight

  options(warn=-1)

  p <- dis_func_occ_loop(functional_profile, sample, min_tpm, max_tpm, tpm_step, min_thre, max_thre, thre_step)
  p
}
