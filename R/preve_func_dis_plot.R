##' @title preve_func_dis_plot
##'
##' @description Plots average occupancy and average TPM between core and non-core functional otu.
##' @param functional_profile a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_tpm a value indicating whether the function is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param core_col the color of the core,
##' @param noncore_col the color of the non-core,
##' @param ...
##' @details Plots average occupancy and average TPM between core and non-core functional otu.
##' @return Two density plots showing the average occupancy and average TPM between core and non-core functional otu.
##' @examples
##'  preve_func_dis_plot(functional_profile, sample, mini_abun=1000, threshold=0.5, core_col = "deeppink4", noncore_col = "darkseagreen3")
##' @export

preve_func_dis_plot <- function(functional_profile, sample, mini_tpm, threshold, core_col, noncore_col){
  
  fun_tpm_sum <- function_core(functional_profile, sample, mini_tpm, threshold)[[1]]
  
  p1 <- ggdensity(fun_tpm_sum, x = "fun_occ",
            add = "mean", rug = TRUE,
            color = "fill", fill = "fill",
            palette = c(core_col, noncore_col)) +
    labs(title = 'Avg occupancy of otu') +
    xlab('avg occupancy') +
    ylab('') +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, color = "black", hjust = 0.5, vjust = 1, lineheight = 0.2),
          axis.title.x = element_text(size = 14, color = "black", hjust = 0.5),
          axis.title.y = element_text(size = 14,color = "black", hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 14),
          text = element_text(size =14),
          legend.position = "right",
          legend.title = element_text(colour = "black", size = 14),
          legend.text = element_text(colour ="black", size = 14),
          panel.background = element_rect(colour = "black", size = 1)
    )
  
  p2 <- ggdensity(fun_tpm_sum, x = "fun_rel",
            add = "mean", rug = TRUE,
            color = "fill", fill = "fill",
            palette = c(core_col, noncore_col)) +
    labs(title = 'Avg TPM of otu') +
    xlab('avg TPM') +
    ylab('') +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, color = "black", hjust = 0.5, vjust = 1, lineheight = 0.2),
          axis.title.x = element_text(size = 14, color = "black", hjust = 0.5),
          axis.title.y = element_text(size = 14,color = "black", hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 14),
          text = element_text(size =14),
          legend.position = "right",
          legend.title = element_text(colour = "black", size = 14),
          legend.text = element_text(colour ="black", size = 14),
          panel.background = element_rect(colour = "black", size = 1)
    )
  
  list(p1, p2)
}