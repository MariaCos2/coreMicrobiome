##' @title func_rank_plot
##'
##' @description Plots TPM vs ranked otu and avg occ vs avg TPM.
##' @param functional_profile a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_tpm a value indicating whether the function is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param ...
##' @details Plots TPM vs ranked otu and avg occ vs avg TPM.
##' @return A scatter plot showing the relationship between TPM and ranked otu with color denoting and core or non-core functional otu and a line indicating the threshold and a scatter plot showing the relationship between average occurancy and average TPM with color denoting and core or non-core functional otu and a line indicating the threshold.
##' @examples
##'  func_rank_plot(functional_profile, sample, mini_tpm=10000, threshold=0.5)
##' @export

func_rank_plot <- function(functional_profile, sample, mini_tpm, threshold){

  func_core <- function_core(functional_profile, sample, mini_tpm, threshold)
  fun_tpm_sum <-  func_core[[1]]
  lastcall <- func_core[[3]]

  p1 <- ggplot() +
    geom_point(data = fun_tpm_sum, aes(x = dim, y = tpm_s, size = fun_rel, fill = fill), pch = 21, color = 'black', alpha = .8) +
    geom_vline(xintercept = lastcall, lty = 2, col = 'black', cex = 2) +
    labs(x = 'ranked OTUs', y = 'TPM') +
    scale_fill_manual(values = c("deeppink4", "darkseagreen3")) +
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
    ) + labs(size = "Avg TPM", fill = "Functional core")

  p2 <- ggplot() +
    geom_point(data = fun_tpm_sum, aes(x = fun_rel, y = fun_occ, size = fun_rel, fill = fill), pch = 21, color = 'black', alpha = .8) +
    labs(x = 'Avg TPM', y = 'Avg Occ') +
    scale_fill_manual(values = c("deeppink4", "darkseagreen3")) +
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
    ) + labs(size = "Avg TPM", fill = "Functional core")

  list(p1, p2)
}
