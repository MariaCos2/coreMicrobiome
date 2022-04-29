##' @title species_acc_plot
##'
##' @description Plots species accumulation curve with specaccum function in vegan package
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param ...
##' @details Plots species accumulation curve with specaccum function in vegan package.
##' @return A species accumulation curve with boxplot indicating the 95% CI.
##' @examples
##'  species_acc_plot(otu)
##' @export

species_acc_plot <- function(OTU_input){
  data <- t(OTU_input)
  data_specacc <- specaccum(data, method = "random", permutations = 100)
  plot(data_specacc, ci.type = "poly", col = "coral4", lwd = 5, ci.lty = 0,
       ci.col = "lightgrey", main = "Species accumulation curve", xlab = "Samples", ylab = "Species")
  boxplot(data_specacc, col = "coral2", add = TRUE, pch = "+")
}
