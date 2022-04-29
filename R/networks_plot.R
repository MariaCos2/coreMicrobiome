##' @title networks_plot
##'
##' @description Constructs and plots co-occurrence networks based on specific occurrence frequency (pre_threshold), permutation number (permutation), fdr threshold (fdr_threshold), association threshold (cor_threshold) with propr, sparcc and cclasso.
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param pre_threshold occurrence frequency,
##' @param fdr_threshold fdr threshold,
##' @param cor_threshold association threshold,
##' @param permutation permutation number of the association methods,
##' @param propr a logical value indicating whether to run propr,
##' @param sparcc a logical value indicating whether to run sparcc,
##' @param cclasso a logical value indicating whether to run cclasso,
##' @param spieceasi a logical value indicating whether to run spieceasi,
##' @param propr_col the color of propr network,
##' @param sparcc_col the color of sparcc network,
##' @param cclasso_col the color of cclasso network,
##' @param spieceasi_col the color of spieceasi network,
##' @param ...
##' @details Constructs and plots networks based on specific occurrence frequency (pre_threshold), permutation number (permutation), fdr threshold (fdr_threshold), association threshold (cor_threshold) with propr, sparcc and cclasso.
##' @return Co-occurrence network plots showing the relationship between otu.
##' @examples
##'  networks_plot(otu, taxa, sample, pre_threshold=0.5, fdr_threshold=0.1, cor_threshold=0.6, permutation=100, propr=TRUE, sparcc=TRUE, cclasso=TRUE, spieceasi=TRUE, propr_col='darkcyan', sparcc_col='firebrick3', cclasso_col='goldenrod3', spieceasi_col = 'dodgerblue3')
##' @export


networks_plot <- function(otu, taxa, sample, pre_threshold, fdr_threshold, cor_threshold, permutation, propr, sparcc, cclasso, spieceasi, propr_col, sparcc_col, cclasso_col, spieceasi_col){

  options(warn=-1)

  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa))
  sample = sample_data(sample)
  mydata = phyloseq(OTU, TAX, sample)

  mydata_filter = filter_prevalence(mydata, pre_threshold)
  otu <- as.data.frame(otu_table(mydata_filter))

  if (!propr && !sparcc && !cclasso && !spieceasi){
    stop("Must provide at least one compositionary correlation method!!!")
  }

  if (propr) {
    p1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)$p
  }

  if (sparcc) {
    p2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)$p
  }

  if (cclasso) {
    p3 <- get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, cclasso_col)$p
  }

  if (spieceasi) {
    p4 <- get_spieceasi_statistics(otu, cor_threshold, permutation, spieceasi_col)$p
  }

  if (propr && sparcc && cclasso && spieceasi){
    grid.arrange(p1, p2, p3, p4, nrow = 4)
  }else if(propr && sparcc && cclasso && !spieceasi){
    grid.arrange(p1, p2, p3, nrow = 3)
  }else if(!propr && sparcc && cclasso && spieceasi){
    grid.arrange(p2, p3, p4, nrow = 3)
  }else if(propr && !sparcc && cclasso && spieceasi){
    grid.arrange(p1, p3, p4, nrow = 3)
  }else if(propr && sparcc && !cclasso && spieceasi){
    grid.arrange(p1, p2, p4, nrow = 3)
  }else if(!propr && sparcc && cclasso && !spieceasi){
    grid.arrange(p2, p3, nrow = 2)
  }else if(propr && !sparcc && cclasso && !spieceasi){
    grid.arrange(p1, p3, nrow = 2)
  }else if(propr && sparcc && !cclasso && !spieceasi){
    grid.arrange(p1, p2, nrow = 2)
  }else if(propr && !sparcc && !cclasso && spieceasi){
    grid.arrange(p1, p4, nrow = 2)
  }else if(!propr && sparcc && !cclasso && spieceasi){
    grid.arrange(p2, p4, nrow = 2)
  }else if(!propr && !sparcc && cclasso && spieceasi){
    grid.arrange(p3, p4, nrow = 2)
  }else if(propr && !sparcc && !cclasso  && !spieceasi){
    p1
  }else if(!propr && sparcc && !cclasso  && !spieceasi){
    p2
  }else if(!propr && !sparcc && cclasso  && !spieceasi){
    p3
  }else if(!propr && !sparcc && !cclasso  && spieceasi){
    p4
  }else(stop("Must provide at least one compositionary correlation method!!!"))


  # p1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)$p
  # p2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)$p
  # p3 <- get_cclasso_statistics(otu, fdr_threshold=1, cor_threshold=0.1, permutation, cclasso_col)$p
  #
  # options(repr.plot.width = 15, repr.plot.height = 20)
  # grid.arrange(p1, p2, p3, nrow = 3)

}
