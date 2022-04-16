# library(propr)
# library(gridExtra)
# library(UpSetR)
# library(ggpubr)
# source("functions/utilizes.R")
# source("functions/networks.R")

##' @title networks_shared_ne_plot
##' 
##' @description Constructs the co-occurrence networks and calculates the common/shared nodes/edges among networks from different methods and visualization. 
##' @param otu a community count data matrix,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu, 
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param pre_threshold occurrence frequency,
##' @param fdr_threshold fdr threshold,
##' @param cor_threshold association threshold,
##' @param permutation permutation number of the association methods,
##' @param propr a logical value indicating whether to run propr,
##' @param sparcc a logical value indicating whether to run sparcc,
##' @param cclasso a logical value indicating whether to run cclasso,
##' @param propr_col the color of the properties of propr network,
##' @param sparcc_col the color of the properties of sparcc network,
##' @param cclasso_col the color of the properties of cclasso network,
##' @param ...
##' @details Constructs the co-occurrence networks and calculates the common/shared nodes/edges among networks from different methods and visualization. 
##' @return Upset plots showing the common/shared nodes/edges between networks.
##' @examples
##'  networks_shared_ne_plot(otu, taxa, sample, pre_threshold=0.5, fdr_threshold=0.1, cor_threshold=0.6, permutation=100, propr=TRUE, sparcc=TRUE, cclasso=TRUE, propr_col='darkcyan', sparcc_col='firebrick3', cclasso_col='goldenrod3')
##' @export

networks_shared_ne_plot <- function(otu, taxa, sample, pre_threshold, fdr_threshold, cor_threshold, permutation, propr = TRUE, sparcc = TRUE, cclasso = TRUE, propr_col, sparcc_col, cclasso_col){
  
  otu = read.csv(otu$datapath, sep = ",", header = T, row.names = 1)
  taxa = read.csv(taxa$datapath, sep = ",", header = T, row.names = 1)
  sample = read.csv(sample$datapath, sep = ",", header = T, row.names = 1)
  
  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa))
  sample = sample_data(sample)
  mydata = phyloseq(OTU, TAX, sample)
  
  mydata_filter = filter_prevalence(mydata, pre_threshold)
  otu <- as.data.frame(otu_table(mydata_filter))
  
  if (propr) {
    s1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)
    propr_nodes <- unique(union(s1[[1]]$from, s1[[1]]$to))
    propr_edge <- paste(s1[[1]]$from, s1[[1]]$to, sep = '-')
  }
  
  if (sparcc) {
    s2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)
    sparcc_nodes <- unique(union(s2[[1]]$from, s2[[1]]$to))
    sparcc_edge <- paste(s2[[1]]$from, s2[[1]]$to, sep = '-')
  }
  
  if (cclasso) {
    s3 <- get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, cclasso_col)
    cclasso_nodes <- unique(union(s3[[1]]$from, s3[[1]]$to))
    cclasso_edge <- paste(s3[[1]]$from, s3[[1]]$to, sep = '-')
  }
  
  if (propr && sparcc && cclasso){
    nodes_input <- c(
      propr = length(propr_nodes),
      sparcc = length(sparcc_nodes),
      cclasso = length(cclasso_nodes),
      "propr&sparcc&cclasso" = length(Reduce(intersect, list(propr_nodes, sparcc_nodes, cclasso_nodes))),
      "propr&sparcc" = length(intersect(propr_nodes, sparcc_nodes)),
      "propr&cclasso" = length(intersect(propr_nodes, cclasso_nodes)),
      "sparcc&cclasso" = length(intersect(sparcc_nodes, cclasso_nodes)))
    edges_input <- c(
      propr = length(propr_edge),
      sparcc = length(sparcc_edge),
      cclasso = length(cclasso_edge),
      "propr&sparcc&cclasso" = length(Reduce(intersect, list(propr_edge, sparcc_edge, cclasso_edge))),
      "propr&sparcc" = length(intersect(propr_edge, sparcc_edge)),
      "propr&cclasso" = length(intersect(propr_edge, cclasso_edge)),
      "sparcc&cclasso" = length(intersect(sparcc_edge, cclasso_edge)))
    p1 <- upset(fromExpression(nodes_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, sparcc_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Node number",
                query.legend = "bottom",
    ) 
    p2 <- upset(fromExpression(edges_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, sparcc_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Edge number",
                query.legend = "bottom",
    ) 
  }else if(!propr && sparcc && cclasso){
    nodes_input <- c(
      sparcc = length(sparcc_nodes),
      cclasso = length(cclasso_nodes),
      "sparcc&cclasso" = length(intersect(sparcc_nodes, cclasso_nodes)))
    edges_input <- c(
      sparcc = length(sparcc_edge),
      cclasso = length(cclasso_edge),
      "sparcc&cclasso" = length(intersect(sparcc_edge, cclasso_edge)))
    p1 <- upset(fromExpression(nodes_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(sparcc_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Node number",
                query.legend = "bottom",
    ) 
    p2 <- upset(fromExpression(edges_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(sparcc_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Edge number",
                query.legend = "bottom",
    ) 
  }else if(propr && !sparcc && cclasso){
    nodes_input <- c(
      propr = length(propr_nodes),
      cclasso = length(cclasso_nodes),
      "propr&cclasso" = length(intersect(propr_nodes, cclasso_nodes)))
    edges_input <- c(
      propr = length(propr_edge),
      cclasso = length(cclasso_edge),
      "propr&cclasso" = length(intersect(propr_edge, cclasso_edge)))
    p1 <- upset(fromExpression(nodes_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Node number",
                query.legend = "bottom",
    ) 
    p2 <- upset(fromExpression(edges_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, cclasso_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Edge number",
                query.legend = "bottom",
    ) 
  }else if(propr && sparcc && !cclasso){
    nodes_input <- c(
      propr = length(propr_nodes),
      sparcc = length(sparcc_nodes),
      "propr&sparcc" = length(intersect(propr_nodes, sparcc_nodes)))
    edges_input <- c(
      propr = length(propr_edge),
      sparcc = length(sparcc_edge),
      "propr&sparcc" = length(intersect(propr_edge, sparcc_edge)))
    p1 <- upset(fromExpression(nodes_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, sparcc_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Node number",
                query.legend = "bottom",
    ) 
    p2 <- upset(fromExpression(edges_input),
                text.scale = 2, 
                point.size = 3, 
                line.size = 2,
                sets.bar.color = c(propr_col, sparcc_col),
                main.bar.color = "coral2",
                mainbar.y.label = "Shared nodes",
                sets.x.label = "Edge number",
                query.legend = "bottom",
    ) 
  }else if(propr && !sparcc && !cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else if(!propr && sparcc && !cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else if(!propr && !sparcc && cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else(stop("Must provide at least two compositionary correlation method!!!"))
  
  
  # s1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)
  # s2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)
  # s3 <- get_cclasso_statistics(otu, fdr_threshold=1, cor_threshold=0.1, permutation, cclasso_col)
  # 
  # propr_nodes <- unique(union(s1[[1]]$from, s1[[1]]$to))
  # sparcc_nodes <- unique(union(s2[[1]]$from, s2[[1]]$to))
  # cclasso_nodes <- unique(union(s3[[1]]$from, s3[[1]]$to))
  # 
  # nodes_input <- c(
  #   propr = length(propr_nodes),
  #   sparcc = length(sparcc_nodes),
  #   cclasso = length(cclasso_nodes),
  #   "propr&sparcc&cclasso" = length(Reduce(intersect, list(propr_nodes, sparcc_nodes, cclasso_nodes))),
  #   "propr&sparcc" = length(intersect(propr_nodes, sparcc_nodes)),
  #   "propr&cclasso" = length(intersect(propr_nodes, cclasso_nodes)),
  #   "sparcc&cclasso" = length(intersect(sparcc_nodes, cclasso_nodes)))
  
  # p1 <- upset(fromExpression(nodes_input),
  #             text.scale = 2, 
  #             point.size = 3, 
  #             line.size = 2,
  #             sets.bar.color = c(propr_col, sparcc_col, cclasso_col),
  #             main.bar.color = "coral2",
  #             mainbar.y.label = "Shared nodes",
  #             sets.x.label = "Node number",
  #             query.legend = "bottom",
  # ) 
  p1 <- cowplot::plot_grid(NULL, 
                           p1$Main_bar, p1$Sizes, p1$Matrix,
                           nrow = 2, align = 'hv', rel_heights = c(3,1), rel_widths = c(2,3))
  
  # propr_edge <- paste(s1[[1]]$from, s1[[1]]$to, sep = '-')
  # sparcc_edge <- paste(s2[[1]]$from, s2[[1]]$to, sep = '-')
  # cclasso_edge <- paste(s3[[1]]$from, s3[[1]]$to, sep = '-')
  # 
  # edges_input <- c(
  #   propr = length(propr_edge),
  #   sparcc = length(sparcc_edge),
  #   cclasso = length(cclasso_edge),
  #   "propr&sparcc&cclasso" = length(Reduce(intersect, list(propr_edge, sparcc_edge, cclasso_edge))),
  #   "propr&sparcc" = length(intersect(propr_edge, sparcc_edge)),
  #   "propr&cclasso" = length(intersect(propr_edge, cclasso_edge)),
  #   "sparcc&cclasso" = length(intersect(sparcc_edge, cclasso_edge)))
  
  # p2 <- upset(fromExpression(edges_input),
  #             text.scale = 2, 
  #             point.size = 3, 
  #             line.size = 2,
  #             sets.bar.color = c(propr_col, sparcc_col, cclasso_col),
  #             main.bar.color = "coral2",
  #             mainbar.y.label = "Shared nodes",
  #             sets.x.label = "Edge number",
  #             query.legend = "bottom",
  # ) 
  p2 <- cowplot::plot_grid(NULL, 
                           p2$Main_bar, p2$Sizes, p2$Matrix,
                           nrow = 2, align = 'hv', rel_heights = c(3,1), rel_widths = c(2,3))
  
  grid.arrange(p1, p2, nrow = 2)
}