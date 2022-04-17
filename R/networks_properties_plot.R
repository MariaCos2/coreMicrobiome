##' @title networks_properties_plot
##'
##' @description Constructs the co-occurrence networks and calculates the node-level and network-level properties and performs statistical testing and visualization.
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
##' @details Constructs the co-occurrence networks and calculates the node-level and network-level properties and performs statistical testing and visualization.
##' @return Bar chart and violin plots showing the comparison results between networks.
##' @examples
##'  networks_properties_plot(otu, taxa, sample, pre_threshold=0.5, fdr_threshold=0.1, cor_threshold=0.6, permutation=100, propr=TRUE, sparcc=TRUE, cclasso=TRUE, propr_col='darkcyan', sparcc_col='firebrick3', cclasso_col='goldenrod3')
##' @export

networks_properties_plot <- function(otu, taxa, sample, pre_threshold, fdr_threshold, cor_threshold, permutation, propr, sparcc, cclasso, propr_col, sparcc_col, cclasso_col){

  options(warn=-1)

  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa))
  sample = sample_data(sample)
  mydata = phyloseq(OTU, TAX, sample)

  mydata_filter = filter_prevalence(mydata, pre_threshold)
  otu <- as.data.frame(otu_table(mydata_filter))

  if (propr) {
    s1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)
  }

  if (sparcc) {
    s2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)
  }

  if (cclasso) {
    s3 <- get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, cclasso_col)
  }

  if (propr && sparcc && cclasso){
    f.global <- data.frame(
      value = c(s1[[2]]$smallworld, s2[[2]]$smallworld, s3[[2]]$smallworld, s1[[2]]$density, s2[[2]]$density, s3[[2]]$density, s1[[2]]$mean_dis, s2[[2]]$mean_dis, s3[[2]]$mean_dis, s1[[2]]$trans, s2[[2]]$trans, s3[[2]]$trans, s1[[2]]$trans_local, s2[[2]]$trans_local, s3[[2]]$trans_local),
      index = c(rep('smallworld', 3), rep('density', 3), rep('mean_dis', 3), rep('trans', 3), rep('trans_local', 3)),
      group = c(rep(c('propr', 'sparcc', 'cclasso'), 5))
    )
    f.local <- rbind(s1[[2]]$nodes, s2[[2]]$nodes, s3[[2]]$nodes)
    my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  }else if(!propr && sparcc && cclasso){
    f.global <- data.frame(
      value = c(s2[[2]]$smallworld, s3[[2]]$smallworld, s2[[2]]$density, s3[[2]]$density, s2[[2]]$mean_dis, s3[[2]]$mean_dis, s2[[2]]$trans, s3[[2]]$trans, s2[[2]]$trans_local, s3[[2]]$trans_local),
      index = c(rep('smallworld', 2), rep('density', 2), rep('mean_dis', 2), rep('trans', 2), rep('trans_local', 2)),
      group = c(rep(c('sparcc', 'cclasso'), 5))
    )
    f.local <- rbind(s2[[2]]$nodes, s3[[2]]$nodes)
    my_comparisons <- list(c("sparcc", "cclasso") )
  }else if(propr && !sparcc && cclasso){
    f.global <- data.frame(
      value = c(s1[[2]]$smallworld, s3[[2]]$smallworld, s1[[2]]$density, s3[[2]]$density, s1[[2]]$mean_dis, s3[[2]]$mean_dis, s1[[2]]$trans, s3[[2]]$trans, s1[[2]]$trans_local, s3[[2]]$trans_local),
      index = c(rep('smallworld', 2), rep('density', 2), rep('mean_dis', 2), rep('trans', 2), rep('trans_local', 2)),
      group = c(rep(c('propr', 'cclasso'), 5))
    )
    f.local <- rbind(s1[[2]]$nodes, s3[[2]]$nodes)
    my_comparisons <- list(c("propr", "cclasso"))
  }else if(propr && sparcc && !cclasso){
    f.global <- data.frame(
      value = c(s1[[2]]$smallworld, s2[[2]]$smallworld, s1[[2]]$density, s2[[2]]$density, s1[[2]]$mean_dis, s2[[2]]$mean_dis, s1[[2]]$trans, s2[[2]]$trans, s1[[2]]$trans_local, s2[[2]]$trans_local),
      index = c(rep('smallworld', 2), rep('density', 2), rep('mean_dis', 2), rep('trans', 2), rep('trans_local', 2)),
      group = c(rep(c('propr', 'sparcc'), 5))
    )
    f.local <- rbind(s1[[2]]$nodes, s2[[2]]$nodes)
    my_comparisons <- list(c("propr", "sparcc"))
  }else if(propr && !sparcc && !cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else if(!propr && sparcc && !cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else if(!propr && !sparcc && cclasso){
    stop("Must provide at least two compositionary correlation method!!!")
  }else(stop("Must provide at least two compositionary correlation method!!!"))

  p_glo <- ggplot(f.global, aes(index, value, fill = group)) + geom_bar(stat="identity",  position=position_dodge()) +
    scale_fill_manual(values=c(propr_col, sparcc_col, cclasso_col)) +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  p_de <- ggviolin(f.local, x = "group", y = "degree", fill = "group",
                   palette = c(propr_col, sparcc_col, cclasso_col),
                   add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  p_be <- ggviolin(f.local, x = "group", y = "betweenness", fill = "group",
                   palette = c(propr_col, sparcc_col, cclasso_col),
                   add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  p_cl <- ggviolin(f.local, x = "group", y = "closeness", fill = "group",
                   palette = c(propr_col, sparcc_col, cclasso_col),
                   add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  p_ei <- ggviolin(f.local, x = "group", y = "eigenvector", fill = "group",
                   palette = c(propr_col, sparcc_col, cclasso_col),
                   add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  p_hub <- ggviolin(f.local, x = "group", y = "hub", fill = "group",
                    palette = c(propr_col, sparcc_col, cclasso_col),
                    add = "boxplot", add.params = list(fill = "white")) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
    theme_minimal() +
    labs(y = 'hub score') +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
          axis.text.x = element_text(colour = "black", face = "bold", size = 12),
          legend.text = element_text(size = 12, face ="bold", colour ="black"),
          legend.position = "right",
          axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
          legend.title = element_text(size = 14, colour = "black", face = "bold"),
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  grid.arrange(p_glo, p_de, p_be, p_cl, p_ei, p_hub, nrow = 6)
  # s1 <- get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)
  # s2 <- get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)
  # s3 <- get_cclasso_statistics(otu, fdr_threshold=1, cor_threshold=0.1, permutation, cclasso_col)

  # f.global <- data.frame(
  #   value = c(s1[[2]]$smallworld, s2[[2]]$smallworld, s3[[2]]$smallworld, s1[[2]]$density, s2[[2]]$density, s3[[2]]$density, s1[[2]]$mean_dis, s2[[2]]$mean_dis, s3[[2]]$mean_dis, s1[[2]]$trans, s2[[2]]$trans, s3[[2]]$trans, s1[[2]]$trans_local, s2[[2]]$trans_local, s3[[2]]$trans_local),
  #   index = c(rep('smallworld', 3), rep('density', 3), rep('mean_dis', 3), rep('trans', 3), rep('trans_local', 3)),
  #   group = c(rep(c('propr', 'sparcc', 'cclasso'), 5))
  # )
  # p_glo <- ggplot(f.global, aes(index, value, fill = group)) + geom_bar(stat="identity",  position=position_dodge()) +
  #   scale_fill_manual(values=c(propr_col, sparcc_col, cclasso_col)) +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))

  # f.local <- rbind(s1[[2]]$nodes, s2[[2]]$nodes, s3[[2]]$nodes)
  # my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  # p_de <- ggviolin(f.local, x = "group", y = "degree", fill = "group",
  #                  palette = c(propr_col, sparcc_col, cclasso_col),
  #                  add = "boxplot", add.params = list(fill = "white")) +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
  #
  # my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  # p_be <- ggviolin(f.local, x = "group", y = "betweenness", fill = "group",
  #                  palette = c(propr_col, sparcc_col, cclasso_col),
  #                  add = "boxplot", add.params = list(fill = "white")) +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
  #
  # my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  # p_cl <- ggviolin(f.local, x = "group", y = "closeness", fill = "group",
  #                  palette = c(propr_col, sparcc_col, cclasso_col),
  #                  add = "boxplot", add.params = list(fill = "white")) +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
  #
  # my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  # p_ei <- ggviolin(f.local, x = "group", y = "eigenvector", fill = "group",
  #                  palette = c(propr_col, sparcc_col, cclasso_col),
  #                  add = "boxplot", add.params = list(fill = "white")) +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
  #
  # my_comparisons <- list( c("propr", "sparcc"), c("propr", "cclasso"), c("sparcc", "cclasso") )
  # p_hub <- ggviolin(f.local, x = "group", y = "hub", fill = "group",
  #                   palette = c(propr_col, sparcc_col, cclasso_col),
  #                   add = "boxplot", add.params = list(fill = "white")) +
  #   stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #   theme_minimal() +
  #   theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  #         axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  #         legend.text = element_text(size = 12, face ="bold", colour ="black"),
  #         legend.position = "right",
  #         axis.title.y = element_text(face = "bold", size = 14, colour = "black"),
  #         axis.title.x = element_text(face = "bold", size = 0, colour = "black"),
  #         legend.title = element_text(size = 14, colour = "black", face = "bold"),
  #         panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2))
  #
  # options(repr.plot.width = 15, repr.plot.height = 20)
  # grid.arrange(p_glo, p_de, p_be, p_cl, p_ei, p_hub, nrow = 6)
}
