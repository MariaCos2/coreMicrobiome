##' @title count_func_compare
##'
##' @description Plots arrow plots showing the procruste rotation comparing otu and functions and the overlap between otu and functions representations using spls.
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param functional_profile a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu,
##' @param rank a taxonomic rank (Phylum or Class or Order or Family or Genus or Species),
##' @param ...
##' @details Plots arrow plots showing the procruste rotation comparing otu and functions and the overlap between otu and functions representations using spls.
##' @return Arrow plot showing the procruste rotation comparing otu and functions and the overlap between otu and functions representations using spls.
##' @examples
##'  count_func_compare(otu, funct, taxa, rank='Class')
##' @export

count_func_compare <- function(otu, funct, taxa, rank){

  otu <- as.data.frame(otu)

  funct <- as.data.frame(funct)
  funct_n <- funct[, !grepl("fun", colnames(funct))]
  funct_n <- funct_n %>%
    group_by(genome) %>%
    summarise_each(list(sum))
  names <- funct_n$genome
  funct_n <- funct_n[, !grepl("genome", colnames(funct_n))]
  funct_n <- funct_n[, colnames(otu)]
  funct_n <- as.data.frame(funct_n)
  row.names(funct_n) <- names
  funct <- funct_n[row.names(otu), ]

  taxa <- as.data.frame(taxa)
  tax <- taxa[, grepl(rank, colnames(taxa))]

  # mantel test
  otu_dist <- vegdist(otu, method = 'bray')
  func_dist <- vegdist(scale(funct), "euclid")
  mantel.results <- mantel(func_dist, otu_dist)
  print('mantel test...')
  print(mantel.results)

  # procrustes test
  otu_dist <- vegdist(otu, method = 'bray')
  otu_pcoa <- cmdscale(otu_dist, k = 2, eig = TRUE)
  func_hel <- decostand(funct, method = 'standardize')
  func_pca <- rda(func_hel)
  procrustes.results <- procrustes(X = otu_pcoa, Y = func_pca, symmetric = TRUE)
  o_pro <- as.data.frame(procrustes.results$Yrot)
  o_pro$pair <- seq(1, dim(o_pro)[1])
  o_pro$type <- 'taxa'
  colnames(o_pro) <- c('D1', 'D2', 'pair', 'type')
  o_pro <- cbind(o_pro, tax)
  p_pro <- as.data.frame(procrustes.results$X)
  p_pro$pair <- seq(1, dim(p_pro)[1])
  p_pro$type <- 'function'
  colnames(p_pro) <- c('D1', 'D2', 'pair', 'type')
  p_pro <- cbind(p_pro, tax)
  n_pro <- rbind(o_pro, p_pro)

  p1 <- n_pro %>%
    ggplot(aes(D1, D2)) +
    geom_point(aes(color = type, shape = tax), size = 3) +
    geom_line(aes(group = pair), color = "grey"#,arrow = arrow(type = "closed",length=unit(0.075, "inches"))
    ) +
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
    ) + labs(title = 'Procruste rotation comparing otu and functions', x = 'Comp1', y = 'Comp2')

  # protest test
  protest.results <- protest(X = otu_pcoa, Y = func_pca, scores = "sites", permutations = 999)
  print('protest.results...')
  print(protest.results)

  # spls
  MyResult.spls <- mixOmics::spls(otu, funct)
  o_spls <- MyResult.spls$variates$X
  o_spls <- as.data.frame(cbind(o_spls, tax))
  o_spls$comp1 <- as.numeric(o_spls$comp1)
  o_spls$comp2 <- as.numeric(o_spls$comp2)
  o_spls$pair <- seq(1, dim(o_spls)[1])
  o_spls$type <- 'taxa'
  p_spls <- MyResult.spls$variates$Y
  p_spls <- as.data.frame(cbind(p_spls, tax))
  p_spls$comp1 <- as.numeric(p_spls$comp1)
  p_spls$comp2 <- as.numeric(p_spls$comp2)
  p_spls$pair <- seq(1, dim(p_spls)[1])
  p_spls$type <- 'function'
  n_spls <- rbind(o_spls, p_spls)

  p2 <- n_spls %>%
    ggplot(aes(comp1, comp2)) +
    geom_point(aes(color = type, shape = tax), size = 3) +
    geom_line(aes(group = pair), color = "grey"#,arrow = arrow(type = "closed",length=unit(0.075, "inches"))
    ) +
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
    ) + labs(title = 'Overlap between otu and functions representations', x = 'Comp1', y = 'Comp2')

  list(p1, p2)
}


