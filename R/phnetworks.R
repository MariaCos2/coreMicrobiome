# library(propr)
# library(gridExtra)
# library(ggtree)
# source("functions/utilizes.R")
# source("functions/cores.R")
# source("functions/networks.R")

##' @title phnetworks
##' 
##' @description Constructs a plot combining the phylogenetic tree, the heatmaps showing average occurrence frequency and average relative abundance of each taxa, and the co-occurrence networks from different methods. 
##' @param otu a community count data matrix,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu, 
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param tree a phylogenetic tree. The tip labels must match the otu names in the otu,
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
##' @param nscore the node-level property (degree, betweenness and closeness),
##' @param hcol the color of heatmap,
##' @param ttip the node size of the tree,
##' @param tlab the font size of the tree tip label,
##' @param hlab the font size of the heatmap label,
##' @param offset the offset between each subfigure,
##' @param width the width of each heatmap cell,
##' @param ...
##' @details Constructs a plot combining the phylogenetic tree, the heatmaps showing average occurrence frequency and average relative abundance of each taxa, and the co-occurrence networks from different methods. 
##' @return Combined plot showing the phylogenetic tree, the heatmap and co-occurrence networks.
##' @examples
##'  phnetworks(otu, taxa, sample, tree, propr = TRUE, sparcc = TRUE, cclasso = TRUE, pre_threshold, fdr_threshold, cor_threshold, permutation, nscore = 'degree', hcol = 'red', propr_col = 'darkcyan', sparcc_col = 'firebrick3', cclasso_col = 'goldenrod3', ttip = 3, tlab = 2, hlab = 4, offset = 0.3, width = 1)
##' @export


phnetworks <- function(otu, taxa, sample, tree, propr = TRUE, sparcc = TRUE, cclasso = TRUE, pre_threshold, fdr_threshold, cor_threshold, permutation, nscore = 'degree', hcol = 'red', propr_col = 'darkcyan', sparcc_col = 'firebrick3', cclasso_col = 'goldenrod3', ttip = 3, tlab = 2, hlab = 4, offset = 0.3, width = 1){
  
  otu = read.csv(otu$datapath, sep = ",", header = T, row.names = 1)
  taxa = read.csv(taxa$datapath, sep = ",", header = T, row.names = 1)
  sample = read.csv(sample$datapath, sep = ",", header = T, row.names = 1)
  tree = read.tree(tree$datapath)
  
  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa))
  sample = sample_data(sample)
  mydata = phyloseq(OTU, TAX, sample, tree)
  
  mydata_filter = filter_prevalence(mydata, pre_threshold)
  
  otu <- as.data.frame(otu_table(mydata_filter))
  tree <- phy_tree(mydata_filter) 

  
  if (!propr && !sparcc && !cclasso){
    stop("Must provide at least one compositionary correlation method!!!")
  }

  if (propr) {
    propr.results = get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, propr_col)
    propr.n <- propr.results$n
    propr.g <- propr.results$g
    if (nscore == 'degree'){
      propr.scores <- data.frame(names = V(propr.g)$name,
                                 local_score = degree(propr.g))
    }else if (nscore == 'betweenness'){
      propr.scores <- data.frame(names = V(propr.g)$name,
                                 local_score = betweenness(propr.g, normalized = T))
    }else if (nscore == 'closeness'){
      propr.scores <- data.frame(names = V(propr.g)$name,
                                 local_score = closeness(propr.g))
    }else (stop("Must be one of the degree, betweenness and closeness!!!"))
  }

  if (sparcc) {
    sparcc.results = get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, sparcc_col)
    sparcc.n <- sparcc.results$n
    sparcc.g <- sparcc.results$g
    if (nscore == 'degree'){
      sparcc.scores <- data.frame(names = V(sparcc.g)$name,
                                  local_score = degree(sparcc.g))
    }else if (nscore == 'betweenness'){
      sparcc.scores <- data.frame(names = V(sparcc.g)$name,
                                  local_score = betweenness(sparcc.g, normalized = T))
    }else if (nscore == 'closeness'){
      sparcc.scores <- data.frame(names = V(sparcc.g)$name,
                                  local_score = closeness(sparcc.g))
    }else (stop("Must be one of the degree, betweenness and closeness!!!"))
  }

  if (cclasso) {
    cclasso.results <- get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, cclasso_col)
    cclasso.n <- cclasso.results$n
    cclasso.g <- cclasso.results$g
    if (nscore == 'degree'){
      cclasso.scores <- data.frame(names = V(cclasso.g)$name,
                                   local_score = degree(cclasso.g))
    }else if (nscore == 'betweenness'){
      cclasso.scores <- data.frame(names = V(cclasso.g)$name,
                                   local_score = betweenness(cclasso.g, normalized = T))
    }else if (nscore == 'closeness'){
      cclasso.scores <- data.frame(names = V(cclasso.g)$name,
                                   local_score = closeness(cclasso.g))
    }else (stop("Must be one of the degree, betweenness and closeness!!!"))
  }


  cor_m <- c()
  cor_d <- c()
  cor_s <- c()

  if (propr && sparcc && cclasso){
    cor_m <- list('propr', 'sparcc', 'cclasso')
    cor_d <- list(dat1 = propr.n, dat2 = sparcc.n, dat3 = cclasso.n)
    cor_s <- list(dat1 = propr.scores, dat2 = sparcc.scores, dat3 = cclasso.scores)
  }else if(!propr && sparcc && cclasso){
    cor_m <- list('sparcc', 'cclasso')
    cor_d <- list(dat1 = sparcc.n, dat2 = cclasso.n)
    cor_s <- list(dat1 = sparcc.scores, dat2 = cclasso.scores)
  }else if(propr && !sparcc && cclasso){
    cor_m <- list('propr', 'cclasso')
    cor_d <- list(dat1 = propr.n, dat2 = cclasso.n)
    cor_s <- list(dat1 = propr.scores, dat2 = cclasso.scores)
  }else if(propr && sparcc && !cclasso){
    cor_m <- list('propr', 'sparcc')
    cor_d <- list(dat1 = propr.n, dat2 = sparcc.n)
    cor_s <- list(dat1 = propr.scores, dat2 = sparcc.scores)
  }else if(propr && !sparcc && !cclasso){
    cor_m <- list('propr')
    cor_d <- list(dat1 = propr.n)
    cor_s <- list(dat1 = propr.scores)
  }else if(!propr && sparcc && !cclasso){
    cor_m <- list('sparcc')
    cor_d <- list(dat1 = sparcc.n)
    cor_s <- list(dat1 = sparcc.scores)
  }else if(!propr && !sparcc && cclasso){
    cor_m <- list('cclasso')
    cor_d <- list(dat1 = cclasso.n)
    cor_s <- list(dat1 = cclasso.scores)
  }else(stop("Must provide at least one compositionary correlation method!!!"))

  BC_ranked <- common_core(otu, sample)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  otu_ranked <- BC_ranked$otu_ranked
  annotation_row = data.frame(ID = BC_ranked_abun$otu, Type = BC_ranked_abun$fill)
  cols <- c('core' = 'deeppink4', 'non-core' = 'darkseagreen3')

  p0 <- ggtree(tree)
  p1 <- p0  %<+%  annotation_row + geom_tippoint(aes(color = Type), size = ttip) +
    scale_color_manual(values = cols) + geom_tiplab(size = tlab, align = TRUE) + coord_fixed(ratio = 0.05, expand = TRUE) + theme(legend.position="none") #
  
  otu_relabun_plot = get_avg_rel_abun(otu, sample, otu_ranked, BC_ranked_abun)
  otu_relaocc_plot = get_avg_occ(otu, sample, otu_ranked, BC_ranked_abun)
  heatmap_coordinates = get_heatmap_coordinates(p1, otu_relabun_plot, otu_relaocc_plot, offset = offset, width = width)
  heatmap_coordinates_abun = heatmap_coordinates[1:(dim(heatmap_coordinates)[1]/2), ]
  heatmap_coordinates_occ = heatmap_coordinates[(dim(heatmap_coordinates)[1]/2+1):(dim(heatmap_coordinates)[1]), ]
  heatmap_la_coordinates_abun = data.frame(
    x = unique(heatmap_coordinates_abun$x),
    y = rep(min(heatmap_coordinates_abun$y), length(unique(heatmap_coordinates_abun$x))),
    lab = unique(heatmap_coordinates_abun$variable))
  heatmap_la_coordinates_occ = data.frame(
    x = unique(heatmap_coordinates_occ$x),
    y = rep(min(heatmap_coordinates_occ$y), length(unique(heatmap_coordinates_occ$x))),
    lab = unique(heatmap_coordinates_occ$variable))

  y_title <- max(as.data.frame(ggplot_build(p1)$data[length(ggplot_build(p1)$data)])$y) + (max(as.data.frame(ggplot_build(p1)$data[length(ggplot_build(p1)$data)])$y) - min(as.data.frame(ggplot_build(p1)$data[length(ggplot_build(p1)$data)])$y))/length(as.data.frame(ggplot_build(p1)$data[length(ggplot_build(p1)$data)])$y)*1.8
  x_htitle1 <- min(heatmap_coordinates_abun$x) + (max(heatmap_coordinates_abun$x) - min(heatmap_coordinates_abun$x))/2
  x_htitle2 <- min(heatmap_coordinates_occ$x) + (max(heatmap_coordinates_occ$x) - min(heatmap_coordinates_occ$x))/2

  p2 <- p1 + geom_tile(data = heatmap_coordinates_abun, aes(x, y, fill = value), color = "black") +
    geom_tile(data = heatmap_coordinates_occ, aes(x, y, fill = value), color = "black") +
    labs(fill = "") +
    geom_text(data = heatmap_la_coordinates_abun, aes(x, y, label = lab), size = 3, inherit.aes = FALSE, angle = 270, nudge_y = -1,  hjust = 0) +
    geom_text(data = heatmap_la_coordinates_occ, aes(x, y, label = lab), size = 3, inherit.aes = FALSE, angle = 270, nudge_y = -1,  hjust = 0) +
    scale_fill_gradient(guide = guide_legend(), low = "white", high = hcol) +
    annotate("text", x = x_htitle1, y = y_title, label = "Avg rel abun", size = hlab) +
    annotate("text", x = x_htitle2, y = y_title, label = "Avg occ", size = hlab)

  x_r <- max(as.data.frame(ggplot_build(p2)$data[length(ggplot_build(p2)$data)-2])$x) + offset #0.05
  links <- get_networks_coodinates(p2, x_r, cor_d$dat1)
  f <- data.frame(row.names = as.data.frame(ggplot_build(p2)$data[4])$label, y = as.integer(as.data.frame(ggplot_build(p2)$data[4])$y))
  f$names <- row.names(f)
  points <- data.frame(y = c(1:length(as.data.frame(ggplot_build(p2)$data[4])$label)), x = rep(x_r, length(as.data.frame(ggplot_build(p2)$data[4])$label)))
  points <- points %>% left_join(f, by = 'y') %>% left_join(cor_s$dat1, by = 'names')
  p <- p2 + geom_point(data = points, aes(x = x, y = y, size = local_score)) + labs(size = "Node score") +
    geom_curve(data = subset(links, ymin > ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = -0.5, ncp = 100, lineend = 'butt', color = propr_col) +
    geom_curve(data = subset(links, ymin < ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = 0.5, ncp = 100, lineend = 'butt', color = propr_col) +
    scale_alpha('Correlation', n.breaks = 20) +
    guides(alpha = guide_legend(reverse = TRUE, override.aes = list(shape = 22))) + theme(legend.position="top") +
    xlim(0, x_r + max(links$yma_ymi)*0.05*0.5*0.5) +
    ylim(-2, y_title + 1) +
    annotate("text", x = x_r + max(links$yma_ymi)*0.05*0.5*0.5*0.5, y = y_title, label = cor_m[1]) 

  if (length(cor_m) >= 2){
    x_r <- x_r + max(links$yma_ymi)*0.05*0.5*0.5 + offset # 0.05
    links <- get_networks_coodinates(p, x_r, cor_d$dat2)
    f <- data.frame(row.names = as.data.frame(ggplot_build(p)$data[4])$label, y = as.integer(as.data.frame(ggplot_build(p)$data[4])$y))
    f$names <- row.names(f)
    points <- data.frame(y = c(1:length(as.data.frame(ggplot_build(p)$data[4])$label)), x = rep(x_r, length(as.data.frame(ggplot_build(p)$data[4])$label)))
    points <- points %>% left_join(f, by = 'y') %>% left_join(cor_s$dat2, by = 'names')
    p <- p + geom_point(data = points, aes(x = x, y = y, size = local_score)) + labs(size = "Node score") +
      geom_curve(data = subset(links, ymin > ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = -0.5, ncp = 100, lineend = 'butt', color = sparcc_col) +
      geom_curve(data = subset(links, ymin < ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = 0.5, ncp = 100, lineend = 'butt', color = sparcc_col) +
      scale_alpha('Correlation', n.breaks = 20) +
      guides(alpha = guide_legend(reverse = TRUE, override.aes = list(shape = 22))) + theme(legend.position="top") +
      xlim(0, x_r + max(links$yma_ymi)*0.05*0.5*0.5) +
      ylim(-2, y_title + 1) +
      annotate("text", x = x_r + max(links$yma_ymi)*0.05*0.5*0.5*0.5, y = y_title, label = cor_m[2])
  } 

  if (length(cor_m) == 3){
    x_r <- x_r + max(links$yma_ymi)*0.05*0.5*0.5 + offset # 0.05
    links <- get_networks_coodinates(p, x_r, cor_d$dat3)
    f <- data.frame(row.names = as.data.frame(ggplot_build(p)$data[4])$label, y = as.integer(as.data.frame(ggplot_build(p)$data[4])$y))
    f$names <- row.names(f)
    points <- data.frame(y = c(1:length(as.data.frame(ggplot_build(p)$data[4])$label)), x = rep(x_r, length(as.data.frame(ggplot_build(p)$data[4])$label)))
    points <- points %>% left_join(f, by = 'y') %>% left_join(cor_s$dat3, by = 'names')
    p <- p + geom_point(data = points, aes(x = x, y = y, size = local_score)) + labs(size = "Node score") +
      geom_curve(data = subset(links, ymin > ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = -0.5, ncp = 100, lineend = 'butt', color = cclasso_col) +
      geom_curve(data = subset(links, ymin < ymax),aes(x = xmin, y = ymin, xend = xmax, yend = ymax, alpha = cor), angle = 90, curvature = 0.5, ncp = 100, lineend = 'butt', color = cclasso_col) +
      scale_alpha('Correlation', n.breaks = 20) +
      guides(alpha = guide_legend(reverse = TRUE, override.aes = list(shape = 22))) + theme(legend.position="top") +
      xlim(0, x_r + max(links$yma_ymi)*0.05*0.5*0.5) +
      ylim(-2, y_title + 1) +
      annotate("text", x = x_r + max(links$yma_ymi)*0.05*0.5*0.5*0.5, y = y_title, label = cor_m[3])
  } 


  p
  
}
  
  
  
  
  