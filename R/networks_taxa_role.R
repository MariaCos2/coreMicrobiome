##' @title networks_taxa_role
##'
##' @description Determine the role of the taxa in a network.
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param pre_threshold occurrence frequency,
##' @param fdr_threshold fdr threshold,
##' @param cor_threshold association threshold,
##' @param cluster_method module clustering methods, including cluster_walktrap, cluster_edge_betweenness, cluster_fast_greedy and cluster_spinglass,
##' @param permutation permutation number of the association methods,
##' @param propr a logical value indicating whether to run propr,
##' @param sparcc a logical value indicating whether to run sparcc,
##' @param cclasso a logical value indicating whether to run cclasso,
##' @param ...
##' @details Determine the role of the taxa in a network.
##' @return A dataframe containing the role of the taxa in a network.
##' @examples
##'  networks_taxa_role(otu, taxa, sample, pre_threshold=0.5, fdr_threshold=0.1, cor_threshold=0.6, permutation=100, cluster_method="cluster_walktrap", propr=TRUE, sparcc=TRUE, cclasso=TRUE)
##' @export

networks_taxa_role <- function(otu, taxa, sample, pre_threshold, fdr_threshold, cor_threshold, permutation, cluster_method, propr = TRUE, sparcc = TRUE, cclasso = TRUE){

  options(warn=-1)

  OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa))
  sample = sample_data(sample)
  mydata = phyloseq(OTU, TAX, sample)

  mydata_filter = filter_prevalence(mydata, pre_threshold)
  otu <- as.data.frame(otu_table(mydata_filter))

  if (propr && sparcc && cclasso){
    propr_gt = get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    propr_t = ZiPivalue(igraph = propr_gt, method = cluster_method)
    propr_t <- propr_t %>% dplyr::select(propr_roles = roles) %>% mutate(taxa = row.names(propr_t))
    sparcc_gt = get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    sparcc_t = ZiPivalue(igraph = sparcc_gt, method = cluster_method)
    sparcc_t <- sparcc_t %>% dplyr::select(sparcc_roles = roles) %>% mutate(taxa = row.names(sparcc_t))
    cclasso_gt = get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    cclasso_t = ZiPivalue(igraph = cclasso_gt, method = cluster_method)
    cclasso_t <- cclasso_t %>% dplyr::select(cclasso_roles = roles) %>% mutate(taxa = row.names(cclasso_t))
    df_list <- list(propr_t, sparcc_t, cclasso_t)
    ros <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  }else if (propr && sparcc && !cclasso){
    propr_gt = get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    propr_t = ZiPivalue(igraph = propr_gt, method = cluster_method)
    propr_t <- propr_t %>% dplyr::select(propr_roles = roles) %>% mutate(taxa = row.names(propr_t))
    sparcc_gt = get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    sparcc_t = ZiPivalue(igraph = sparcc_gt, method = cluster_method)
    sparcc_t <- sparcc_t %>% dplyr::select(sparcc_roles = roles) %>% mutate(taxa = row.names(sparcc_t))
    df_list <- list(propr_t, sparcc_t)
    ros <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  }else if (propr && !sparcc && cclasso){
    propr_gt = get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    propr_t = ZiPivalue(igraph = propr_gt, method = cluster_method)
    propr_t <- propr_t %>% dplyr::select(propr_roles = roles) %>% mutate(taxa = row.names(propr_t))
    cclasso_gt = get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    cclasso_t = ZiPivalue(igraph = cclasso_gt, method = cluster_method)
    cclasso_t <- cclasso_t %>% dplyr::select(cclasso_roles = roles) %>% mutate(taxa = row.names(cclasso_t))
    df_list <- list(propr_t, cclasso_t)
    ros <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  }else if(!propr && sparcc && cclasso){
    sparcc_gt = get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    sparcc_t = ZiPivalue(igraph = sparcc_gt, method = cluster_method)
    sparcc_t <- sparcc_t %>% dplyr::select(sparcc_roles = roles) %>% mutate(taxa = row.names(sparcc_t))
    cclasso_gt = get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    cclasso_t = ZiPivalue(igraph = cclasso_gt, method = cluster_method)
    cclasso_t <- cclasso_t %>% dplyr::select(cclasso_roles = roles) %>% mutate(taxa = row.names(cclasso_t))
    df_list <- list(sparcc_t, cclasso_t)
    ros <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  }else if(propr && !sparcc && !cclasso){
    propr_gt = get_propr_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    propr_t = ZiPivalue(igraph = propr_gt, method = cluster_method)
    ros <- propr_t %>% dplyr::select(propr_roles = roles) %>% mutate(taxa = row.names(propr_t))
  }else if (!propr && sparcc && !cclasso){
    sparcc_gt = get_sparcc_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    sparcc_t = ZiPivalue(igraph = sparcc_gt, method = cluster_method)
    ros <- sparcc_t %>% dplyr::select(sparcc_roles = roles) %>% mutate(taxa = row.names(sparcc_t))
  }else if (!propr && !sparcc && cclasso){
    cclasso_gt = get_cclasso_statistics(otu, fdr_threshold, cor_threshold, permutation, 'red')$g
    cclasso_t = ZiPivalue(igraph = cclasso_gt, method = cluster_method)
    ros <- cclasso_t %>% dplyr::select(cclasso_roles = roles) %>% mutate(taxa = row.names(cclasso_t))
  }else(stop("Must provide at least one compositionary correlation method!!!"))

  print('taxa roles...')
  ros
}
