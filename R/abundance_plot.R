##' @title abundance_plot
##'
##' @description Plots abundance among samples and groups (sample_group) at specific taxonomic level.
##' @param otu a community count data matrix with samples in rows and OTUs/taxa in column,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param rank a taxonomic rank (Phylum or Class or Order or Family or Genus or Species),
##' @param sample_group the name of a column that containing the sample type/group,
##' @param ...
##' @details Plots abundance among samples and groups (sample_group) at specific taxonomic level.
##' @return two subfigures containing abundance among samples and groups (sample_group) at specific taxonomic level.
##' @examples
##'  abundance_plot(otu, taxa, sample, sample_group, 'Genus')
##' @export

abundance_plot <- function(otu, taxa, sample, sample_group, rank){

  mydata = phyloseq(otu, taxa, sample)

  mydata_abund = phyloseq::transform_sample_counts(mydata, function(x){x / sum(x)})
  p1 <- phyloseq::plot_bar(mydata_abund, fill = rank) +
    geom_bar(aes(color = rank), stat = "identity", position = "stack") +
    labs(title = 'Relative Abundance', x = '', y = '') +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 11, colour ="black"),
          legend.title = element_text(size = 14, colour = "black"),
          plot.title = element_text(color = "black", size = 14, face = "bold.italic"))

  p2 <- phyloseq::plot_bar(mydata, x = sample_group, fill = rank) +
    geom_bar(aes(color = rank), stat = "identity", position = "stack") +
    labs(title = 'Abundance by Sample Type', x = '', y = '') +
    theme(panel.background = element_blank(),
          axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 11, colour ="black"),
          legend.title = element_text(size = 14, colour = "black"),
          plot.title = element_text(color = "black", size = 14, face = "bold.italic"))

  grid.arrange(p1, p2, nrow = 1)
  #ggplotly(p) %>% layout(height = hi, width = wi)
}
