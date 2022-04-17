# library(microbiome)
# library(phyloseq)
# library(pheatmap)

##' @title heatmap_plot
##' 
##' @description Plots heatmap based on count data, relative abundance data, clr transformation data or log10 transformation data.
##' @param otu a community count data matrix,
##' @param taxa a taxa information matrix. The rownames must match the OTU names (taxa_names) of the otu, 
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param rank a taxonomic rank (Phylum or Class or Order or Family or Genus or Species),
##' @param ...
##' @details Plots abundance among samples and groups (SampleType) at specific taxonomic level.
##' @return four heatmap with clustering abundance based on different transformation methods.
##' @examples
##'  heatmap_plot(otu, taxa, sample, 'Genus')
##' @export

heatmap_plot <- function(otu, taxa, sample, rank){
  # otu = read.csv(otu, sep = ",", header = T, row.names = 1)
  # taxa = read.csv(taxa, sep = ",", header = T, row.names = 1)
  # sample = read.csv(sample, sep = ",", header = T, row.names = 1)
  # 
  # OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
  # TAX = tax_table(as.matrix(taxa))
  # sample = sample_data(sample)
  mydata = phyloseq(otu, taxa, sample)
  
  mydata_r = aggregate_taxa(mydata, rank)
  otu = otu_table(mydata_r)
  p1 <- pheatmap(otu, cellwidth = 6, cellheight = 4, fontsize = 5, main = 'Count')
  
  mydata_r = aggregate_taxa(mydata, rank)
  mydata_r = transform(mydata_r, "compositional")
  otu = otu_table(mydata_r)
  p2 <- pheatmap(otu, cellwidth = 6, cellheight = 4, fontsize = 5, main = 'Relative abundance')
  
  mydata_r = aggregate_taxa(mydata, rank)
  mydata_r = transform(mydata_r, "clr")
  otu = otu_table(mydata_r)
  p3 <- pheatmap(otu, cellwidth = 6, cellheight = 4, fontsize = 5, main = 'clr transformation')
  
  mydata_r = aggregate_taxa(mydata, rank)
  mydata_r = transform(mydata_r, "log10")
  otu = otu_table(mydata_r)
  p4 <- pheatmap(otu, cellwidth = 6, cellheight = 4, fontsize = 5, main = 'log10 transformation')
  
  grid.arrange(p1[[4]], p2[[4]], p3[[4]], p4[[4]], nrow = 2)
  # plot_heatmap(mydata, method = "NMDS", distance = "bray", taxa.label = rank, low="#000033", high="#CCFF66") +
  #   theme(panel.background = element_blank(),
  #         axis.text.x = element_text(colour = "black", size = 12),
  #         axis.text.y = element_text(colour = "black", size = 5))
    }
