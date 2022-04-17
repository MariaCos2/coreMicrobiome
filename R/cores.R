# library(dplyr) #add_rownames
# library(tidyr)  #spread
get_heatmap_coordinates <- function(p, data, data1, offset=0, width=1) {
  
  variable <- value <- lab <- y <- NULL
  
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
  
  isTip <- x <- from <- to <- custom_labels <- NULL
  
  df <- p$data
  nodeCo <- intersect(df %>% dplyr::filter(is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist(), 
                      df %>% dplyr::filter(!is.na(x)) %>% 
                        select(.data$parent, .data$node) %>% unlist())
  labCo <- df %>% dplyr::filter(.data$node %in% nodeCo) %>% 
    select(.data$label) %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- df$label %in% selCo
  
  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset
  
  dd <- as.data.frame(data)
  dd1 <- as.data.frame(data1)
  i <- order(df$y)
  
  i <- i[!is.na(df$y[i])]
  
  lab <- df$label[i]
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
  dd1 <- dd1[match(lab, rownames(dd1)), , drop = FALSE]
  
  dd$y <- sort(df$y)
  dd$lab <- lab
  dd <- gather(dd, variable, value, -c(lab, y))
  dd1$y <- sort(df$y)
  dd1$lab <- lab
  dd1 <- gather(dd1, variable, value, -c(lab, y))
  
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
    dd1$value[i] <- NA
  }
  dd$variable <- factor(dd$variable, levels=colnames(data))
  dd1$variable <- factor(dd1$variable, levels=colnames(data))
  
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from=dd$variable, to=V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  
  V3 <- V2 + max(V2) - min(V2) + offset - 0.1
  mapping <- data.frame(from=dd1$variable, to=V3)
  mapping <- unique(mapping)
  dd1$x <- V3
  dd1$width <- width
  dd1[[".panel"]] <- factor("Tree")
  
  rbind(dd, dd1)
}

common_core <- function(otu, map, mini_abun=0, threshold=0.01, sample_name, sample_group) {
  
  dims = dim(otu)[1]
  if (dims > 1000){
    dims = 300
  }
  
  otu_PA <- 1*((otu>mini_abun) == 1)                                         # presence: otu>mini_abun
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                    # occupancy calculation
  otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)    # mean relative abundance
  occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu')     # combine occupancy and abundance data frame
  occ_abun$rank <- as.factor(occ_abun$otu)
  
  map$Sample_ID <- map$sample_name
  map$SampleType <- map$sample_group
  
  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
    gather(Sample_ID, abun, -otu) %>%
    left_join(map, by = 'Sample_ID') %>%
    group_by(otu, SampleType) %>%
    summarise(plot_freq = sum(abun > mini_abun)/length(abun),        # occurrence frequency in each sampletype
              coreSite = ifelse(plot_freq == 1, 1, 0), 
              detect = ifelse(plot_freq > 0, 1, 0)) %>%    
    group_by(otu) %>%
    summarise(sumF = sum(plot_freq),
              sumG = sum(coreSite),
              nS = length(SampleType)*2,
              Index = (sumF + sumG)/nS) 
  
  otu_ranked <- occ_abun %>%
    left_join(PresenceSum, by = 'otu') %>%
    transmute(otu = otu,
              rank = Index) %>%
    arrange(desc(rank))
  
  BCaddition <- NULL
  
  otu_start = otu_ranked$otu[1]                   
  start_matrix <- as.matrix(otu[otu_start,])
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[,x[2]]))/(sum(start_matrix[,x[1]] + start_matrix[,x[2]])))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
  
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1 
  BCaddition <- rbind(BCaddition,df_s) 
  
  for(i in 2:dims){                              
    otu_add <- otu_ranked$otu[i]                       
    add_matrix <- as.matrix(otu[otu_add,])
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]] - start_matrix[,x[2]]))/(sum(start_matrix[,x[1]] + start_matrix[,x[2]])))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse = ' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i 
    BCaddition <- left_join(BCaddition, df_a, by = c('x_names'))
  }
  
  rownames(BCaddition) <- BCaddition$x_names
  temp_BC <- BCaddition
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)
  
  BC_ranked <- data.frame(rank = as.factor(otu_ranked$otu[1:dims]), t(temp_BC_matrix)) %>% 
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC = mean(BC)) %>%
    arrange(-desc(MeanBC)) %>%
    mutate(proportionBC = MeanBC/max(MeanBC))
  Increase <- BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC = c(0,(Increase)))
  BC_ranked <- cbind(BC_ranked, increaseDF)
  
  lastCall <- as.numeric(which(BC_ranked$rank == last(BC_ranked$rank[(BC_ranked$IncreaseBC >= (1+threshold))]))) + 1
  if (length(lastCall)==0){(stop("Should reduce the threshold!!!"))}
  
  BC_ranked$dim <- c(1:dim(BC_ranked)[1])
  BC_ranked$fill <- 'non-core'
  BC_ranked$fill <- c(rep('core', times = lastCall), rep('non-core', times = dim(BC_ranked)[1] - lastCall))
  BC_ranked_abun <- BC_ranked %>% left_join(occ_abun, by = 'rank')
  
  list(otu_ranked = otu_ranked, BC_ranked_abun = BC_ranked_abun, lastCall = lastCall)
}

get_avg_rel_abun <- function(otu, map, otu_ranked, BC_ranked_abun, sample_name, sample_group) {
  
  map$Sample_ID <- map$sample_name
  map$SampleType <- map$sample_group
  
  otu_relabun <- decostand(otu, method = "total", MARGIN = 2)
  
  avg_rel_abun <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
    gather(Sample_ID, relabun, -otu) %>%
    left_join(map, by = 'Sample_ID') %>%
    left_join(otu_ranked, bu = 'otu') %>%
    group_by(SampleType, otu) %>%
    summarise(plot_freq = mean(relabun)) %>% 
    spread(key = "SampleType", value = "plot_freq") %>% 
    filter(otu %in% BC_ranked_abun$rank)
  
  names <- avg_rel_abun$otu
  avg_rel_abun = avg_rel_abun[, !grepl("otu", colnames(avg_rel_abun))]
  row.names(avg_rel_abun) <- names
  
  avg_rel_abun
}

get_avg_occ <- function(otu, map, otu_ranked, BC_ranked_abun, sample_name, sample_group) {
  
  map$Sample_ID <- map$sample_name
  map$SampleType <- map$sample_group
  
  otu_relabun <- decostand(otu, method = "total", MARGIN = 2)
  
  avg_occ <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
    gather(Sample_ID, relabun, -otu) %>%
    left_join(map, by = 'Sample_ID') %>%
    left_join(otu_ranked, bu = 'otu') %>%
    group_by(SampleType, otu) %>%
    summarise(plot_freq = sum(relabun>0)/length(relabun)) %>% 
    spread(key = "SampleType", value = "plot_freq") %>% 
    filter(otu %in% BC_ranked_abun$rank)
  
  names <- avg_occ$otu
  avg_occ = avg_occ[, !grepl("otu", colnames(avg_occ))]
  row.names(avg_occ) <- names
  
  avg_occ
}

dis_occ <- function(otu, map, mini_abun, threshold, sample_name, sample_group) {
  otu_PA <- 1*((otu>mini_abun) == 1)                                               
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                
  otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)     
  occ_abun <- add_rownames(as.data.frame(cbind(otu_occ, otu_rel)),'otu') 
  occ_abun$rank <- as.factor(occ_abun$otu)
  
  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun
  
  occ_abun_noncore <- occ_abun %>% filter(!rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  occ_abun_core <- occ_abun %>% filter(rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  occ <- abs(mean(occ_abun_core$otu_occ) - mean(occ_abun_noncore$otu_occ))
  
  occ
}

dis_occ_loop <- function(otu, sample, min_abun, max_abun, abun_step, min_thre, max_thre, thre_step) {
  
  dis_occ_df <- c()
  print('dd')
  for (i in seq(min_thre, max_thre, thre_step)){ 
    dis_occ <- c()
    print('cc')
    for (j in seq(min_abun, max_abun, abun_step)){
      print('ee')
      tr <- tryCatch({
        dis_occ(otu, sample, j, i)
      },error = function(e){
        0
      })
      dis_occ = append(dis_occ, tr)
    }
    dis_occ_df = cbind(dis_occ_df, dis_occ)
  }
  
  dis_occ_df = as.data.frame(dis_occ_df)
  colnames(dis_occ_df) = seq(min_thre, max_thre, thre_step)
  row.names(dis_occ_df) = seq(min_abun, max_abun, abun_step)
  dis_occ_df$name <- row.names(dis_occ_df)
  dis_occ_df = melt(dis_occ_df, id = c("name"))
  
  p <- ggplot(dis_occ_df, aes(name, variable)) +
    geom_tile(aes(fill = value, width = 0.4, height = 0.9)) +
    labs(fill = "mean difference") +
    scale_fill_gradient2(low = "white", high = 'red') +
    theme_minimal() +
    labs(title = 'Mean occurrence frequency difference', 
         subtitle = 'between core and non-core',
         x = 'Minimal abundance', y = 'Threshold',
         # caption = "*mean difference = 0 means NA"
    ) +
    theme(plot.title = element_text(size = 16, color = "black", hjust = 0, vjust = 1, lineheight = 0.2, face = "bold.italic"),
          axis.title.x = element_text(size = 15, color = "black", hjust = 0.5),
          axis.title.y = element_text(size = 15,color = "black", hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 14),
          text = element_text(size = 14),
          legend.position = "right",
          legend.title = element_text(colour = "black", size = 14),
          legend.text = element_text(colour ="black", size = 14),
          panel.background = element_rect(colour = "black", size = 1)
    )
  p
}

get_networks_coodinates <- function(p, x_r, cor_matrix) {
  
  f <- data.frame(row.names = as.data.frame(ggplot_build(p)$data[4])$label, val = as.data.frame(ggplot_build(p)$data[4])$y)
  
  xmin = c()
  xmax = c()
  ymin = c()
  ymax = c()
  cor = c()
  for (i in c(1:dim(cor_matrix)[1])) {
    ymin = append(ymin, f[cor_matrix[i,]$from,])
    ymax = append(ymax, f[cor_matrix[i,]$to,])
    xmin = append(xmin, x_r)
    xmax = append(xmax, x_r)
    cor = append(cor, cor_matrix[i,]$cor)
  }
  
  links <- data.frame(ymin = ymin,ymax = ymax,xmin = xmin,xmax = xmax,cor = cor)
  links$yma_ymi <- abs(links$ymax - links$ymin)
  
  links
}