##' @title random_forest_plot
##'
##' @description Plots the importance of otu on predicting/classifying the sample_group and the accuracy using core otu, non-core otu and all otu.
##' @param otu a community count data matrix,
##' @param sample a sample information data.frame. The rownames must match the sample names in the otu,
##' @param mini_abun a value indicating whether the otu is present,
##' @param threshold a value indicating the final percent increase in beta-diversity,
##' @param sample_name the name of a column that containing the sample name,
##' @param sample_group the name of a column that containing the sample type/group,
##' @param ...
##' @details Plots the importance of otu on predicting/classifying the sample_group and the accuracy using core otu, non-core otu and all otu.
##' @return A bar plot showing importance of otu on predicting/classifying the sample_group and a bar plot showing the accuracy using core otu, non-core otu and all otu.
##' @examples
##'  random_forest_plot(otu, sample, mini_abun=0, threshold=0.02, sample_name, sample_group)
##' @export
##'
random_forest_plot <- function(otu, sample, mini_abun, threshold, sample_name, sample_group){

  BC_ranked <- common_core(otu, sample, mini_abun = mini_abun, threshold = threshold, sample_name, sample_group)
  BC_ranked_abun <- BC_ranked$BC_ranked_abun

  otu_rf <- data.frame(otu)
  otu_rf$rank <- row.names(otu)

  otu_rf <- otu_rf %>% dplyr::filter(rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  row.names(otu_rf) <- otu_rf$rank
  otu_rf <- otu_rf[, !grepl("rank", colnames(otu_rf))]
  otu_rf <- as.data.frame(t(otu_rf))
  otu_rf$Sample_ID <- row.names(otu_rf)

  sample <- data.frame(sample)
  sample$Sample_ID <- data.frame(sample[, which(sample_name == colnames(sample))])[, 1]
  sample$SampleType <- data.frame(sample[, which(sample_group == colnames(sample))])[, 1]
  map_sub <- sample %>% dplyr::select(Sample_ID, SampleType)

  data_rf <- otu_rf %>% left_join(map_sub, by = 'Sample_ID')
  data_rf <- data_rf[, !grepl("Sample_ID", colnames(data_rf))]
  data_rf$SampleType <- as.factor(data_rf$SampleType)
  row.names(data_rf) <- 1:nrow(data_rf)

  set.seed(222)
  ind <- sample(2, nrow(data_rf), replace = TRUE, prob = c(0.7, 0.3))
  train <- data_rf[ind == 1,]
  test <- data_rf[ind == 2,]

  rf <- randomForest(y=train$SampleType, x=train[,1:ncol(train)-1], proximity = TRUE)
  d_rf <- randomForest::importance(rf)
  importance_rf <- as.data.frame(d_rf) %>% arrange(desc(MeanDecreaseGini))
  importance_rf$id <- as.factor(row.names(importance_rf))
  p_core_train <- predict(rf, train[,1:ncol(train)-1])
  p_core_train_acc <- confusionMatrix(p_core_train, train$SampleType)$overall[1]
  p_core_test <- predict(rf, test[,1:ncol(test)-1])
  p_core_test_acc <-  confusionMatrix(p_core_test, test$SampleType)$overall[1]

  otu_rf <- otu
  otu_rf <- data.frame(otu)
  otu_rf$rank <- row.names(otu)
  otu_rf <- otu_rf %>% dplyr::filter(!rank %in% BC_ranked_abun$otu[BC_ranked_abun$fill == 'core'])
  row.names(otu_rf) <- otu_rf$rank
  otu_rf <- otu_rf[, !grepl("rank", colnames(otu_rf))]
  otu_rf <- as.data.frame(t(otu_rf))
  otu_rf$Sample_ID <- row.names(otu_rf)
  map_sub <- sample %>% dplyr::select(Sample_ID, SampleType)
  data_rf <- otu_rf %>% left_join(map_sub, by = 'Sample_ID')
  data_rf <- data_rf[, !grepl("Sample_ID", colnames(data_rf))]
  data_rf$SampleType <- as.factor(data_rf$SampleType)
  row.names(data_rf) <- 1:nrow(data_rf)
  set.seed(222)
  ind <- sample(2, nrow(data_rf), replace = TRUE, prob = c(0.7, 0.3))
  train <- data_rf[ind == 1,]
  test <- data_rf[ind == 2,]
  rf <- randomForest(y=train$SampleType, x=train[,1:ncol(train)-1], proximity = TRUE)
  p_noncore_train <- predict(rf, train[,1:ncol(train)-1])
  p_noncore_train_acc <- confusionMatrix(p_noncore_train, train$SampleType)$overall[1]
  p_noncore_test <- predict(rf, test[,1:ncol(test)-1])
  p_noncore_test_acc <-  confusionMatrix(p_noncore_test, test$SampleType)$overall[1]

  otu_rf <- otu
  otu_rf <- data.frame(otu)
  otu_rf$rank <- row.names(otu)
  row.names(otu_rf) <- otu_rf$rank
  otu_rf <- otu_rf[, !grepl("rank", colnames(otu_rf))]
  otu_rf <- as.data.frame(t(otu_rf))
  otu_rf$Sample_ID <- row.names(otu_rf)
  map_sub <- sample %>% dplyr::select(Sample_ID, SampleType)
  data_rf <- otu_rf %>% left_join(map_sub, by = 'Sample_ID')
  data_rf <- data_rf[, !grepl("Sample_ID", colnames(data_rf))]
  data_rf$SampleType <- as.factor(data_rf$SampleType)
  row.names(data_rf) <- 1:nrow(data_rf)
  set.seed(222)
  ind <- sample(2, nrow(data_rf), replace = TRUE, prob = c(0.7, 0.3))
  train <- data_rf[ind == 1,]
  test <- data_rf[ind == 2,]
  rf <- randomForest(y=train$SampleType, x=train[,1:ncol(train)-1], proximity = TRUE)
  p_all_train <- predict(rf, train[,1:ncol(train)-1])
  p_all_train_acc <- confusionMatrix(p_all_train, train$SampleType)$overall[1]
  p_all_test <- predict(rf, test[,1:ncol(test)-1])
  p_all_test_acc <-  confusionMatrix(p_all_test, test$SampleType)$overall[1]

  print('importance_rf')
  print(importance_rf)
  p1 <- ggbarplot(importance_rf, x = "id", y = "MeanDecreaseGini",
                  fill = 'dodgerblue3',
                  color = "black",
                  sort.val = "desc",
                  sort.by.groups = FALSE,
                  x.text.angle = 90,
                  font.label = list(color = "white", size = 12,
                                    vjust = 0.5),
                  title = 'Importance (Mean Decrease Gini) of otu',
                  ylab = FALSE,
                  xlab = "OTU",
  ) +
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
    )

  acc_rf <- data.frame (
    type  = c("core", "core", "non-core", "non-core", "all", "all"),
    train_test = c("train", "test", "train", "test", "train", "test"),
    acc = c(p_core_train_acc, p_core_test_acc, p_noncore_train_acc, p_noncore_test_acc, p_all_train_acc, p_all_test_acc),
    id = c('core_train_acc', 'core_test_acc', 'noncore_train_acc', 'noncore_test_acc', 'all_train_acc', 'all_test_acc')
  )
  acc_rf$id <- factor(acc_rf$id, levels = c('core_train_acc', 'core_test_acc', 'noncore_train_acc', 'noncore_test_acc', 'all_train_acc', 'all_test_acc'))
  print('accuracy')
  print(acc_rf)
  p2 <- ggbarplot(acc_rf, x = "id", y = "acc",
                  fill = 'type',
                  color = "black",
                  palette = c("purple", "deeppink4", "darkseagreen3"),
                  sort.by.groups = FALSE,
                  x.text.angle = 90,
                  font.label = list(color = "white", size = 12,
                                    vjust = 0.5),
                  title = 'Accuracy of prediction',
                  ylab = FALSE,
                  xlab = FALSE,
  ) +
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
    )


  ggarrange(p1, p2, ncol = 2, nrow = 1)
  #p <- subplot(p1, p2, nrows = 2)

  #ggplotly(p, height = height, width = weight, showlegend = T)

  #list(p1, p2)

}

