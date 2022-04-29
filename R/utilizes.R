filter_prevalence <- function(x, prevalence) {
  
  if ((prevalence < 0) | (prevalence > 1)) {
    stop("The prevalence argument should be in [0, 1].")
  }
  otu <- otu_table(x)
  otu_PA <- 1*((otu > 0) == 1)                             
  otu_pre <- data.frame(pre = rowSums(otu_PA)/ncol(otu_PA), taxa = row.names(otu))
  taxa <- otu_pre$taxa[otu_pre$pre > prevalence]
  prune_taxa(taxa, x)
}

delet_repeat <- function(x, cor) {
  delet_r <- c()
  for (i in c(1:(dim(x)[1]-1))){
    for (j in c((i+1):dim(x)[1])){
      if (x$cor[i] == x$cor[j]){
        delet_r = append(delet_r, j)
      }
    }
  }
  x <- x[-c(delet_r), ] 
  x
}

triu <- function(x) x[upper.tri(x)]

triu2diag <- function(x, diagval=0) {
  stopifnot(is.null(dim(x)), TRUE)
  e <- length(x)
  n <- .5 * (sqrt(8*e + 1)+1)
  mat <- matrix(0, n, n)
  mat[upper.tri(mat)] <- x
  mat <- mat + t(mat)
  diag(mat) <- diagval
  mat
}

norm_diric <- function(x, rep=1) {
  dmat <- VGAM::rdiric(rep, x+1)
  norm_to_total(colMeans(dmat))
}

norm_to_total <- function(x) x/sum(x)

cor2cov <- function(cor, sds) {
  if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
  cor * sds * rep(sds, each=nrow(cor))
}