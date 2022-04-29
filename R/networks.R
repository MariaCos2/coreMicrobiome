# library(propr)
# library(tidygraph)
# library(reshape2) # melt
# library(ggraph)
# library(igraph)
# library(qgraph) # smallworld
# library(compositions) # clr
# source("functions/utilizes.R")

#' This function runs SparCC in R
#' adapted from https://github.com/zdk123/SpiecEasi/blob/master/R/spaRcc.R
sparcc <- function(data, iter=20, inner_iter=10, th=.1) {

  sparccs <- lapply(1:iter, function(i)
    sparccinner(t(apply(data, 1, norm_diric)),
                iter=inner_iter, th=th))
  cors <- array(unlist(lapply(sparccs, function(x) x$Cor)),
                c(ncol(data),ncol(data),iter))
  corMed <- apply(cors, 1:2, median)
  covs <- array(unlist(lapply(sparccs, function(x) x$Cov)),
                c(ncol(data),ncol(data),iter))
  covMed <- apply(covs, 1:2, median)
  covMed <- cor2cov(corMed, sqrt(diag(covMed)))
  list(Cov=covMed, Cor=corMed)
}

#' This function gets bootstrapped estimates of SparCC correlation coefficients
#' adapted from https://github.com/zdk123/SpiecEasi/blob/master/R/spaRcc.R
sparccboot <- function(data, sparcc.params=list(),
                       statisticboot=function(data, indices) triu(do.call("sparcc",
                                                                          c(list(data[indices,,drop=FALSE]), sparcc.params))$Cor),
                       statisticperm=function(data, indices) triu(do.call("sparcc",  c(list(apply(data[indices,], 2, sample)), sparcc.params))$Cor),
                       R, ncpus=1, ...) {

  if (!requireNamespace('boot', quietly=TRUE))
    stop('\'boot\' package is not installed')

  res     <- boot::boot(data, statisticboot, R=R, parallel="multicore", ncpus=ncpus, ...)
  null_av <- boot::boot(data, statisticperm, sim='permutation', R=R, parallel="multicore", ncpus=ncpus)
  class(res) <- 'list'
  structure(c(res, list(null_av=null_av)), class='sparccboot')
}

#' This function get empirical p-values from bootstrap SparCC output
#' adapted from https://github.com/zdk123/SpiecEasi/blob/master/R/spaRcc.R
pval.sparccboot <- function(x, sided='both') {

  if (sided != "both") stop("only two-sided currently supported")
  nparams  <- ncol(x$t)
  tmeans   <- colMeans(x$null_av$t)
  niters   <- nrow(x$t)
  ind95    <- max(1,round(.025*niters)):round(.975*niters)
  boot_ord <- apply(x$t, 2, sort)
  boot_ord95 <- boot_ord[ind95,]
  outofrange <- unlist(lapply(1:length(x$t0), function(i) {
    aitvar <- x$t0[i]
    range  <- range(boot_ord95[,i])
    range[1] > aitvar || range[2] < aitvar
  }))

  bs_above <- unlist(lapply(1:nparams, function(i)
    length(which(x$t[, i] > tmeans[i]))))
  is_above <- bs_above > x$R/2
  cors <- x$t0

  pvals    <- ifelse(is_above, 2*(1-bs_above/x$R), 2*bs_above/x$R)
  pvals[pvals > 1]  <- 1
  pvals[outofrange] <- NaN
  list(cors=cors, pvals=pvals)
}

#' @noRd
sparccinner <- function(data.f, T=NULL, iter=10, th=0.1) {

  if (is.null(T))   T  <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor

  excluded <- NULL
  for (i in 1:iter) {
    res.excl <- exclude_pairs(Cor, M, th, excluded)
    M <- res.excl$M
    excluded <- res.excl$excluded
    if (res.excl$break_flag) break
    res.bv <- basis_var(T, M=M, excluded=excluded)
    Vbase  <- res.bv$Vbase
    M      <- res.bv$M
    K <- M
    diag(K) <- 1
    cbase  <- C_from_V(T, Vbase)
    Cov    <- cbase$Cov
    Cor    <- cbase$Cor
  }
  list(Cov=Cov, Cor=Cor, i=i, M=M, excluded=excluded)
}

#' @noRd
exclude_pairs <- function(Cor, M, th=0.1, excluded=NULL) {

  break_flag <- FALSE
  C_temp <- abs(Cor - diag(diag(Cor)) )
  if (!is.null(excluded)) C_temp[excluded] <- 0
  exclude <- which(abs(C_temp - max(C_temp)) < .Machine$double.eps*100)[1:2]
  if (max(C_temp) > th)  {
    i <- na.exclude(arrayInd(exclude, c(nrow(M), ncol(M)))[,1])
    M[i,i] <- M[i,i] - 1
    excluded_new <- c(excluded, exclude)
  } else {
    excluded_new <- excluded
    break_flag   <- TRUE
  }
  list(M=M, excluded=excluded_new, break_flag=break_flag)
}

#' @noRd
basis_cov <- function(data.f) {

  T <- av(data.f)
  res.bv <- basis_var(T)
  Vbase  <- res.bv$Vbase
  M      <- res.bv$M
  cbase  <- C_from_V(T, Vbase)
  Cov    <- cbase$Cov
  Cor    <- cbase$Cor
  list(Cov=Cov, M=M)
}

#' @noRd
basis_var <- function(T, CovMat = matrix(0, nrow(T), ncol(T)),
                      M = matrix(1, nrow(T), ncol(T)) + (diag(ncol(T))*(ncol(T)-2)),
                      excluded = NULL, Vmin=1e-4) {

  if (!is.null(excluded)) {
    T[excluded] <- 0
  }
  Ti     <- matrix(rowSums(T))
  CovVec <- matrix(rowSums(CovMat - diag(diag(CovMat)))) # row sum of off diagonals
  M.I <- tryCatch(solve(M), error=function(e) MASS::ginv(M))
  Vbase <- M.I %*% (Ti + 2*CovVec)
  Vbase[Vbase < Vmin] <- Vmin
  list(Vbase=Vbase, M=M)
}

#' @noRd
C_from_V <- function(T, Vbase) {

  J      <- matrix(1, nrow(T), ncol(T))
  Vdiag  <- diag(c(Vbase))
  CovMat <- .5*((J %*% Vdiag) + (Vdiag %*% J) - T)
  CovMat <- (CovMat + t(CovMat))/2

  CorMat <- cov2cor(CovMat)
  CorMat[abs(CorMat) > 1] <- sign(CorMat[abs(CorMat) > 1])
  CovMat <- cor2cov(CorMat, sqrt(as.vector(Vbase)))
  list(Cov=CovMat, Cor=CorMat)
}

#' @noRd
av <- function(data) {

  cov.clr <- cov(clr(data))
  J <- matrix(1, ncol(data), ncol(data))
  (J %*% diag(diag(cov.clr))) + (diag(diag(cov.clr)) %*% J) - (2*cov.clr)
}


#' This function runs CCLasso in R
#' adapted from https://github.com/stefpeschel/NetCoMi/blob/master/R/cclasso.R
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3,
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {

  n <- nrow(x);
  p <- ncol(x);

  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);

  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);

  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) /
                   (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;

  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1];
  b1 <- lam_int2[2];

  lams <- NULL;
  fvals <- NULL;

  a2 <- a1 + 0.382 * (b1 - a1);
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv,
                         sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
                         wd2 = wd2);
  lams <- c(lams, b2);
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv,
                         sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
                         wd2 = wd2);
  lams <- c(lams, a2);
  fvals <- c(fvals, fa2$cv_loss);

  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;

  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);

    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv,
                             sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
                             wd2 = wd2);

      lams <- c(lams, b2);
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv,
                             sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd,
                             wd2 = wd2);

      lams <- c(lams, a2);
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);

    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2,
                  lam_int = 10^c(a1, b1));
  if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",",
        10^b1, "])\n", sep = "");
  }

  lambda <- 10^((a2 + b2)/2);

  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2,
                            n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w,
              p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}

#' This function gets bootstrapped estimates of CCLasso correlation coefficients
#' adapted from https://github.com/stefpeschel/NetCoMi/blob/master/R/cclasso.R
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma,
                            wd, u_f, u0_wd, d0_wd, wd2) {

  n <- nrow(x);
  p <- ncol(x);

  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);

    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2,
                         wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }

  return(list(cv_loss = cv_loss, sigma = sigma));
}

#' This function runs CCLasso for only one lambda
#' adapted from https://github.com/stefpeschel/NetCoMi/blob/master/R/cclasso.R
cclasso_sub <- function(sigma, vx, lambda2,
                        wd, u_f, u0_wd, d0_wd,
                        k_max = 200, x_tol = 1e-4) {

  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;

  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {

    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) *
                                   d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);

    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) *
      (A + lambda2);

    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);

    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1),
                abs(sigma2_new - sigma2)/(abs(sigma2) + 1) );
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }

  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max,
        "&& Relative error:", err, "!\n");
  }

  return(sigma);
}

#' This function runs CCLasso for only one lambda
#' adapted from https://github.com/stefpeschel/NetCoMi/blob/master/R/cclasso.R
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20,
                         wd, u_f, u0_wd, d0_wd) {

  n <- nrow(x);
  p <- ncol(x);

  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);

  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T),
                     ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]),
                          lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);

  vars2 <- rowMeans(vars_boot);

  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art),
                           lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];

  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) *
                                                           (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);

  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) +
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;

  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));
}

networks_statistics <- function(g, group) {

  names=V(g)$name
  de = igraph::degree(g)
  be = igraph::betweenness(g, normalized = T)
  cl = igraph::closeness(g)
  ei = igraph::eigen_centrality(g)$vector
  hub = igraph::hub_score(g)$vector
  nodes = data.frame(node.name = names, degree = de, betweenness = be, closeness = cl, eigenvector = ei, hub = hub, group = group)

  d <- degree_distribution(g, mode = "all", cumulative = T)
  scalefree = fit_power_law(d)$KS.p
  smallworld = smallworldIndex(g)$index

  density = edge_density(g)
  mean_dis = mean_distance(g)
  trans = transitivity(g, "global")
  trans_local = transitivity(g, "localaverage")

  stat = list(nodes = nodes, scalefree = scalefree, smallworld = smallworld, density = density, mean_dis = mean_dis, trans = trans, trans_local = trans_local)
}

get_propr_statistics <- function(otu, fdr_threshold, cor_threshold, permutation, propr_col) {

  pr <- propr(t(otu), metric = "rho", ivar = "clr", alpha = NA, p = permutation)
  cutoffs <- updateCutoffs(pr,
                           cutoff = seq(0, 1, .05),
                           ncores = 1)
  fdrs <- cutoffs@fdr
  fdr1 <- 0

  for (i in c(1:(dim(fdrs)[1]-1))){
    if (fdrs$FDR[i] < fdr_threshold && fdr1 > fdr_threshold){
      ccutoff <- fdrs$cutoff[i]
    }
    fdr1 <- fdrs$FDR[i]
  }

  propr.results <- getResults(pr) %>%
    filter(propr > ccutoff) %>%
    dplyr::select(Partner, Pair, propr)
  colnames(propr.results) <- c('from', 'to', 'cor')

  if (dim(propr.results)[1] == 0){
    stop("Propr cannot generate any co-occurrence. Please raise the fdr threshold or/and reduce cor threshold!!!!")
  }

  propr.n <- propr.results %>% filter(cor > cor_threshold)
  propr.g = graph.data.frame(propr.n, directed = FALSE)
  propr.f = networks_statistics(propr.g, group = 'propr')

  propr.n <- propr.results %>% filter(cor > cor_threshold)
  propr.nn = propr.n %>% dplyr::select(from = from, to = to)
  graph_gt <- as_tbl_graph(propr.nn) %>%
    dplyr::mutate(betweenness = centrality_betweenness(), degree = centrality_degree(mode = 'in'))
  p <- ggraph(graph_gt,layout = 'fr') +
    geom_edge_link(color = "grey")+
    geom_node_point(aes(size = degree, alpha = degree), color = propr_col) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    theme_graph() +
    ggtitle("Propr") +
    theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"),
          panel.border = element_rect(fill = NA,color = "black", size = 0.5, linetype = "solid"))

  list(n = propr.n, f = propr.f, g = propr.g, p = p)

}

get_sparcc_statistics <- function(otu, fdr_threshold, cor_threshold, permutation, sparcc_col) {

  sparcc.boot = sparccboot(t(otu), R = permutation)
  sparcc.p = pval.sparccboot(sparcc.boot)

  sparcc.stat <- as.data.frame(triu2diag(sparcc.p$pvals))
  row.names(sparcc.stat) <- row.names(otu)
  colnames(sparcc.stat) <- row.names(otu)
  sparcc.stat$taxa <- row.names(sparcc.stat)
  sparcc.stat.l <- melt(sparcc.stat, id = c("taxa"))
  sparcc.stat.l <- sparcc.stat.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  sparcc.stat.l$variable <- as.character(sparcc.stat.l$variable)
  colnames(sparcc.stat.l) <- c('p', 'from', 'to')

  sparcc.cor <- as.data.frame(triu2diag(sparcc.p$cors))
  row.names(sparcc.cor) <- row.names(otu)
  colnames(sparcc.cor) <- row.names(otu)
  sparcc.cor$taxa <- row.names(sparcc.cor)
  sparcc.cor.l <- melt(sparcc.cor, id = c("taxa"))
  sparcc.cor.l <- sparcc.cor.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  sparcc.cor.l$variable <- as.character(sparcc.cor.l$variable)
  colnames(sparcc.cor.l) <- c('cor', 'from', 'to')

  sparcc.results <- cbind(sparcc.stat.l, cor = sparcc.cor.l$cor) %>%
    mutate(fdr = p.adjust(p, method = "BH"))

  sparcc.n <- sparcc.results %>% filter(fdr < fdr_threshold) %>% filter(cor > cor_threshold)
  if (dim(sparcc.n)[1] == 0){
    stop("Sparcc cannot generate any co-occurrence. Please raise the fdr threshold or/and reduce cor threshold!!!!")
  }

  sparcc.n <- delet_repeat(sparcc.n, 'cor')
  sparcc.n <- sparcc.n %>% dplyr::select(from, to, cor)
  sparcc.g = graph.data.frame(sparcc.n, directed = FALSE)
  sparcc.f = networks_statistics(sparcc.g, group = 'sparcc')

  sparcc.nn = sparcc.n %>% dplyr::select(from = from, to = to)
  graph_gt <- as_tbl_graph(sparcc.nn) %>%
    dplyr::mutate(betweenness = centrality_betweenness(), degree = centrality_degree(mode = 'in'))
  p <- ggraph(graph_gt,layout = 'fr') +
    geom_edge_link(color = "grey")+
    geom_node_point(aes(size = degree, alpha = degree), color = sparcc_col) +
    geom_node_text(aes(label = name),size = 3, repel = TRUE) +
    theme_graph() +
    ggtitle("Sparcc") +
    theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"),
          panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

  list(n = sparcc.n, f = sparcc.f, g = sparcc.g, p = p)
}

get_cclasso_statistics <- function(otu, fdr_threshold, cor_threshold, permutation, cclasso_col) {

  cclasso.results <- cclasso(t(otu), counts = T, n_boot = permutation)

  cclasso.stat <- as.data.frame(cclasso.results$p_vals)
  row.names(cclasso.stat) <- row.names(otu)
  colnames(cclasso.stat) <- row.names(otu)
  cclasso.stat$taxa <- row.names(cclasso.stat)
  cclasso.stat.l <- melt(cclasso.stat, id = c("taxa"))
  cclasso.stat.l <- cclasso.stat.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  cclasso.stat.l$variable <- as.character(cclasso.stat.l$variable)
  colnames(cclasso.stat.l) <- c('p', 'from', 'to')

  cclasso.cor <- as.data.frame(cclasso.results$cor_w)
  row.names(cclasso.cor) <- row.names(otu)
  colnames(cclasso.cor) <- row.names(otu)
  cclasso.cor$taxa <- row.names(cclasso.cor)
  cclasso.cor.l <- melt(cclasso.cor, id = c("taxa"))
  cclasso.cor.l <- cclasso.cor.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  cclasso.cor.l$variable <- as.character(cclasso.cor.l$variable)
  colnames(cclasso.cor.l) <- c('cor', 'from', 'to')

  cclasso.results <- cbind(cclasso.stat.l, cor = cclasso.cor.l$cor) %>%
    mutate(fdr = p.adjust(p, method = "BH"))

  fdr_threshold = 1
  cor_threshold = 0.2
  cclasso.n <- cclasso.results  %>% filter(fdr < fdr_threshold) %>% filter(cor > cor_threshold)
  if (dim(cclasso.n)[1] == 0){
    stop("CClasso cannot generate any co-occurrence. Please raise the fdr threshold or/and reduce cor threshold!!!!")
  }

  cclasso.n <- delet_repeat(cclasso.n, 'cor')
  cclasso.n <- cclasso.n %>% dplyr::select(from, to, cor)
  cclasso.g = graph.data.frame(cclasso.n, directed = FALSE)
  cclasso.f = networks_statistics(cclasso.g, group = 'cclasso')

  cclasso.nn = cclasso.n %>% dplyr::select(from = from, to = to)
  graph_gt <- as_tbl_graph(cclasso.nn) %>%
    dplyr::mutate(betweenness = centrality_betweenness(), degree = centrality_degree(mode = 'in'))
  p <- ggraph(graph_gt,layout = 'fr') +
    geom_edge_link(color = "grey")+
    geom_node_point(aes(size = degree, alpha = degree), color = cclasso_col) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE) +
    theme_graph() +
    ggtitle("CClasso") +
    theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"),
          panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

  list(n = cclasso.n, f = cclasso.f, g = cclasso.g, p = p)

}

get_spieceasi_statistics <- function(otu, weight_threshold, sampling, spieceasi_col) {

  dat = as.matrix(t(otu))
  se.mb.amgut <- spiec.easi(dat, method='mb', lambda.min.ratio=1e-2,
                            nlambda=20, pulsar.params=list(rep.num=sampling))

  spieceasi.cor <- as.data.frame(as.matrix(symBeta(getOptBeta(se.mb.amgut), mode = "maxabs")))
  row.names(spieceasi.cor) <- row.names(otu)
  colnames(spieceasi.cor) <- row.names(otu)
  spieceasi.cor$taxa <- row.names(spieceasi.cor)
  spieceasi.cor.l <- melt(spieceasi.cor, id = c("taxa"))
  spieceasi.cor.l <- spieceasi.cor.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  spieceasi.cor.l$variable <- as.character(spieceasi.cor.l$variable)
  colnames(spieceasi.cor.l) <- c('cor', 'from', 'to')

  spieceasi.net <- as.data.frame(as.matrix(getRefit(se.mb.amgut)))
  row.names(spieceasi.net) <- row.names(otu)
  colnames(spieceasi.net) <- row.names(otu)
  spieceasi.net$taxa <- row.names(spieceasi.net)
  spieceasi.net.l <- melt(spieceasi.net, id = c("taxa"))
  spieceasi.net.l <- spieceasi.net.l %>% filter(taxa != variable) %>% select(value, taxa, variable)
  spieceasi.net.l$variable <- as.character(spieceasi.net.l$variable)
  colnames(spieceasi.net.l) <- c('n', 'from', 'to')

  spieceasi.results <- cbind(spieceasi.cor.l, n = spieceasi.net.l$n)
  spieceasi.n <- spieceasi.results %>% filter(n == 1) %>% filter(cor > weight_threshold)
  spieceasi.n <- spieceasi.n %>% dplyr::select(from, to, cor)
  spieceasi.n <- delet_repeat(spieceasi.n, 'cor')

  if (dim(spieceasi.n)[1] == 0){
    stop("spieceasi cannot generate any co-occurrence. Please reduce the weight threshold!!!!")
  }

  spieceasi.g = graph.data.frame(spieceasi.n, directed = FALSE)
  spieceasi.f = networks_statistics(spieceasi.g, group = 'spieceasi')

  spieceasi.nn = spieceasi.n %>% dplyr::select(from = from, to = to)
  graph_gt <- as_tbl_graph(spieceasi.nn) %>%
    dplyr::mutate(betweenness = centrality_betweenness(), degree = centrality_degree(mode = 'in'))
  p <- ggraph(graph_gt,layout = 'fr') +
    geom_edge_link(color = "grey")+
    geom_node_point(aes(size = degree, alpha = degree), color = spieceasi_col) +
    geom_node_text(aes(label = name),size = 3, repel = TRUE) +
    theme_graph() +
    ggtitle("spieceasi") +
    theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic"),
          panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"))

  list(n = spieceasi.n, f = spieceasi.f, g = spieceasi.g, p = p)

}

ZiPivalue = function(igraph = igraph,method = "cluster_fast_greedy"){

  if (method == "cluster_walktrap" ) {
    fc <- igraph::cluster_walktrap(igraph) # ,weights =  abs(igraph::E(igraph)$weight)
    print(fc)
  }
  if (method == "cluster_edge_betweenness" ) {
    print('yes')
    fc <- igraph::cluster_edge_betweenness(igraph) # ,weights =  abs(igraph::E(igraph)$weight)
  }
  if (method == "cluster_fast_greedy" ) {
    fc <- igraph::cluster_fast_greedy(igraph) # ,weights =  abs(igraph::E(igraph)$weight)
  }
  if (method == "cluster_spinglass" ) {
    fc <- igraph::cluster_spinglass(igraph) # ,weights =  abs(igraph::E(igraph)$weight)
  }

  modularity <- igraph::modularity(igraph,igraph::membership(fc))

  comps <- igraph::membership(fc)

  igraph::V(igraph)$module <- as.character(comps)

  taxa.roles <- module.roles(igraph)

  taxa.roles$label = row.names(taxa.roles)
  for (i in 1:nrow(taxa.roles))if(taxa.roles[i,3]> 0.62|taxa.roles[i,1]> 2.5) {
    taxa.roles[i,5]=taxa.roles[i,5]
  }else{
    taxa.roles[i,5]= ""
  }
  taxa.roles$role_7 = taxa.roles$roles

  taxa.roles <- na.omit(taxa.roles)   # remove NA values
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Peripherals'
  taxa.roles[which(taxa.roles$z < 2.5 & taxa.roles$p > 0.62),'roles'] <- 'Connectors'
  taxa.roles[which(taxa.roles$z > 2.5 & taxa.roles$p < 0.62),'roles'] <- 'Module hubs'
  taxa.roles[which(taxa.roles$z > 2.5 & taxa.roles$p > 0.62),'roles'] <- 'Network hubs'

  #return(list(p, taxa.roles))
  taxa.roles
}

module.roles <- function(comm_graph){

  td <- network_degree(comm_graph)
  wmd <- within_module_degree(comm_graph)
  z <- zscore(wmd)
  amd <- among_module_connectivity(comm_graph)
  pc <- participation_coeffiecient(amd, td)
  zp <- data.frame(z,pc)
  nod.roles <- assign_module_roles(zp)
  return(nod.roles)
}

network_degree <- function(comm_graph){

  ki_total <-NULL
  net_degree <- igraph::degree(comm_graph)
  for(i in 1:length(V(comm_graph))){
    ki <- net_degree[i]
    tmp <- data.frame(taxa=names(ki), total_links=ki)
    if(is.null(ki_total)){ki_total<-tmp} else{ki_total <- rbind(ki_total, tmp)}
  }
  return(ki_total)
}

within_module_degree <- function(comm_graph){

  mods <- get.vertex.attribute(comm_graph, "module")
  vs <- as.list(V(comm_graph))
  modvs <- data.frame("taxon"= names(vs), "mod"=mods)
  sg1 <- decompose.graph(comm_graph,mode="strong")
  df <- data.frame()
  for(mod in unique(modvs$mod)){
    mod_nodes <- subset(modvs$taxon,modvs$mod==mod)
    neighverts <- unique(unlist(sapply(sg1,FUN=function(s){if(any(V(s)$name %in% mod_nodes)) V(s)$name else NULL})))
    g3 <- induced.subgraph(graph=comm_graph,vids=neighverts)
    mod_degree <- igraph::degree(g3)
    for(i in mod_nodes){
      ki <- mod_degree[which(names(mod_degree)==i)]
      tmp <- data.frame(module=mod, taxa=names(ki), mod_links=ki)
      df <- rbind(df,tmp)
    }
  }
  return(df)
}

zscore <- function(mod.degree){

  ksi_bar <- aggregate(mod_links ~ module, data=mod.degree, FUN = mean)
  ksi_sigma <- aggregate(mod_links ~ module, data=mod.degree, FUN = sd)
  z <- NULL
  for(i in 1:dim(mod.degree)[1]){
    mod_mean <- ksi_bar$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    mod_sig <- ksi_sigma$mod_links[which(ksi_bar$module == mod.degree$module[i])]
    z[i] <- (mod.degree$mod_links[i] - mod_mean)/mod_sig
  }
  z <- data.frame(row.names=rownames(mod.degree), z, module=mod.degree$module)
  return(z)
}

among_module_connectivity <- function(comm_graph){

  mods <- get.vertex.attribute(comm_graph, "module")
  vs <- as.list(V(comm_graph))
  modvs <- data.frame("taxa"= names(vs), "mod"=mods)
  df <- data.frame()
  for(i in modvs$taxa){
    for(j in modvs$taxa){
      if(are_adjacent(graph=comm_graph, v1=i , v2=j)){
        mod <- subset(modvs$mod, modvs$taxa==j)
        tmp <- data.frame(taxa=i, taxa2=j, deg=1, mod_links=mod)
        df <- rbind(df, tmp)
      }
    }
  }
  out <- aggregate(list(mod_links=df$deg), by=list(taxa=df$taxa, module=df$mod), FUN=sum)
  return(out)
}

participation_coeffiecient <- function(mod.degree, total.degree){

  p <- NULL
  for(i in total.degree$taxa){
    ki <- subset(total.degree$total_links, total.degree$taxa==i)
    taxa.mod.degree <- subset(mod.degree$mod_links, mod.degree$taxa==i)
    p[i] <- 1 - (sum((taxa.mod.degree)**2)/ki**2)
  }
  p <- as.data.frame(p)
  return(p)
}

assign_module_roles <- function(zp){

  zp <- na.omit(zp)
  zp$roles <- rep(0, dim(zp)[1])
  outdf <- NULL
  for(i in 1:dim(zp)[1]){
    df <- zp[i, ]
    if(df$z < 2.5){ #non hubs
      if(df$p < 0.05){
        df$roles <- "ultra peripheral"
      }
      else if(df$p < 0.620){
        df$roles <- "peripheral"
      }
      else if(df$p < 0.80){
        df$roles <- "non hub connector"
      }
      else{
        df$roles <- "non hub kinless"
      }
    }
    else { # module hubs
      if(df$p < 0.3){
        df$roles <- "provincial hub"
      }
      else if(df$p < 0.75){
        df$roles <- "connector hub"
      }
      else {
        df$roles <- "kinless hub"
      }
    }
    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
  }
  return(outdf)
}

