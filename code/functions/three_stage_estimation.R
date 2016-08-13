###############################################################################
# Multi-resolution blockmodel
#
# file: three_stage_estimation.R
#
# This file contains functions to perform three-stage estimation of the model
# First stage: SBM.
# Second stage: Latent space VB for each block
# Third stage: MVN MLE of latent space parameters
#
# We use this estimation procedure to initialize our MCMC
#
# Author: tedwestling
###############################################################################


#############################################################
# Function: three_stage_est
#
# Performs the three state estimation 
#############################################################

three_stage_est <- function(Y, D, K=NULL, plot_folder=NULL,  ls_method='vb', sigma=NULL, nstart=2) { #tau=mean(colSums(Y)), normalize=TRUE, threshold=TRUE, gamma=ifelse(threshold, 1,0),
  N <- nrow(Y)
  miss_inds <- which(is.na(Y), arr.ind = TRUE)
  miss_inds <- miss_inds[miss_inds[,1] < miss_inds[,2],]
  # 
  # #impute missing values using observed node-wise edge frequencies for initialization
  # obs_prob <- colMeans(Y, na.rm=TRUE)
  # exp_prob <- outer(obs_prob, obs_prob)
  #if(sum(is.na(Y)) > 0) Y[is.na(Y)] <- exp_prob[is.na(Y)]
  if(is.null(K)) {
    all_clusts <- spectral_cluster(Y, Krange=(2:floor(nrow(Y)/10)), tau=tau, normalize=normalize, threshold=threshold, gamma=gamma, nstart = nstart, degree_correct = FALSE, assortative=TRUE)
    aics <- unlist(lapply(all_clusts, function(c) c$sbm_aic))
    clust <- all_clusts[[which.min(aics)]]  
  } else {
    clust <- spectral_cluster(Y, Krange=K, degree_correct = FALSE, assortative=TRUE)#, tau=tau, normalize=normalize, threshold=threshold, gamma=gamma, nstart=nstart)
    clust <- clust[[1]]
  }
  if(sum(is.na(Y)) > 0) Y[is.na(Y)] <- clust$pred[is.na(Y)]
  memb <- clust$clusters
  if(!is.null(plot_folder)) {
    png(paste0(plot_folder, "/spectral_fit.png"), height=6, width=6, units='in', res=300)
    plot_blocked_matrix(Y, memb)
    dev.off()
  }
  out <- list(K=clust$K,
              gamma=memb,
              B=clust$Bhat)
  K <- length(unique(memb))
  
  Z <- matrix(nrow=N, ncol=D)
  theta <- matrix(NA, nrow=K, ncol=2)
  beta <- rep(NA, K)
  sigma <- rep(NA, K)
  latent_spaces <- for(k in 1:K) {
    k_membs <- which(memb == k)
    block_graph <- Y[k_membs, k_membs]
    if(ls_method=='vb') {
      ls <- latent_space_vb(block_graph, D = D)
    } else if(tolower(ls_method)=="mds") {
      inv.logit <- function(x) exp(x) / (1 + exp(x))
      logit <- function(x) log(x / (1-x))
      N <- nrow(block_graph)
      ls <= list()
      require(igraph)
      dists <- shortest.paths(graph.adjacency(block_graph, mode='undirected', diag=FALSE))
      dists[dists==Inf] <- max(dists[dists < Inf]) + 1 # for nodes in different components, set distance to n
      # Z is initialized using mutidimensional scaling of the initial distance matrix
      ls$Z <- cmdscale(D_start, D)
      ls$Z <- ls$Z+ runif(nrow(ls$Z) * ncol(ls$Z), min=-.1, max=.1) # in case two nodes have the same position
      ls$Z[,1] <- ls$Z[,1] - mean(ls$Z[,1])
      ls$Z[,2] <- ls$Z[,2] - mean(ls$Z[,2])
      ls$beta <- logit(mean(block_graph[lower.tri(block_grpah, diag=FALSE)])) + mean(dists[lower.tri(dists,diag=FALSE)])
      ls$sigma <- sd(c(ls$Z))# / sqrt(nrow(block_graph))
      ls$theta <- c(ls$beta, log(ls$sigma))
    } else {
      ls <- latent_space_mle(block_graph, D = D)
    }
    Z[k_membs,] <- ls$Z
    beta[k] <- ls$beta
    sigma[k] <- ls$sigma
    theta[k,] <- ls$theta
    if(!is.null(plot_folder)) {
      p <- plot_ls(block_graph, ls$Z) + theme(panel.grid=element_blank()) + ggtitle(paste0("Block ", k, "; theta = (", round(ls$theta[1],3), ", ", round(ls$theta[2],3), ")"))
      ggsave(paste0(plot_folder, "/latent_space_", k, ".png"), p, height=6, width=6, units='in')     
    }
  }
  out$Z <- Z
  out$beta <- beta
  out$sigma <- sigma
  out$theta <- theta

  out$mu <- mu <- colMeans(theta)
  out$Sigma <- Sigma <- cov(theta)
  
  if(!is.null(plot_folder)) {
    require(ellipse)
    thetas_df <- data.frame(theta1=theta[,1], theta2=theta[,2])
    l <- data.frame(ellipse(cov2cor(Sigma), scale=sqrt(diag(Sigma)), centre=mu))
    p <- ggplot(thetas_df) + geom_point(aes(x=theta1, y=theta2)) + geom_point(x=mu[1], y=mu[2], col='red') + geom_path(data=l, aes(x,y)) + theme_bw() + theme(panel.grid=element_blank()) + xlab(expression(theta[1])) + ylab(expression(theta[2])) 
    ggsave(paste0(plot_folder, "/thetas.png"),p, height=4, width=4, units='in')
  }
  
  pred <- out$B[out$gamma, out$gamma]
  for(k in 1:length(out$beta)) {
    membs <- which(out$gamma == k)
    if(length(membs) > 1) pred[membs, membs] <- logit.inv(out$beta[k] - as.matrix(dist(out$Z[membs,])))
  }
  diag(pred) <- 0
  out$pred <- pred
  return(out)
}


