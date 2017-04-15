###############################################################################
# Multi-resolution blockmodel
#
# file: sbm_mcmc

# This file contains a function for performing Gibbs sampling for the stochastic
# block model (SBM). It allows missing edges.
# 
# Author: tedwestling
###############################################################################

#################################################################
# Function: smb_mcmc
# 
# Performs Gibbs sampling for the stochastic blockmodel.
# Y contains the adjacency matrix
# K is the number of blocks (if NULL chooses K with best
#   AIC from spectral clustering)
# a0_within, b0_within specify the prior beta dist. within blocks
# a0_between, b0_between are beta param. between blocks
# other parameters control spectral clustering initialization
#################################################################

sbm_mcmc <- function(Y, K=NULL, alpha0=NULL, a0_within=4, b0_within=1, a0_between=1, b0_between=10, burnin=1000, samples=1000, thin=10, tau=mean(colSums(Y)), normalize=TRUE, threshold=TRUE, gamma=ifelse(threshold, 1,0), true_gamma=NULL, print_samples=FALSE, allow_zero=FALSE) {
  require(gtools)
  # Initialize with spectral clustering
  N <- nrow(Y)
  any_miss <- any(is.na(Y[lower.tri(Y, diag=FALSE)]))
  if(any_miss) {
    miss_inds <- which(is.na(Y), arr.ind = TRUE)
    miss_inds <- miss_inds[miss_inds[,1] < miss_inds[,2],]
    #impute missing values using observed node-wise edge frequencies for initialization
    obs_prob <- colMeans(Y, na.rm=TRUE)
    exp_prob <- outer(obs_prob, obs_prob)
    if(nrow(miss_inds) > 0) Y[is.na(Y)] <- exp_prob[is.na(Y)]
    diag(Y) <- NA
  }
  
  print("Initializing...")
  if(is.null(K)) {
    all_clusts <- spectral_cluster(Y, Krange=(2:floor(nrow(Y)/10)), tau=tau, normalize=normalize, threshold=threshold, gamma=gamma)
    aics <- unlist(lapply(all_clusts, function(c) c$sbm_aic))
    clust <- all_clusts[[which.min(aics)]]  
  } else {
    clust <- spectral_cluster(Y, Krange=K, tau=tau, normalize=normalize, threshold=threshold, gamma=gamma)
    clust <- clust[[1]]
  }
  par_t <- list(memb=clust$clusters, Y=Y, K=length(unique(clust$clusters)), B=clust$Bhat)
  if(print_samples) plot_blocked_matrix(Y, par_t$memb)
  if(is.null(alpha0)) alpha0 <- K
  
  # Create matrices to store samples
  tot_samp <- burnin + samples * thin
  memb_samples <- matrix(nrow=samples, ncol=N)
  pi_samples <- matrix(nrow=samples, ncol=par_t$K)
  B_samples <- array(dim=c(samples, par_t$K, par_t$K))
  if(any_miss) miss_samples <- matrix(nrow=samples, ncol=nrow(miss_inds))
  logliks <- rep(NA, tot_samp)
  
  print("Sampling...")
  k <- 1
  for(t in 1:tot_samp) {
    if(t %% 100 == 0) print(paste0("Iteration ", t, " of ", tot_samp))
    if(print_samples & (t %% 100 == 0)) plot_blocked_matrix(par_t$Y, par_t$memb, sort=FALSE)
    #plot_blocked_matrix(par_t$Y, par_t$memb)
    par_t$pi <- sample_pi(par_t$memb, alpha0, K)
    par_t$B <- sample_B(par_t, a0_within, b0_within, a0_between, b0_between)
    # Ordering the blocks by within-block density enforces identifiability, see Nowicki & Snijders
    ord <- order(diag(par_t$B))
    # par_t$B <- par_t$B[ord, ord]
    # par_t$pi <- par_t$pi[ord]
    
    if(any_miss) par_t$Y[miss_inds] <- par_t$Y[miss_inds[,2:1]] <- sample_miss(par_t, miss_inds)
    for(i in 1:N) par_t$memb[i] <- sample_memb(par_t, i, allow_zero)
    if(t > burnin & (t - burnin) %% thin == 0) {
      pi_samples[k,] <- par_t$pi
      B_samples[k,,] <- par_t$B
      memb_samples[k,] <- par_t$memb
      if(any_miss) miss_samples[k,] <- par_t$Y[miss_inds]
      k <- k + 1
    }
    logliks[t] <- loglik(par_t, alpha0, a0_within, b0_within, a0_between, b0_between)
  }
  if(any_miss) return(list(memb_samples=memb_samples, B_samples=B_samples, pi_samples=pi_samples, miss_samples=miss_samples, logliks=logliks))
  else return(list(memb_samples=memb_samples, B_samples=B_samples, pi_samples=pi_samples, logliks=logliks))
}

sample_pi <- function(membs_t, alpha0, K) {
  tab <- tabulate(membs_t, nbins=K) #sapply(1:K, function(k) sum(membs_t == k))
  c(rdirichlet(1, alpha0 + tab))
}

sample_memb <- function(par_t, i, allow_zero) {
  if(!allow_zero & table(par_t$memb)[par_t$memb[i]] == 1) return(par_t$memb[i])
  Yi <- par_t$Y[i,-i]

  p <- par_t$pi * apply(par_t$B[par_t$memb[-i],]^Yi * (1 - par_t$B[par_t$memb[-i],])^(1-Yi), 2, prod)
  p <- p / sum(p)
  newblock <- sample(1:par_t$K, 1, prob=p)
  return(newblock)
}


sample_B <- function(par_t, a0_within, b0_within, a0_between, b0_between) {
  B <- matrix(nrow=par_t$K, ncol=par_t$K)
  for(k in 1:par_t$K) {
    for(l in 1:k) {
      membk <- which(par_t$memb == k)
      membl <- which(par_t$memb == l)
      Skl <- sum(par_t$Y[membk, membl], na.rm=TRUE)
      Tkl <- sum(1-par_t$Y[membk, membl], na.rm=TRUE)
      if(k == l) B[k,l] <- B[l,k] <- rbeta(1, a0_within + Skl, b0_within + Tkl)
      if(k != l) B[k,l] <- B[l,k] <- rbeta(1, a0_between + Skl, b0_between + Tkl)
    }
  }
  return(B)
}

sample_miss <- function(par_t, miss_inds) {
  prob <- par_t$B[cbind(par_t$memb[miss_inds[,1]], par_t$memb[miss_inds[,2]])]
  rbinom(nrow(miss_inds), 1, prob)
}

loglik <- function(par_t, alpha0, a0_within, b0_within, a0_between, b0_between) {
  prob <- par_t$B[par_t$memb, par_t$memb]
  sum((par_t$Y * prob + (1 - par_t$Y) * (1 - prob))[lower.tri(par_t$Y, diag=FALSE)]) + 
        sum(dbeta(diag(par_t$B), a0_within, b0_within, log=TRUE)) + 
        sum(dbeta(par_t$B[lower.tri(par_t$B, diag=FALSE)], a0_between, b0_between, log=TRUE)) + 
        sum(log(par_t$pi[par_t$Z])) +
        log(ddirichlet(par_t$pi, rep(alpha0, par_t$K)))
}