###############################################################################
# Multi-resolution blockmodel
#
# file: mcmc_sampler
#
# This file contains functions for sampling from the posterior distribution of 
# our model
# Author: tedwestling
###############################################################################

#######################################
# function: block_latent_MCMC
#  
# -------------------------------------
# Driver function for the sampling algorithm
#######################################

block_latent_MCMC <- function(Y, D, K, burn_in, n_samples, thin, v, epsilon=.1, rZ=1, Atheta=matrix(c(2,1,1,1), nrow=2), alpha0=.1, a0=NULL, b0=NULL, m0=c(0,0), s0=0.01, psi0=matrix(c(5.1,0,0,5.1), nrow=2), nu0=5.1, verbose=TRUE, memb_start=NULL, plot_init=FALSE, likelihood=FALSE, true_gamma=NULL, postprocess=TRUE, sample_membs=TRUE, record_acc_probs=FALSE, debug_output=FALSE, perturb_init=TRUE) {
  require(gtools)
  require(MASS)
  require(MCMCpack)
  N <- nrow(Y)
  miss_inds <- which(is.na(Y), arr.ind = TRUE)
  miss_inds <- miss_inds[miss_inds[,1] < miss_inds[,2],]
  
  #impute missing values using observed node-wise edge frequencies for initialization
  # obs_prob <- colMeans(Y, na.rm=TRUE)
  # exp_prob <- outer(obs_prob, obs_prob)
  # Yimp <- Y
  # if(sum(is.na(Y)) > 0) Yimp[is.na(Yimp)] <- exp_prob[is.na(Yimp)]
  

  # if no starting values are provided, initialize using the three-stage procedure
  # par_t will contain the latest parameter values, which will be used for proposal distributions
  if(verbose) cat("Initalizing...")
  par_t <- three_stage_est(Y, D, K=K, plot_folder = NULL)
  if(perturb_init) {
    perturb_nodes <- sample(1:N, ceiling(N/10), replace=FALSE)
    for(node in perturb_nodes) {
      par_t$gamma[node] <- sample((1:par_t$K)[-par_t$gamma[node]], 1)
    }
  }
  par_t$Y <- Y
  par_t$Y[is.na(par_t$Y)] <- par_t$pred[is.na(par_t$Y)]
  par_t$block_n <- c(table(par_t$gamma))
  par_t$s <- matrix(NA, nrow=par_t$K, ncol=par_t$K)
  for(i in 1:par_t$K) for(j in 1:i) par_t$s[i,j] <- par_t$s[j,i] <- sum(par_t$Y[par_t$gamma == i, par_t$gamma == j])
  par_t$pi <- par_t$block_n/N
  K <- par_t$K
  dens <- sum(Y, na.rm=TRUE)/(N^2 - N)
  if(is.null(a0)) a0 <- 10 * dens
  if(is.null(b0)) b0 <- 10 * (1-dens)
  
  par_t$N <- N
  par_t$dists <- lapply(1:K, function(k) as.matrix(dist(par_t$Z[which(par_t$gamma == k),])))
  if(plot_init) plot_blocked_matrix(Y, par_t$gamma)
  if(verbose) {
    cat("Done initializing. K =", K, "\n")
    for(k in 1:K) cat("theta", k, "=", par_t$theta[k,], "\n")
  }
  cD <- 2 * gamma((D + 1)/2) / gamma(D/2)
  expit <- function(x) 1/(1 + exp(-x))
  
  minval <- par_t$mu[1] - cD * exp(par_t$mu[2])
  if(digamma(a0) - digamma(b0) > minval) par_t$mu[1] <- 1 + digamma(a0) - digamma(b0) + cD * exp(par_t$mu[2])
    
  # Create list to store our samples in
  if(nrow(miss_inds) > 0) Ymiss <- matrix(NA, nrow=n_samples, ncol=nrow(miss_inds))
  beta <- matrix(NA, nrow=n_samples, ncol=K)
  sigma <- matrix(NA, nrow=n_samples, ncol=K)
  gamma <- matrix(NA, nrow=n_samples, ncol=N)
  pi <- matrix(NA, nrow=n_samples, ncol=K)
  mu <- matrix(NA, nrow=n_samples, ncol=2)
  B <- array(dim=c(n_samples,K,K))
  Z <- array(dim=c(n_samples,N,D))
  Sigma <- array(dim=c(n_samples,2,2))
  if(likelihood) loglik <- matrix(NA, nrow=n_samples, ncol=9)
  
  if(record_acc_probs) {
    joint_acc_probs <- matrix(NA, nrow=n_samples, ncol=N)
    pos_acc_probs <- matrix(NA, nrow=n_samples, ncol=N)
    theta_acc_probs <- matrix(NA, nrow=n_samples, ncol=K)
  }
  if(verbose) cat("Starting sampling... \n")
  start_time <- Sys.time()
  total_samples <- v*(burn_in + n_samples * thin)
  sample_num <- 0
  for(t in 1:total_samples) {
    keeping_sample <- (t > v * burn_in & t %% (v * thin) == 0)
    if(keeping_sample) sample_num <- sample_num + 1 #(t - v * burn_in) / (v * thin)
    if(verbose & (t %% 100 == 0)) {  
      time_left <- (c(difftime(Sys.time(), start_time, units="secs")) / t) * (total_samples - t)
      cat("Iteration", t, "of", total_samples, ",", round(c(time_left)/60, 2), "minutes left\n")
      cat("Current mu value",round(par_t$mu,2),"constraints enforced\n")
    }
    if(plot_init & (t %% 100 == 0))  plot_blocked_matrix(par_t$Y, par_t$gamma)
    # Sample missing values from current parameter estimates first because we need Y with no missing values for the other updates
    if(nrow(miss_inds) > 0) {
      if(debug_output) cat("Missing edge updates\n")
      par_t$Y[miss_inds] <- par_t$Y[miss_inds[,c(2,1)]] <- Ymiss[t,] <- missing_edge_update(par_t, miss_inds)
      for(i in 1:par_t$K) for(j in 1:i) par_t$s[i,j] <- par_t$s[j,i] <- sum(par_t$Y[par_t$gamma == i, par_t$gamma == j])
    }
    
    if(t %% v == 0) {
      #if(verbose) cat("  Performing joint proposal.\n")
      for(i in 1:N) {
        joint_prop <- joint_proposal(i, epsilon, rZ, par_t)
        if(keeping_sample & record_acc_probs) joint_acc_probs[sample_num, i] <- joint_prop$prob
        #if(verbose) cat("    Proposal acceptance probability", joint_prop$prob, "\n")
        if(runif(1)  < joint_prop$prob) {
          prev_block <- par_t$gamma[i]
          new_block <- joint_prop$gamma_i
          par_t$gamma[i] <- new_block
          par_t$Z[i,] <- joint_prop$Zi
          par_t$dists[[new_block]] <- as.matrix(dist(par_t$Z[par_t$gamma == new_block,]))
          if(prev_block != new_block) {
            par_t$dists[[prev_block]] <- as.matrix(dist(par_t$Z[par_t$gamma == prev_block,]))
            par_t$block_n[prev_block] <- par_t$block_n[prev_block] - 1
            par_t$block_n[new_block] <- par_t$block_n[new_block] + 1
            for(j in 1:K) {
              par_t$s[prev_block,j] <- par_t$s[j,prev_block] <- sum(par_t$Y[par_t$gamma == prev_block, par_t$gamma == j])
              par_t$s[new_block, j] <- par_t$s[j,new_block] <- sum(par_t$Y[par_t$gamma == new_block, par_t$gamma == j])
            }
          }
        }
      }
    }
    
    # MH position steps
    #if(verbose) cat("  Proposing positions.\n")
    if(debug_output) cat("Position proposals\n")
    for(i in 1:N) {
      pos_prop <- position_proposal(i, rZ, par_t)
      if(keeping_sample & record_acc_probs) pos_acc_probs[sample_num,i] <- pos_prop$prob
      if(runif(1) < pos_prop$prob) {
        par_t$Z[i,] <- pos_prop$Zi
        par_t$dists[[par_t$gamma[i]]] <- as.matrix(dist(par_t$Z[par_t$gamma == par_t$gamma[i],]))
      } 
    }
    #if(verbose) cat("    Av proposal acceptance probability", mean(pmin(pos_acc_probs[t,], 1)), "\n")
    
    # Other updates
    # randomize order
    #if(verbose) cat("  Updating others.\n")
    ord <- sample(5, 5, replace=FALSE)
    for(up in ord) {
      if(up == 1) {
        if(debug_output) cat("pi update\n") 
        par_t$pi <- pi_update(par_t, alpha0)
      }
      if(up == 2) {
        if(debug_output) cat("B update\n")
        par_t$B <- B_update(par_t, a0, b0)
      }
      if(up == 3) {
        if(debug_output) cat("beta, sigma updates\n")
        for(k in 1:K) {
          theta_prop <- theta_proposal(k, Atheta, par_t)
          if(keeping_sample & record_acc_probs) theta_acc_probs[sample_num,k] <- theta_prop$prob
          if(runif(1) < theta_prop$prob) {
            par_t$beta[k] <- theta_prop$thetak[1]
            par_t$sigma[k] <- exp(theta_prop$thetak[2])
            par_t$theta[k,] <- theta_prop$thetak
          }
        }
        #if(verbose) cat("    Av beta acceptance probability", mean(pmin(beta_acc_probs[t,], 1)), "\n")
      }
      if(up == 4) {
        if(debug_output) cat("mu update\n")
        par_t$mu <- mu_update(par_t, m0, s0, a0, b0)
      }
      if(up == 5) {
        if(debug_output) cat("Sigma update \n")
        par_t$Sigma <- Sigma_update(par_t, m0, s0, psi0, nu0)
      }
    }
    
    # Save samples
    if(keeping_sample) {
      beta[sample_num,] <- par_t$beta
      sigma[sample_num,] <- par_t$sigma
      gamma[sample_num,] <- par_t$gamma
      pi[sample_num,] <- par_t$pi
      mu[sample_num,] <- par_t$mu
      Sigma[sample_num,,] <- par_t$Sigma
      B[sample_num,,] <- par_t$B
      Z[sample_num,,] <- par_t$Z
      if(likelihood) loglik[sample_num,] <- log_likelihood(Y, par_t, alpha0, a0, b0, m0, s0, psi0, nu0)
    }
  }
  output <- list(beta=beta, sigma=sigma, pi=pi, mu=mu, Sigma=Sigma, B=B, Z=Z, gamma=gamma)
  if(likelihood) output$loglik <- loglik
  if(record_acc_probs) {
    output$joint_acc_probs <- joint_acc_probs
    output$pos_acc_probs <- pos_acc_probs
    output$theta_acc_probs <- theta_acc_probs
  }
  if(nrow(miss_inds) > 0) {
    output$Ymiss <- Ymiss
    output$miss_inds <- miss_inds
  }
  if(postprocess) return(postprocess_MCMC(output, Y, true_gamma=true_gamma))
  else return(output)
}


log_likelihood <- function(Y, par_t, alpha0, a0, b0, m0, s0, psi0, nu0) {
  K <- length(par_t$pi)
  d <- ncol(par_t$Z)
  cD <- 2 * gamma((d+1)/2) / gamma(d/2)
  pi_ll <- sum(par_t$block_n * log(par_t$pi))
  Yoff_ll <- sum(par_t$s[lower.tri(par_t$s, diag=FALSE)] * log(par_t$B[lower.tri(par_t$B, diag=FALSE)]))
  theta <- cbind(par_t$beta, log(par_t$sigma))
  theta_ll <- - K * log(2*3.1415)- (K/2) * log(det(par_t$Sigma)) + sum(apply(theta, 1, function(thetak) - .5 * t(thetak - par_t$mu) %*% solve(par_t$Sigma) %*% (thetak - par_t$mu)))
  Ydiag_ll <- z_ll <- 0
  for(k in 1:K) {
    k_membs <- which(par_t$gamma == k) 
    etak <- par_t$beta[k] - par_t$dists[[k]]
    mat <- Y[k_membs, k_membs] * etak - log1p(exp(etak))
    Ydiag_ll <- Ydiag_ll + sum(mat[lower.tri(mat, diag=FALSE)]) 
    z_ll <- z_ll+ sum(dnorm(c(par_t$Z[k_membs,]), 0, par_t$sigma[k], log=TRUE))
  }
  piprior_ll <- sum(log(ddirichlet(par_t$pi, rep(alpha0, K))))
  Bprior_ll <- sum(log(dbeta(par_t$B[lower.tri(par_t$B, diag=FALSE)], a0, b0)))
  muprior_ll <- - log(2*3.1415) - (1/2) * log(det(par_t$Sigma/ s0)) - (1/2) * t(par_t$mu - m0) %*% solve(par_t$Sigma / s0) %*% (par_t$mu - m0)
  Sigmaprior_ll <- log(diwish(par_t$Sigma, nu0, psi0))
  c(Ydiag_ll=Ydiag_ll, Yoff_ll=Yoff_ll, pi_ll=pi_ll, z_ll, theta_ll=theta_ll, piprior_ll=piprior_ll, Bprior_ll=Bprior_ll, muprior_ll=muprior_ll, Sigmaprior_ll=Sigmaprior_ll)
}

#######################################
# functions: joint_proposal, etc
#  
# -------------------------------------
# These functions calculate proposals 
# and acceptance probabilities for
# each step in the sampler. To avoid
# passing back and forth large objects, 
# the functions assume that 
# par_t is in the parent environment
# and only take proposal-specific 
# variances, etc as arguments.
#######################################

joint_proposal <- function(i, epsilon, rZ, par_t) {
  Y <- par_t$Y
  Yi <- Y[i,]
  N <- par_t$N
  K <- par_t$K
  D <- ncol(par_t$Z)
  
  gammai_t <- par_t$gamma[i]
  Zi_t <- par_t$Z[i,]
  #propose new memberships and positions
  edge_blocks <- Yi * par_t$gamma
  block_counts <- sapply(1:K, function(k) sum(edge_blocks == k))
  if(par_t$block_n[gammai_t] <= 2) {
    lambdai <- rep(0, K)
    lambdai[gammai_t] <- 1
    gammai_star <- gammai_t
  } else {
    unscaled <- (block_counts+epsilon) / (par_t$block_n +1)
    lambdai <- unscaled/sum(unscaled)
    gammai_star <- sample(K, 1, prob=lambdai)
  }
  
  block_n_star <- par_t$block_n
  block_n_star[gammai_t] <- block_n_star[gammai_t] - 1
  block_n_star[gammai_star] <- block_n_star[gammai_star] + 1
  
  unscaled_rev <- (block_counts+epsilon) / (block_n_star +1)
  lambdai_rev <- unscaled_rev / sum(unscaled_rev)
  
  if(gammai_star == gammai_t) {
    mi <- Zi_t
    Zistar <- rnorm(D, mi, rZ)
    mi_rev <- Zistar
  } else {
    friends_in_block <- which(Yi == 1 & par_t$gamma == gammai_star)
    if(length(friends_in_block) > 1) mi <- colMeans(par_t$Z[friends_in_block,])
    if(length(friends_in_block) == 1) mi <- par_t$Z[friends_in_block,]
    if(length(friends_in_block) == 0) mi <- c(0,0)
    
    Zistar <- rnorm(D, mi, rZ)
    
    # calculate reverse
    friends_in_block <- which(Yi  == 1 & par_t$gamma == gammai_t)
    if(length(friends_in_block) > 1) mi_rev <- colMeans(par_t$Z[friends_in_block,])
    if(length(friends_in_block) == 1) mi_rev <- par_t$Z[friends_in_block,]
    if(length(friends_in_block) == 0) mi_rev <- c(0,0)
  } 
  
  # Prior for Z
  log_acc_prob <- sum(dnorm(Zistar, 0, par_t$sigma[gammai_star], log = TRUE)) - sum(dnorm(Zi_t, 0, par_t$sigma[gammai_t], log = TRUE))
  
  # Prior for gamma
  log_acc_prob <- log_acc_prob + log(par_t$pi[gammai_star]) - log(par_t$pi[gammai_t])
  
  # likelihood
  for(k in 1:K) {
    blockk_t <- which(par_t$gamma == k)
    blockk_t <- blockk_t[blockk_t != i]
    if(k == gammai_star) {
      dij_star <- dists_to_vec(Zistar, par_t$Z[blockk_t,])
      eta_star <- par_t$beta[k] - dij_star
      log_acc_prob <- log_acc_prob + sum(Yi[blockk_t]*eta_star -log1p(exp(eta_star)))
    } else {
      nk <- length(blockk_t)
      sk <- sum(Yi[blockk_t])
      log_acc_prob <- log_acc_prob + sk * log(par_t$B[k,gammai_star]) + (nk - sk)*log(1-par_t$B[k,gammai_star])
    }
    
    if(k == gammai_t) {
      dij_t <- dists_to_vec(Zi_t, par_t$Z[blockk_t,])
      eta_t <- par_t$beta[k] - dij_t
      log_acc_prob <- log_acc_prob - sum(Yi[blockk_t]*eta_t -log1p(exp(eta_t)))
    } else {
      nk <- length(blockk_t)
      sk <- sum(Yi[blockk_t])
      log_acc_prob <- log_acc_prob - sk * log(par_t$B[k,gammai_t]) - (nk - sk)*log(1-par_t$B[k,gammai_t])
    }
  }
  
  # Transition probabilities
  log_acc_prob <- log_acc_prob + log(lambdai_rev[gammai_t]) - log(lambdai[gammai_star])
  
  log_acc_prob <- log_acc_prob + sum(dnorm(Zi_t, mi_rev, rZ, log=TRUE)) - sum(dnorm(Zistar, mi, rZ, log=TRUE))
  
  return(list(gamma_i=gammai_star,
              Zi=Zistar,
              prob=as.numeric(exp(log_acc_prob))))
}

position_proposal <- function(i, rZ, par_t) {
  Y <- par_t$Y
  N <- par_t$N
  K <- par_t$K
  D <- ncol(par_t$Z)
  
  # propose new position and calculate probability of accepting new position, return as list 
  Zi_t <- par_t$Z[i,]
  Zi_star <- rnorm(D, Zi_t, rZ)
  
  # acceptance probability related to prior on Zi
  k <- par_t$gamma[i]
  log_acc_prob <- sum(Zi_t^2 - Zi_star^2)/(2* par_t$sigma[k]^2)
  
  # acceptance probability related to data
  kmembs <- which(par_t$gamma == k)
  kmembs <- kmembs[kmembs != i]
  Zj <- par_t$Z[kmembs,]
  Yij <- Y[i,kmembs]
  
  dij_t <- dists_to_vec(Zi_t, Zj)
  dij_star <- dists_to_vec(Zi_star, Zj)
  
  log_acc_prob <- log_acc_prob + sum((dij_t - dij_star)[Yij==1]) + sum(log(1 + exp(par_t$beta[k] - dij_t)) - log(1 + exp(par_t$beta[k]- dij_star)))
  
  return(list(Zi=Zi_star,
              prob=exp(log_acc_prob)))
}


pi_update <- function(par_t, alpha0) {
  c(rdirichlet(1, alpha0 + par_t$block_n))
}

B_update <- function(par_t, a0, b0) {
  expit <- function(x) 1/(1 + exp(-x))
  N <- par_t$N
  K <- par_t$K
  nt <- c(table(par_t$gamma))
  newB <- matrix(NA, nrow=K, ncol=K)
  D <- ncol(par_t$Z)
  cD <- 2 * gamma((D + 1)/2) / gamma(D/2)
  for(k in 2:K) {
    for(l in 1:(k-1)) {
      #while(TRUE) {
        prop <- rbeta(1, a0 + par_t$s[k, l], b0 + par_t$block_n[k] * par_t$block_n[l] - par_t$s[k, l])
        #if(prop < min(exp_value(par_t$theta[k,1], exp(par_t$theta[k,2])), exp_value(par_t$theta[l,1],  exp(par_t$theta[l,2])))) break
      #}
      newB[k, l] <- newB[l, k] <- prop
    }
  }
  return(newB)
}

theta_proposal <- function(k, Atheta, par_t) {
  expit <- function(x) 1/(1 + exp(-x))
  Y <- par_t$Y
  N <- par_t$N
  K <- par_t$K
  D <- ncol(par_t$Z)
  cD <- 2 * gamma((D + 1)/2) / gamma(D/2)
  
  # make proposal
  betak_t <- par_t$beta[k]
  sigmak_t <- par_t$sigma[k]
  thetak_t <- c(betak_t, log(sigmak_t))
  #while(TRUE) {
    thetak_star <- mvrnorm(1, thetak_t, Atheta)
    betak_star <- thetak_star[1]
    sigmak_star <- exp(thetak_star[2])
    
    # Check that proposal satisfies constraints
  #  if(sum(par_t$B[k, -k] >= exp_value(betak_star, sigmak_star)) == 0) break
  #}

  # part of acceptance prob regarding theta
  log_acc_prob <- -(1/2) * t(thetak_star - par_t$mu) %*% solve(par_t$Sigma) %*% (thetak_star - par_t$mu) + (1/2) * t(thetak_t - par_t$mu) %*% solve(par_t$Sigma) %*% (thetak_t - par_t$mu)
  
  # part of  acceptance prob involving block Y
  etak_t <- betak_t -  par_t$dists[[k]]
  etak_star <- betak_star -  par_t$dists[[k]]
  mat <- - log1p(exp(etak_star)) + log1p(exp(etak_t))
  
  log_acc_prob <- log_acc_prob + (par_t$s[k,k]/2) * (betak_star - betak_t)  + sum(mat[lower.tri(mat, diag=FALSE)])
  
  # Part involving z
  this_z <- c(par_t$Z[which(par_t$gamma==k),])
  log_acc_prob <- log_acc_prob +  (D * par_t$block_n[k]) * (-log(sigmak_star) + log(sigmak_t)) + (-sigmak_star^{-2} + sigmak_t^{-2}) * sum(this_z^2) / (2)
  
  list(thetak=thetak_star,
       prob=c(exp(log_acc_prob)))
}

mu_update <- function(par_t,  m0, s0, a0, b0) {
	require(msm)
#new strategy is to update the two mu's marginally so we enforce the constraint one one dimension at a time
	#sample size 
  N <- par_t$N
  	#number of blocks
  K <- par_t$K
  	#Joint update would be a sample from a MVN with mean m and sd matrix co
  	#m is vector of length 2 co is 2 by 2 matrix
  m <- (colSums(par_t$theta) + s0 * m0) / (K + s0)
  co <- par_t$Sigma / (K + s0)
  #print(co)
#dimension of latent space
  D <- ncol(par_t$Z)
#constant in front of constraint
  cD <- 2 * gamma((D + 1)/2) / gamma(D/2)
  mu_star<-c(NA,NA)
  #####print everything before the update
 # print("printing inputs")
 # print("m")
 # print(m)
 # print("co")
 # print(co)
 # print("cD")
 # print(cD)
  #############
  #note that the OLD value of mu2 goes into the constraint, not the current average of thetas that goes into the gibbs update 
  constraintmu1=digamma(a0) - digamma(b0)+(cD * exp(par_t$mu[2]))
  #update mu1
  rho=co[1,2]/sqrt(co[1,1]*co[2,2])
  #print(rho)
  mu1mean=m[1]+((sqrt(co[1,1])/sqrt(co[2,2]))*rho*(par_t$mu[2]-m[2]))
  mu1sd=sqrt((1-rho^2)*co[1,1])
  ######
  #print("mu1 constraint")
  #print(constraintmu1)
  #print("mean of mu1")
  #print(mu1mean)
  #print("mu 1sd")
  #print(mu1sd)
  ######
  mu_star[1]<-rtnorm(1,lower=constraintmu1,mean=mu1mean,sd=mu1sd)
  #note this is an UPPER bound now; sign flips bc of subtraction
  #note now we use the UPDATED value of mu[1]
  constraintmu2=log((digamma(a0)-digamma(b0)-mu_star[1])/(-cD))
  mu2mean=m[2]+((sqrt(co[2,2])/sqrt(co[1,1]))*rho*(mu_star[1]-m[1]))
  mu2sd=sqrt((1-rho^2)*co[2,2])
    ######
  #print("mu2 constraint")
  #print(constraintmu2)
  #print("mean of mu2")
  #print(mu2mean)
  #print("mu 2sd")
  #print(mu2sd)
  ######
  mu_star[2]<-rtnorm(1,upper=constraintmu2,mean=mu2mean,sd=mu2sd)
 #print("mu_star")
 #	print(mu_star)
  return(mu_star)
}

Sigma_update <- function(par_t, m0, s0, psi0, nu0) {
  N <- par_t$N
  K <- par_t$K
  S <- psi0 + cov(par_t$theta) * (K-1) + (K * s0/ (K + s0)) * tcrossprod(colMeans(par_t$theta) - m0)
  v <- K + nu0
  riwish(v, S)
}

missing_edge_update <- function(par_t, miss_inds) {
  apply(miss_inds, 1, function(ind) {
    i <- ind[1]
    j <- ind[2]
    gammai <- par_t$gamma[i]
    gammaj <- par_t$gamma[j]
    p <- ifelse(gammai == gammaj,
                1/(1 + exp(-par_t$beta[gammai] + sqrt(sum((par_t$Z[i,] - par_t$Z[j,])^2)))),
                par_t$B[gammai, gammaj])
    rbinom(1,1,p)
  })
} 

dists_to_vec <- function(vec, mat) {
  if(is.matrix(mat)) apply(mat, 1, function(row) sqrt(sum((row - vec)^2)))
  else sqrt(sum((mat - vec)^2))
}

#######################################
# function: postprocess_MCMC
#  
# -------------------------------------
# Postprocesses the output of a call to
# block_latent_mcmc using the algorithm
# described in the appendix.
# Rotates membership samples toward
# fixed_memb.
######################################

postprocess_MCMC <- function(output, Y, fixed_memb) {
  require(clue)
  K <- ncol(output$pi)
  N <- ncol(output$gamma)
  D <- dim(output$Z)[3]
  print("Postprocessing posteriors samples...")
  nsamp <- nrow(output$gamma)
  
  for(s in 1:nrow(output$gamma)) {
    tab <- table(output$gamma[s,], fixed_memb)
    optimal_perm <- c(solve_LSAP(as.matrix(tab), maximum=TRUE))
    ord <- order(optimal_perm)
    new_gamma <- rep(NA, N)
    for(val in 1:K) new_gamma[output$gamma[s,] == val] <- optimal_perm[val]
    output$gamma[s,] <- new_gamma
    output$beta[s,] <- output$beta[s, ord]
    output$sigma[s,] <- output$sigma[s, ord]
    output$pi[s,] <- output$pi[s, ord]
    output$B[s,,] <- output$B[s, ord, ord]
  }
  
  ## Find marginal posterior memberships probs
  gamma_freqs <- apply(output$gamma, 2, function(col) sapply(1:K, function(k) sum(col == k)/length(col)))
  for(k in 1:K) {
    # Get nodes with any probability of being in this block
    nodes_in <- which(gamma_freqs[k,] > 0)
    # Construct a fixed set of positions
    # First get the distances between the nodes across samples 
    dists <- array(dim=c(nsamp, length(nodes_in), length(nodes_in)))
    for(s in 1:nsamp) {
      # only want distances for nodes actually in the block now
      nodes_in_now <- which(output$gamma[s,] == k)
      subs <- which(nodes_in %in% nodes_in_now)
      if(length(subs) > 1)  dists[s, subs, subs] <- as.matrix(dist(output$Z[s,nodes_in[subs],]))# / (sigma[s,k] * sqrt(length(nodes_in_now)))))
    }
    # Take mean distances
    mean_dists <- apply(dists, c(2,3), mean, na.rm=TRUE)
    # What we put for missing means (i.e. never in block together) doesn't matter since it will have weight 0
    mean_dists[is.na(mean_dists)] <- max(mean_dists, na.rm=TRUE)
    # Find the number of times they were in the block together as a weight matrix
    weights <- apply(dists, c(2,3), function(vec) sum(!is.na(vec)))
    
    # Apply MDS to mean
    require(smacof)
    fit <- try(smacofSym(mean_dists, ndim = min(D, nrow(mean_dists) - 1), weightmat=weights), silent=TRUE)
    fudge <- 1
    while(class(fit)[1] == "try-error" & fudge < 100) {
      weights[weights == fudge - 1] <- fudge
      fit <- try(smacofSym(mean_dists, ndim = min(D, nrow(mean_dists) - 1), weightmat=weights), silent=TRUE)
      fudge <- fudge + 1
    }
    Z0 <- fit$conf * mean(fit$delta / fit$dhat, na.rm=TRUE)
    if(is.null(dim(Z0))) Z0 <- matrix(Z0)
    if(ncol(Z0) < D) Z0 <- cbind(Z0, matrix(0, nrow=nrow(mean_dists), ncol=D-ncol(Z0)))
    #Z0 <- cmdscale(mean_dists, D)
    #Z0 <- Z0 / mean(rowSums(Z0^2))
    
    for(s in 1:nrow(output$gamma)) {
      #print(s)
      nodes_in_now <- which(output$gamma[s,] == k)
      Z0_now <- Z0[which(nodes_in %in% nodes_in_now),]
      Z_now <- output$Z[s,nodes_in_now[nodes_in_now %in% nodes_in],]
      if(length(which(nodes_in %in% nodes_in_now)) == 0) next
      if(length(which(nodes_in %in% nodes_in_now)) == 1) new_Z_now <- Z0_now
      else new_Z_now <- procrustes(Z_now, Z0_now, scale=FALSE)
      output$Z[s,nodes_in_now[nodes_in_now %in% nodes_in],] <- new_Z_now
      # Set the positions of the nodes in this block whose posterior mode is not this block to NA
      output$Z[s,nodes_in_now[!(nodes_in_now %in% nodes_in)],] <- NA
    }
  }
  
  return(output)
}


procrustes <- function(Z, Z0, scale=TRUE){
  D <- ncol(Z)
  #Procrustes transform; gives rotation,reflection,trranslation, and possibly scaling
  #of Z closest to Z0
  Z0_means <- colMeans(Z0)
  Z0 <- apply(Z0, 2, function(col) col - mean(col))
  Z <- apply(Z, 2, function(col) col - mean(col))
  
  if(scale) {
    s0 <- sqrt(sum(Z0^2) / nrow(Z0))
    s <- sqrt(sum(Z^2) / nrow(Z))
    Z <- Z * s0/s
  }
  
  #   A <- crossprod(t(Z0) %*% Z) 
  #   eA <- eigen(A,symmetric=TRUE)
  #   Ahalf <- eA$vec[,1:D]%*%diag(sqrt(eA$val[1:D]))%*%t(eA$vec[,1:D])
  #   
  #   Z_trans <- t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z))
  
  Asvd <- svd(crossprod(Z, Z0))
  Z_trans <- Z %*% Asvd$u %*% t(Asvd$v)
  for(i in 1:D) Z_trans[,i] <- Z_trans[,i] + Z0_means[i]
  Z_trans
}


exp_value <- function(beta, sigma) {
  g <- function(x) expit(x) * dchisq((x - beta)^2 / (2*sigma^2), 2) * (x - beta) / sigma^2 #expit(beta - sqrt(2) * sigma * sqrt(x)) * dchisq(x, 2)
  int <- try(-integrate(g, -Inf, beta)$value, silent=TRUE)
  if(class(int) == "try-error") {
    print(c(beta, sigma))
  } else return(int)
}


