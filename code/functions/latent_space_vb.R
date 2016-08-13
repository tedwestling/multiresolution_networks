###############################################################################
# Multi-resolution blockmodel
#
# file: latent_space_vb

# This file contains functions for fitting a latent space model to a network
# using variational bayes. We do not expect this to provide a consistent estimate of 
# the model parameters or a realistic sense of uncertainty; we only use it
# to initialize our MCMC sampler.
#
# Author: tedwestling
###############################################################################

#######################################
# function: latent_space_vb
#  
#-------------------------------------
# Driver function to perform VB
#######################################

latent_space_vb <- function(Y, D=2, m0=0, t0=.5, a0=2.1, b0=1.1, m_start=NULL, l_start=NULL, b_start=NULL, s_start=NULL, t_start=NULL, node_effects=rep(1, nrow(Y))) {
  if(length(Y) == 1) {
    return(list(Z=matrix(rep(0,D), nrow=1),
                beta=m0,
                sigma=sqrt(b0/(a0-1)),
                theta=c(m0, log(sqrt(b0/(a0-1)))),
                dists=matrix(0, nrow=1),
                elbo=NA,
                var_param=NULL))
  }
  
  inv.logit <- logit.inv <- expit <- function(x) exp(x) / (1 + exp(x))
  logit <- function(x) log(x / (1-x))
  N <- nrow(Y)
  a <- a0 + N/2
  
  if(is.null(l_start)) {
    ## Initialize Z 
    if(D >= N-1) l_start <- matrix(rnorm(N*D, 0, 1), ncol=D)
    else {
      # Distances are initialized to the minimum path-length
      D_start <- shortest.paths(graph.adjacency(Y, mode='undirected', diag=FALSE))
      D_start[D_start==Inf] <- max(D_start[D_start < Inf]) + 1 # for nodes in different components, set distance to n
      # Z is initialized using mutidimensional scaling of the initial distance matrix
      l_start <- cmdscale(D_start, D)
      l_start <- l_start+ runif(nrow(l_start) * ncol(l_start), min=-.01, max=.01) # in case two nodes have the same position
      if(ncol(l_start) < D) l_start <- cbind(l_start, matrix(0, nrow=nrow(l_start), ncol=D-ncol(l_start)))
    }
    l_start[,1] <- l_start[,1] - mean(l_start[,1])
    l_start[,2] <- l_start[,2] - mean(l_start[,2])
    D_start <- as.matrix(dist(l_start))
  }
  
  if(is.null(s_start)) s_start <- rep(1, N)
  
  if(is.null(m_start)) m_start <- logit(mean(Y[lower.tri(Y, diag=FALSE)])) + mean(D_start[lower.tri(D_start,diag=FALSE)])
  if(is.infinite(m_start)) m_start <- m0
  
  if(is.null(b_start)) b_start <- (b0 + .5*sum(D / s_start + rowSums(l_start^2)))#/N
  
  if(is.null(t_start)) t_start <- .2
  
  ## now use optim to find maximizer of ELBO
  par_start <- c(m_start, log(t_start), log(b_start), c(l_start), log(s_start))
  
  # BFGS optimization
  while(TRUE) { 
    optimized <- optim(par_start, latent_space_elbo, gr=latent_space_elbo_grad, Y=Y, D=D, m0=m0, t0=t0, a0=a0, b0=b0, method="BFGS", control=list(fnscale=-1,  maxit=500))
    if(optimized$convergence == 0) break
    par_start <- par_start + rnorm(length(par_start), 0, .3)
  }
  
  m <- optimized$par[1]
  t <- exp(optimized$par[2])
  b <- exp(optimized$par[3])
  sigma <- sqrt(b/(a-1))
  l <- matrix(optimized$par[4:(3+N*D)], ncol=D)
  s <- exp(optimized$par[(4+N*D):(3+N+N*D)])
  dists <- as.matrix(dist(l))
  cD <- 2 * gamma((D+1)/2) / gamma(D/2)
  theta <- c(m, log(sigma))
  
  return(list(Z=l,
              beta=m,
              sigma=sigma,
              theta=theta,
              dists=dists,
              elbo=optimized$value,
              var_param=list(t=t, a=a, b=b, s=s)))
}

#######################################
# functions: latent_space_elbo
#            latent_space_elbo_grad
#-------------------------------------
# Returns the elbo and
# gradient of the elbo
# respectively of the latent space model
# for the given variational paramters,
# priors, and 
# observed undirected binary adjacency
# Y and number of dimensions D
#######################################

latent_space_elbo <- function(par_vec, Y, D, m0, t0, a0, b0, node_effects=rep(nrow(Y), 1)){
  inv.logit <- function(x) exp(x) / (1 + exp(x))
  N <- nrow(Y)
  m <- par_vec[1]; t <- exp(par_vec[2]); b <- exp(par_vec[3]); l <- matrix(par_vec[4:(3+(N*D))],ncol=D); s <- exp(par_vec[(4+(N*D)):(3+N+N*D)])
  a <- a0 + N/2
  
  dists <- as.matrix(dist(l))^2
  norms <- rowSums(l^2)
  
  Epij <- m - sqrt(dists + D * outer(s^{-1}, s^{-1}, "+"))
  data_part <- Y * Epij - log(1 + exp(Epij + 1/(2*t)))
  pos_kl <- .5 * (N * (digamma(a) - log(b) + D) - sum(a * (D / s + norms) / (b) + log(s)))
  prior_kl <- -a0 *log(b) + lgamma(a) + (a0 - a) * digamma(a) - a *(b0 / b - 1) - (1/2)*log(t) - (t0/2) * ((m - m0)^2 + 1/t) + a0 *log(b0)- lgamma(a0) + log(t0)/2
  as.numeric(sum(data_part[lower.tri(data_part, diag=FALSE)]) + pos_kl + prior_kl)
}

latent_space_elbo_grad <- function(par_vec, Y, D, m0, t0, a0, b0) {
  inv.logit <- function(x) exp(x) / (1 + exp(x))
  N <- nrow(Y)
  m <- par_vec[1]; t <- exp(par_vec[2]); b <- exp(par_vec[3]); l <- matrix(par_vec[4:(3+(N*D))],ncol=D); s <- exp(par_vec[(4+(N*D)):(3+N+N*D)])
  a <- a0 + N/2
  dists <- as.matrix(dist(l))^2
  norms <- rowSums(l^2)
  ezizj <- dists + D * outer(s^{-1}, s^{-1}, "+")
  Epij <- m - sqrt(ezizj)
  sum_logit <- sum((inv.logit(Epij + 1/(2*t)))[lower.tri(Y, diag=FALSE)])
  
  grad_m <- sum(Y[lower.tri(Y, diag=FALSE)]) - sum_logit - t0*(m - m0)
  grad_t <- (1/(2*t^2)) * (t0 + sum_logit) - (1/(2*t))
  grad_b <- a * (-b^{-1} + b^{-2}*(b0 + .5 * sum(D / s + norms)))
  resid <- Y - logit.inv(Epij + (1/(2*t)))
  diag(resid) <- 0
  grad_l <- sapply(1:D, function(k) {
    lk_mat <- rep(1,N) %*% t(l[,k])
    diff_mat <- t(lk_mat) - lk_mat # produces l_ik - l_jk
    full_mat <- -diff_mat * resid * ezizj^{-1/2}
    diag(full_mat) <- 0
    #full_mat <- full_mat * upper.tri(full_mat, diag=FALSE)
    rowSums(full_mat) - (a/b) * l[,k]
  })
  grad_s <- rowSums(.5*D * s^{-2} * ezizj^(-1/2) * resid) + a*D / (2*b*s^2) - 1/(2*s)
  as.numeric(c(grad_m, t * grad_t, b * grad_b, c(grad_l), s * grad_s ))
}


#######################################
# function: plot_ls
#  
#-------------------------------------
# Plots a latent space estimate using
# ggplot2
#######################################

plot_ls <- function(Y, Z) {
  pt_df <- data.frame(x=Z[,1], y=Z[,2])
  edges <- which(Y==1, arr.ind=TRUE)
  edge_df <- data.frame(x=pt_df$x[edges[,1]], xend=pt_df$x[edges[,2]], y=pt_df$y[edges[,1]], yend=pt_df$y[edges[,2]])
  ggplot(pt_df) + geom_point(aes(x,y, position="jitter")) + geom_segment(data=edge_df, aes(x=x, xend=xend, y=y, yend=yend))+ theme_bw()
}
