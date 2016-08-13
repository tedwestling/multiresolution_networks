###############################################################################
# Multi-resolution blockmodel
#
# file: spectral_clustering

# This file contains functions for performing spectral clustering.
# Author: tedwestling
###############################################################################

#############################################################
# Function: spectral_cluster
# 
# Performs regularized spectral clustering on the set S 
# (algorithm reference: http://papers.nips.cc/paper/5099-regularized-spectral-clustering-under-the-degree-corrected-stochastic-blockmodel.pdf)
# on the provided adjacency matrix Y for the provided range of
# number of blocks Krange. 
# tau is the regularization parameter, defaults to mean degree. 
#  tau = 0 corresponds to non-regularized.
# normalize (boolean) normalizes the eigenvectors, defaults to TRUE
# threshold (boolean) only clusters nodes with large enough leverage scores
#  (defaults to FALSE)
# gamma controls thresholding, defaults to 0 if threshold=FALSE,
#  1 if threshold=TRUE (only nodes with threshold scores >= gamma/sqrt(N))
#  are used to cluster
#############################################################

spectral_cluster <- function(Y, Krange=(2:floor(nrow(Y)/10)), tau=mean(colSums(Y)), normalize=TRUE, threshold=TRUE, gamma=ifelse(threshold, 1,0), plot=FALSE, nstart=5, niter=5, verbose=FALSE, degree_correct=FALSE, assortative=TRUE) {
  require(rARPACK)
  
  miss_inds <- which(is.na(Y), arr.ind = TRUE)
  miss_inds <- miss_inds[miss_inds[,1] < miss_inds[,2],]
  
  #impute missing values using observed node-wise edge frequencies for initialization
  #if(degree_correct) {
    obs_prob <- colMeans(Y, na.rm=TRUE)
    exp_prob <- outer(obs_prob, obs_prob)
  #} else exp_prob <- mean(Y[lower.tri(Y, diag=FALSE)], na.rm=TRUE)
  if(sum(is.na(Y)) > 0) Y[is.na(Y)] <- exp_prob[is.na(Y)]
  
  # Construct matrices
  N <- nrow(Y)
  dc <- degree_correct
  if(!assortative) {
    dtau_sqrtinv <- 1/sqrt(colSums(Y) + tau)
    #Dtau_sqrtinv <- diag()
    L <- outer(dtau_sqrtinv, dtau_sqrtinv) * Y #Dtau_sqrtinv %*% Y %*% Dtau_sqrtinv
    eig <- eigs(L, k = max(Krange), which='LR')
    ordered_eigvals <- order(eig$values, decreasing=TRUE)
    vectors <- eig$vectors[,ordered_eigvals]
  } else {
    degrees <- colSums(Y)
    D <- diag(degrees)
    r <- sqrt(mean(degrees))
    H <- (r^2 - 1) * matrix(1, nrow=N, ncol=N) - r * Y + D
    eig <- eigen(H, symmetric=TRUE)
    negs <- which(eig$values < 0)
    vectors <- eig$vectors[,negs]
  }

  lapply(Krange, function(K) {
    #print(K)
    if(verbose) print(paste0("K = ", K))
    
    if(!assortative) eigvecs <- vectors[,1:K]
    else eigvecs <- vectors
    eigvec_lengths <- sqrt(rowSums(eigvecs^2))
    if(any(eigvec_lengths == 0)) eigvec_lengths[eigvec_lengths == 0] <- 1
    # normalize vectors
    if(normalize) eigvecs <- eigvecs / eigvec_lengths
    # run k-means on thresholded eigenvectors
    made_threshold <- eigvec_lengths >= gamma/sqrt(N)
    kmeans_fit <- kmeans(eigvecs[made_threshold,], centers=K, nstart=nstart)
    # If the eigenvectors were thresholded, need to assign ones below threshold
    if(length(kmeans_fit$cluster) < N) {
      clusters <- rep(NA, N)
      clusters[made_threshold] <- kmeans_fit$cluster
      clusters[!made_threshold] <- apply(subset(eigvecs, !made_threshold),1, function(xi_norm) {
        # assign thresholded to closest center found by kmeans
        as.numeric(which.min(colSums((t(kmeans_fit$centers) - xi_norm)^2)))
      })
    } else clusters <- kmeans_fit$cluster
    # compute estimated B matrix
    cluster_sizes <- as.numeric(table(clusters))
    poss_rel <- cluster_sizes %*% t(cluster_sizes) - diag(cluster_sizes)
    obs_rel <- matrix(NA, nrow=K, ncol=K)
    for(i in 1:K) for(j in 1:i) obs_rel[i,j] <- obs_rel[j,i] <- sum(Y[clusters == i, clusters == j])
    Bhat <-  obs_rel / poss_rel
    gap_stat <- min(diag(Bhat)) - max(Bhat[lower.tri(Bhat)])
    Bhat[Bhat == 0] <- .000001
    Bhat[Bhat == 1] <- .999999
    Bhat[is.na(Bhat)] <- sum(Y) / (N^2 - N) # For within-block probs for blocks of size one, just set to overall density as a prior since no data (0/0)
    
    if(dc) {
      degrees <- colSums(Y)
      block_deg_sums <- sapply(1:K, function(k) sum(degrees[clusters == k]))
      norm_degrees <- degrees / block_deg_sums[clusters]
      newpred <- outer(norm_degrees, norm_degrees) * obs_rel[clusters, clusters]
    } else {
      newpred <- Bhat[clusters, clusters]
    }
    newpred[newpred > 1]  <- .99999
    if(niter > 0) {
      Y[miss_inds] <- newpred[miss_inds]
      for(iter in 1:niter) {

        if(!assortative) {
          dtau_sqrtinv <- 1/sqrt(colSums(Y) + tau)
          #Dtau_sqrtinv <- diag()
          L <- outer(dtau_sqrtinv, dtau_sqrtinv) * Y #Dtau_sqrtinv %*% Y %*% Dtau_sqrtinv
          eig <- eigs(L, k = max(Krange), which='LR')
          ordered_eigvals <- order(eig$values, decreasing=TRUE)
          vectors <- Re(eig$vectors[,ordered_eigvals])
        } else {
          degrees <- colSums(Y)
          D <- diag(degrees)
          r <- sqrt(mean(degrees))
          H <- (r^2 - 1) * matrix(1, nrow=N, ncol=N) - r * Y + D
          eig <- eigen(H, symmetric=TRUE)
          ordered_eigvals <- order(eig$values, decreasing=FALSE)
          #negs <- which(eig$values < 0)
          vectors <- Re(eig$vectors[,ordered_eigvals])
        }

        if(!assortative) eigvecs <- vectors[,1:K]
        else eigvecs <- vectors[,1:K]
        eigvec_lengths <- sqrt(rowSums(eigvecs^2))
        if(any(eigvec_lengths == 0)) eigvec_lengths[eigvec_lengths == 0] <- 1
        # normalize vectors
        if(normalize) eigvecs <- eigvecs / eigvec_lengths
        # run k-means on thresholded eigenvectors
        made_threshold <- eigvec_lengths >= gamma/sqrt(N)
        kmeans_fit <- kmeans(eigvecs[made_threshold,], centers=K, nstart=nstart)
        # If the eigenvectors were thresholded, need to assign ones below threshold
        if(length(kmeans_fit$cluster) < N) {
          clusters <- rep(NA, N)
          clusters[made_threshold] <- kmeans_fit$cluster
          clusters[!made_threshold] <- apply(subset(eigvecs, !made_threshold),1, function(xi_norm) {
            # assign thresholded to closest center found by kmeans
            as.numeric(which.min(colSums((t(kmeans_fit$centers) - xi_norm)^2)))
          })
        } else clusters <- kmeans_fit$cluster
        # compute estimated B matrix

        cluster_sizes <- as.numeric(table(clusters))
        poss_rel <- cluster_sizes %*% t(cluster_sizes) - diag(cluster_sizes)
        obs_rel <- matrix(NA, nrow=K, ncol=K)
        for(i in 1:K) for(j in 1:i) obs_rel[i,j] <- obs_rel[j,i] <- sum(Y[clusters == i, clusters == j])
        Bhat <-  obs_rel / poss_rel
        gap_stat <- min(diag(Bhat)) - max(Bhat[lower.tri(Bhat)])
        Bhat[Bhat == 0] <- .000001
        Bhat[Bhat == 1] <- .999999
        Bhat[is.na(Bhat)] <- sum(Y) / (N^2 - N) # For within-block probs for blocks of size one, just set to overall density as a prior since no data (0/0)
        
        if(dc) {
          degrees <- colSums(Y)
          block_deg_sums <- sapply(1:K, function(k) sum(degrees[clusters == k]))
          norm_degrees <- degrees / block_deg_sums[clusters]
          newpred <- outer(norm_degrees, norm_degrees) * obs_rel[clusters, clusters]
        } else {
          newpred <- Bhat[clusters, clusters]
        }
        newpred[newpred > 1]  <- .99999
        Y[miss_inds] <- newpred[miss_inds]
      }
    }
    newpred[newpred == 0] <- .00001
    newpred[newpred == 1] <- .99999
    if(!dc) norm_degrees <- NULL
    return(list(K=K,
                clusters=clusters,
                Bhat=Bhat,
                pred=newpred,
                norm_degrees=norm_degrees))
  })
  
}

#############################################################
# Function: choose_K
# 
# Performs "n_splits" repititions of "n_folds"-fold CV 
# using spectral clustering. Returns the CV AUC, MSE, and MPI
# in a data frame. Additional arguments passed on to spectral
# clustering
#############################################################


choose_K <- function(Y, Krange=(2:floor(nrow(Y)/5)), n_splits = 5, n_folds = 10, verbose=TRUE, ...) {
  library(AUC)
  N <- nrow(Y)
  allinds <- as.matrix(expand.grid(row=1:N, col=1:N))
  allinds <- allinds[allinds[,1] < allinds[,2],]
  aucs <- expand.grid(K=Krange, split=1:n_splits)
  aucs$auc <- aucs$mpi <- aucs$mse <- NA
  for(split in 1:n_splits) {
    folds <- rep(1:n_folds, length.out = nrow(allinds))
    folds <- folds[sample(nrow(allinds), nrow(allinds), replace=FALSE)]
    pred <- list(rep(NA, length(Krange))) 
    for(k in 1:length(Krange)) pred[[k]] <- matrix(NA, nrow=N, ncol=N)
    for(fold in 1:n_folds) {
      if(verbose) print(paste0("Split ", split, ", fold ", fold))
      test_inds <- allinds[folds == fold,]
      Ymiss <- Y
      Ymiss[test_inds] <- NA
      Ymiss <- (Ymiss + t(Ymiss)) /2
      spec <- spectral_cluster(Ymiss, K=Krange, verbose=FALSE, ...)
      for(k in 1:length(Krange)) pred[[k]][test_inds] <- spec[[k]]$pred[test_inds]
    }
    for(k in 1:length(Krange)) {
      roc_est <- roc(predictions=pred[[k]][allinds], labels=as.factor(Y[allinds]))
      #plot(roc_est, add=TRUE)
      aucs$auc[aucs$K == Krange[k] & aucs$split == split] <- auc(roc_est)   
      aucs$mse[aucs$K == Krange[k] & aucs$split == split] <- mean((pred[[k]][allinds] - Y[allinds])^2)
      aucs$mpi[aucs$K == Krange[k] & aucs$split == split] <-  mean(Y[allinds] * log(pred[[k]][allinds]) + (1-Y[allinds]) * log(1-pred[[k]][allinds])) + 1
    }
  }

  mean_aucs <- ddply(aucs, .(K), function(df) {
    tt_ts <- t.test(df$auc)
    tt_mse <- t.test(df$mse)
    tt_mpi <- t.test(df$mpi)
    data.frame(mean_auc=tt_ts$estimate, lower95CI_auc=tt_ts$conf.int[1], upper95CI_auc=tt_ts$conf.int[2],
               mean_mse=tt_mse$estimate, lower95CI_mse=tt_mse$conf.int[1], upper95CI_mse=tt_mse$conf.int[2],
               mean_mpi=tt_mpi$estimate, lower95CI_mpi=tt_mpi$conf.int[1], upper95CI_mpi=tt_mpi$conf.int[2])
  }) 
  optK_auc <- mean_aucs$K[which.max(mean_aucs$mean_auc)]
  minoptK_auc <- min(mean_aucs$K[which(mean_aucs$mean_auc >= mean_aucs$lower95CI_auc[mean_aucs$K == optK_auc])])
  optK_mse <- mean_aucs$K[which.min(mean_aucs$mean_mse)]
  minoptK_mse <- min(mean_aucs$K[which(mean_aucs$mean_mse <= mean_aucs$upper95CI_mse[mean_aucs$K == optK_mse])])
  optK_mpi <- mean_aucs$K[which.max(mean_aucs$mean_mpi)]
  minoptK_mpi <- min(mean_aucs$K[which(mean_aucs$mean_mpi >= mean_aucs$lower95CI_mpi[mean_aucs$K == optK_mpi])])
  list(aucs=aucs, mean_aucs=mean_aucs, optK_auc=optK_auc, minoptK_auc=minoptK_auc, optK_mse=optK_mse, minoptK_mse=minoptK_mse, optK_mpi=optK_mpi, minoptK_mpi=minoptK_mpi)
}

#############################################################
# Function: plot_blocked_matrix
# 
# Plots an adjacency matrix (network) sorted by a membership
# vector. Also plots lines to indicate block breaks. 
# argument "sort" indicates whether blocks should be sorted 
# by size (TRUE by default)
#############################################################

plot_blocked_matrix <- function(network, memb, sort=TRUE) {
  N <- nrow(network)
  if(sort) {
    clust_ord <- order(order(table(memb)))
    ord <- order(clust_ord[memb])
  } else {
    clust_ord <- 1:length(memb)
    ord <- order(clust_ord[memb])
  }
  parnow <- par(no.readonly=TRUE)
  par(mai=c(0,0,0,0), omi=rep(0,4), xaxs='i', yaxs='i', mgp=c(0,0,0))
  plot(0:N, 0:N, type='n', axes=FALSE, xlab="", ylab="", asp=1)
  rasterImage(as.raster(1-network[ord,ord]),ybottom=0, xleft=0, ytop=N, xright=N, interpolate=FALSE)
  k <- max(memb)
  for(j in 2:k) {
    breakpt <- min(which(clust_ord[memb][ord] == j))
    segments(x0=c(breakpt-1, 0), y0=c(0,N-breakpt+1), x1=c(breakpt-1,N), y1=c(N,N-breakpt+1))
  }
  par(parnow)
}

