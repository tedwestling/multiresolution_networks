###############################################################################
# Multi-resolution blockmodel
#
# file: indian_village_model_estimation.R
# 
# This file performs all model estimation for the visitcome relation of
# village 59 from the Karnataka village dataset as desribed in the paper.
#
# Author: tedwestling
###############################################################################
##set your working directory to the top level multiresolution_networks folder
#setwd("")

# Source header file (should be in top level of working directory)
source('header.R')

# We set the seed numerous times throughout the document at points at natural reload points
set.seed(22311) 
vilno <- 59
edgetype <- 'visitcome'

# Load in household-level network data for the particular edgetype and village
network <- as.matrix(read.csv(paste0('data/indian_village_raw/1. Network Data/Adjacency Matrices/adj_', edgetype, '_HH_vilno_', vilno, '.csv'), header=FALSE))

# Remove nodes with degree 0
fullnetwork <- network
zros <- which(colSums(network) == 0)
network <- network[-zros, -zros]
N <- nrow(network)

# Choose the number of blocks
Kchoice <- choose_K(network, Krange = 2:ceiling(N/4), n_splits = 20, n_folds = 10, assortative=TRUE, degree_correct=FALSE, threshold=TRUE)
save(Kchoice, file='data/results/village_59_Kchoice.Rdata')

#if you've already run the script to choose K, just load it to make pics
#load('data/results/village_59_Kchoice.Rdata')

# Make some plots 
ggplot(Kchoice$mean_aucs) +
  geom_errorbar(aes(K, ymin=lower95CI_auc, ymax=upper95CI_auc), color='grey') +
  geom_point(aes(K, mean_auc)) +
  geom_line(aes(K, mean_auc)) +
  geom_point(x=Kchoice$minoptK_auc, y=Kchoice$mean_aucs$mean_auc[Kchoice$mean_aucs$K == Kchoice$minoptK_auc], color='red') +
  theme_bw() +
  xlab("K") + ylab("Avg. CV-AUC") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_auc, ")"))

ggplot(Kchoice$mean_aucs) +
  geom_errorbar(aes(K, ymin=lower95CI_mse, ymax=upper95CI_mse), color='grey') +
  geom_point(aes(K, mean_mse)) +
  geom_line(aes(K, mean_mse)) +
  geom_point(x=Kchoice$minoptK_mse, y=Kchoice$mean_aucs$mean_mse[Kchoice$mean_aucs$K == Kchoice$minoptK_mse], color='red') +
  theme_bw() +
  xlab("K") + ylab("Avg. CV-MSE") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_mse, ")"))

ggplot(Kchoice$mean_aucs) +
  geom_errorbar(aes(K, ymin=lower95CI_mpi, ymax=upper95CI_mpi), color='grey') +
  geom_point(aes(K, mean_mpi)) +
  geom_line(aes(K, mean_mpi)) +
  geom_point(x=Kchoice$minoptK_mpi, y=Kchoice$mean_aucs$mean_mpi[Kchoice$mean_aucs$K == Kchoice$minoptK_mpi], color='red') +
  theme_bw() +
  xlab("K") + ylab("Avg. CV-MPI") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_mpi, ")"))

# Set the number of blocks to 6

K <- Khat <- 6

##set the a0 and b0 parameters according to assortativity restriction
dens <- sum(network, na.rm=TRUE)/(N^2 - N)
b0=1
a0=(b0*(10*dens))/(1-(10*dens))

# Perform the MCMC
# In our implementation, we ran four chains in parallel.  
chain1 <- block_latent_MCMC(network, D=2, K=Khat, burn_in=0, n_samples=8000, thin=20, rZ=3, v=3, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = Khat,a0=a0, b0=b0, postprocess = FALSE, plot_init=TRUE)

save(chain1, file='data/results/village_59_mcmc_chain1.Rdata')

