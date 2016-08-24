###############################################################################
# Multi-resolution blockmodel
#
# file: indian_village_postprocess.R
# 
# This file does convergence checks, does post processing, and makes plots for 
#the India village data.  Assuming that we have already run India_village_estimation.R
## and have results.  
#
#
#Setup here to use four chains, but can be run with more/fewer.  To use this file, you'll
#need to be sure that the names of the 
#
#
# Author: tedwestling
###############################################################################
rm(list = ls())
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
K <- Khat <- 6

thinfirst=T
thinfac=4 #keep every fourth
#note on these runs we've already thinned by 10 within the sampler
####for the longer chains, add some thinning

load('data/results/village_59_mcmc_chain1.Rdata')
chain1a=chain1
load('data/results/village_59_mcmc_chain2.Rdata')
chain2=chain1
load('data/results/village_59_mcmc_chain3.Rdata')
chain3=chain1
load('data/results/village_59_mcmc_chain4.Rdata')
chain4=chain1
chain1=chain1a

if(thinfirst==T){
	#list(beta=beta, sigma=sigma, pi=pi, mu=mu, Sigma=Sigma, B=B, Z=Z, gamma=gamma)
	keep=which(c(1:length(chain1$beta[,1]))%%thinfac==0)
	chain1=list(beta=chain1$beta[keep,], sigma=chain1$sigma[keep,], pi=chain1$pi[keep,], mu=chain1$mu[keep,], Sigma=chain1$Sigma[keep,,], B=chain1$B[keep,,], Z=chain1$Z[keep,,], gamma=chain1$gamma[keep,])
	#
	chain2=list(beta=chain2$beta[keep,], sigma=chain2$sigma[keep,], pi=chain2$pi[keep,], mu=chain2$mu[keep,], Sigma=chain2$Sigma[keep,,], B=chain2$B[keep,,], Z=chain2$Z[keep,,], gamma=chain2$gamma[keep,])
	#
	chain3=list(beta=chain3$beta[keep,], sigma=chain3$sigma[keep,], pi=chain3$pi[keep,], mu=chain3$mu[keep,], Sigma=chain3$Sigma[keep,,], B=chain3$B[keep,,], Z=chain3$Z[keep,,], gamma=chain3$gamma[keep,])
	#
	chain4=list(beta=chain4$beta[keep,], sigma=chain4$sigma[keep,], pi=chain4$pi[keep,], mu=chain4$mu[keep,], Sigma=chain4$Sigma[keep,,], B=chain4$B[keep,,], Z=chain4$Z[keep,,], gamma=chain4$gamma[keep,])

}

# Postprocess Samples
all_chains <- list(beta=rbind(chain1$beta, chain2$beta, chain3$beta, chain4$beta),
                sigma=rbind(chain1$sigma, chain2$sigma, chain3$sigma, chain4$sigma),
                 pi=rbind(chain1$pi, chain2$pi, chain3$pi, chain4$pi),
                 mu=rbind(chain1$mu, chain2$mu, chain3$mu, chain4$mu),
                 Sigma=abind(chain1$Sigma, chain2$Sigma, chain3$Sigma, chain4$Sigma, along=1),
                 B=abind(chain1$B, chain2$B, chain3$B, chain4$B, along=1),
                 Z=abind(chain1$Z, chain2$Z, chain3$Z, chain4$Z, along=1),
                 gamma=rbind(chain1$gamma, chain2$gamma, chain3$gamma, chain4$gamma))

set.seed(303854)
# Use asortative spectral clustering estimate as fixed membership vector to rotate towards
spec <- spectral_cluster(network, Krange = Khat, assortative = TRUE, plot=FALSE, degree_correct = FALSE)
#uncomment to plot spectral clusters
#plot_blocked_matrix(network, spec[[1]]$clusters)

 mcmc_samplesc1 <- postprocess_MCMC(chain1, network, fixed_memb = spec[[1]]$clusters)
mcmc_samplesc2 <- postprocess_MCMC(chain2, network, fixed_memb = spec[[1]]$clusters)
mcmc_samplesc3 <- postprocess_MCMC(chain3, network, fixed_memb = spec[[1]]$clusters)
mcmc_samplesc4 <- postprocess_MCMC(chain4, network, fixed_memb = spec[[1]]$clusters)

#save your post-processed samples
#save.image("data/results/village_59_mcmc_strongass_v1_postprocessed.Rdata")
save(network,chain1,chain2,chain3,chain4,mcmc_samplesc1,mcmc_samplesc2,mcmc_samplesc3,mcmc_samplesc4,N,K,Khat,vilno,file="data/results/village_59_mcmc_strongass_v1_postprocessed.Rdata")

