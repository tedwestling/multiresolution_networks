#should htere be a header file? needs to load the other scripts
#logit.inv is defined in header!!!!
###need to load igraph

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

# Source header file (should be in top level of working directory)
#source('header.R')

# We set the seed numerous times throughout the document at points at natural reload points
#set.seed(22311) 

jobnumber=as.numeric(commandArgs(TRUE)[1])
print(jobnumber)
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
#Kchoice <- choose_K(network, Krange = 2:ceiling(N/4), n_splits = 20, n_folds = 10, assortative=TRUE, degree_correct=FALSE, threshold=TRUE)
#save(Kchoice, file='data/results/village_59_Kchoice.Rdata')

load('data/results/village_59_Kchoice.Rdata')
# library(ggplot2)
# # Make some plots 
# ggp1<-ggplot(Kchoice$mean_aucs) +
  # geom_errorbar(aes(K, ymin=lower95CI_auc, ymax=upper95CI_auc), color='grey') +
  # geom_point(aes(K, mean_auc)) +
  # geom_line(aes(K, mean_auc)) +
  # geom_point(x=Kchoice$minoptK_auc, y=Kchoice$mean_aucs$mean_auc[Kchoice$mean_aucs$K == Kchoice$minoptK_auc], color='red') +
  # theme_bw() +
  # xlab("K") + ylab("Avg. CV-AUC") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_auc, ")"))
# ggsave('plots/kchoice1.png', ggp1, width=4.5, height=2.5, units='in', scale=1.5)


# ggp2<-ggplot(Kchoice$mean_aucs) +
  # geom_errorbar(aes(K, ymin=lower95CI_mse, ymax=upper95CI_mse), color='grey') +
  # geom_point(aes(K, mean_mse)) +
  # geom_line(aes(K, mean_mse)) +
  # geom_point(x=Kchoice$minoptK_mse, y=Kchoice$mean_aucs$mean_mse[Kchoice$mean_aucs$K == Kchoice$minoptK_mse], color='red') +
  # theme_bw() +
  # xlab("K") + ylab("Avg. CV-MSE") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_mse, ")"))
# ggsave('plots/kchoice2.png', ggp2, width=4.5, height=2.5, units='in', scale=1.5)


# ggp3<-ggplot(Kchoice$mean_aucs) +
  # geom_errorbar(aes(K, ymin=lower95CI_mpi, ymax=upper95CI_mpi), color='grey') +
  # geom_point(aes(K, mean_mpi)) +
  # geom_line(aes(K, mean_mpi)) +
  # geom_point(x=Kchoice$minoptK_mpi, y=Kchoice$mean_aucs$mean_mpi[Kchoice$mean_aucs$K == Kchoice$minoptK_mpi], color='red') +
  # theme_bw() +
  # xlab("K") + ylab("Avg. CV-MPI") + ggtitle(paste0("Village ", vilno, " (Khat = ", Kchoice$minoptK_mpi, ")"))
# ggsave('plots/kchoice3.png', ggp3, width=4.5, height=2.5, units='in', scale=1.5)

# Set the number of blocks to 6

K <- Khat <- 6

# Perform the four chains of MCMC
# WARNING: each of these chains takes approximately two days to run
#set.seed(3291)
#####################################################################
source("code/functions/mcmc_sampler.R")
source("code/functions/three_stage_estimation.R")
source("code/functions/spectral_clustering.R")
source("code/functions/latent_space_vb.R")
library(igraph)
logit <- function(x) log(x/(1-x))
expit <- inv.logit <- logit.inv <- function(x) 1/ (1 + exp(-x))
#####################################################################

#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0
dens <- sum(network, na.rm=TRUE)/(N^2 - N)
b0=1
a0=(b0*(10*dens))/(1-(10*dens))
#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0
#########define a0 b0

chain1 <- block_latent_MCMC(network, D=2, K=Khat, burn_in=0, n_samples=8000, thin=20, rZ=3, v=3, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = Khat,a0=a0, b0=b0, postprocess = FALSE, plot_init=TRUE)

#random sequence to be sure we're actually getting different chains
#diffseq=round(abs(rnorm(1,5000,1000)))
save(chain1, file=paste('village_59_mcmc_strongass_v1_',jobnumber,'.Rdata',sep=""))
#rm("chain1")

#set.seed(667553)
#chain2 <- block_latent_MCMC(network, D=2, K=Khat, burn_in=20000, n_samples=2500, thin=100, rZ=3, v=5, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = Khat, postprocess = FALSE, plot_init=TRUE)

# save(chain2, file='data/results/village_59_mcmc_chain2.Rdata')
# rm("chain2")

# set.seed(839144)
# chain3 <- block_latent_MCMC(network, D=2, K=Khat, burn_in=20000, n_samples=2500, thin=100, rZ=3, v=5, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = Khat, postprocess = FALSE, plot_init=TRUE)

# save(chain3, file='data/results/village_59_mcmc_chain3.Rdata')
# rm("chain3")

# set.seed(73737)
# chain4 <- block_latent_MCMC(network, D=2, K=Khat, burn_in=20000, n_samples=2500, thin=100, rZ=3, v=5, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = Khat, postprocess = FALSE, plot_init=TRUE)

# save(chain4, file='data/results/village_59_mcmc_chain4.Rdata')
# rm("chain4")

# load('data/results/village_59_mcmc_chain1.Rdata')
# load('data/results/village_59_mcmc_chain2.Rdata')
# load('data/results/village_59_mcmc_chain3.Rdata')
# load('data/results/village_59_mcmc_chain4.Rdata')

# # Postprocess Samples
# all_chains <- list(beta=rbind(chain1$beta, chain2$beta, chain3$beta, chain4$beta),
                 # sigma=rbind(chain1$sigma, chain2$sigma, chain3$sigma, chain4$sigma),
                 # pi=rbind(chain1$pi, chain2$pi, chain3$pi, chain4$pi),
                 # mu=rbind(chain1$mu, chain2$mu, chain3$mu, chain4$mu),
                 # Sigma=abind(chain1$Sigma, chain2$Sigma, chain3$Sigma, chain4$Sigma, along=1),
                 # B=abind(chain1$B, chain2$B, chain3$B, chain4$B, along=1),
                 # Z=abind(chain1$Z, chain2$Z, chain3$Z, chain4$Z, along=1),
                 # gamma=rbind(chain1$gamma, chain2$gamma, chain3$gamma, chain4$gamma))

# set.seed(303854)
# spec <- spectral_cluster(network, Krange = Khat, assortative = TRUE, plot=FALSE, degree_correct = FALSE)
# plot_blocked_matrix(network, spec[[1]]$clusters)
# mcmc_samples <- postprocess_MCMC(all_chains, network, true_gamma = spec[[1]]$clusters)

# save(mcmc_samples, file='data/results/village_59_mcmc_samples.Rdata')
# load('data/results/village_59_mcmc_samples.Rdata')

# set.seed(47676728)
# gamma_mode <- as.numeric(apply(mcmc_samples$gamma, 2, function(memb) names(table(memb))[which.max(table(memb))]))

# # Plot blocked adjacency matrix - Figure 1
# png('plots/blocked_adjacency.png', width=2, height=2, units='in', res=300)
# plot_blocked_matrix(network, gamma_mode, sort=FALSE)
# dev.off()

# ## Make dataset of latent positions
# nsamp <- nrow(mcmc_samples$beta)
# Z <- ldply(1:N, function(memb) {
  # all <- data.frame(mcmc_samples$Z[,memb,])
  # names(all) <- paste0("Z", 1:2)
  # all$node <- memb
  # all$sample <- 1:nsamp
  # all$block <- mcmc_samples$gamma[,memb]
  # return(all)
# })

# all_mean_df <- NULL
# all_edge_df <- NULL
# for(k in 1:K) {
  # nodes_in <- which(gamma_mode == k)
  # edges <- which(network[nodes_in,nodes_in]==1, arr.ind=TRUE)
  # samples_df <- subset(Z, block==k & !is.na(Z1) & node %in% nodes_in)
  # mean_df <- ddply(samples_df, .(node), function(subdf) data.frame(Z1=mean(subdf$Z1), Z2=mean(subdf$Z2)))
  # mean_df$Block <- k
  # all_mean_df <- rbind(all_mean_df, mean_df)
  
  # edge_df <- data.frame(x=mean_df$Z1[edges[,1]], xend=mean_df$Z1[edges[,2]], y=mean_df$Z2[edges[,1]], yend=mean_df$Z2[edges[,2]])
  # if(nrow(edges) != 0) edge_df$Block <- k
  # all_edge_df <- rbind(all_edge_df, edge_df)
# }
# # 
# # ggplot(all_mean_df) + geom_point(aes(Z1, Z2)) + geom_segment(data=all_edge_df, aes(x=x, xend=xend, y=y, yend=yend))+ theme_bw() + facet_wrap(~block, nrow=floor(sqrt(K)), scales='free') +coord_fixed(ratio=1)

# hhold <- read.dta('data/indian_village_raw/2. Demographics and Outcomes/household_characteristics.dta')
# hhold <- subset(hhold, village==vilno)
# hhold$hohreligion <- as.character(hhold$hohreligion)
# names(hhold)[3] <- "node"
# all_mean_df2 <- merge(hhold, all_mean_df, all.y = TRUE )
# all_mean_df2$hohreligion[is.na(all_mean_df2$hohreligion)] <- "Unknown"
# all_mean_df2$leader[is.na(all_mean_df2$leader)] <- "Unknown"

# # Plot of latent positions - Figure 3
# (g <- ggplot(all_mean_df2) + 
  # geom_segment(data=all_edge_df, aes(x=x, xend=xend, y=y, yend=yend), color='grey') + 
  # geom_point(aes(Z1, Z2, color=as.factor(leader), shape=hohreligion)) + 
  # theme_bw() +
  # facet_wrap(~Block, nrow=2, scales='free', labeller=label_both) +
  # coord_fixed(ratio=1) + 
  # scale_shape_manual(name="HH Religion", breaks=c("HINDUISM", "ISLAM", "CHRISTIANITY", "Unknown"), values=c(1, 2, 3, 4), labels=c("Hindu", "Muslim", "Christian", "Unknown")) +
  # scale_color_manual(name="HH Status", values=c("black", "red", "blue"), labels=c("Non-leader", "Leader", "Unknown")) +
  # theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
# ggsave('plots/latent_positions.png', g, width=4.5, height=2.5, units='in', scale=1.5)

# ## Posteriors of global parameters
# hpd_mu <- HPDregionplot(mcmc(mcmc_samples$mu))
# hpd <- data.frame(x=hpd_mu[[1]]$x, y=hpd_mu[[1]]$y)

# hpd_blocks <- NULL
# for(k in 1:K) {
  # hpd_k <- HPDregionplot(mcmc(cbind(mcmc_samples$beta[,k], log(mcmc_samples$sigma[,k]))))
  # for(i in 1:length(hpd_k)) hpd_blocks <- rbind(hpd_blocks, data.frame(x=hpd_k[[i]]$x, y=hpd_k[[i]]$y, Block=paste(k), group=paste(k, i)))
# }

# block_params <- ldply(1:K, function(k) data.frame(beta=mcmc_samples$beta[,k], logsigma=log(mcmc_samples$sigma[,k]), Block=k))
# mean_block_params <- ddply(block_params, .(Block), summarise, beta=mean(beta), logsigma=mean(logsigma))

# ## HPD of block-level parameters and global mean parameters - Figure 2

# (g <- ggplot(hpd_blocks) +
  # geom_path(aes(x=x, y=y, group=group)) +
  # geom_point(data=mean_block_params, aes(beta, logsigma)) +
  # xlab(expression(beta)) +
  # ylab(expression(log(sigma))) +
  # facet_wrap(~Block, labeller=label_both, nrow=2) +
  # geom_path(data=hpd, aes(x,y), linetype=2) +
  # geom_point(x=mean(mcmc_samples$beta), y=mean(log(mcmc_samples$sigma)), shape=2) +
  # theme_bw() +
  # theme(strip.background=element_blank()))
# ggsave('plots/posterior_parameters.png', g, width=4.5, height=2.5, units='in', scale=1.5)

# ## Density plots of block-specific parameters ## NOT USED IN PAPER

# ggplot(block_params) +
  # geom_density_2d(aes(beta, logsigma)) +
  # geom_point(data=mean_block_params, aes(beta, logsigma), color='red') +
  # facet_wrap(~Block, labeller=label_both, nrow=2) +
  # xlab(expression(beta)) +
  # ylab(expression(log(sigma))) +
  # theme_bw() +
  # geom_density_2d(data=data.frame(beta=mcmc_samples$mu[,1], logsigma=mcmc_samples$mu[,2]), aes(beta, logsigma), color='black') +
  # theme(strip.background=element_blank())

# ## Global mean samples and posterior density ## NOT USED IN PAPER

# ggplot(data.frame(mcmc_samples$mu)) +
  # geom_point(data=data.frame(mcmc_samples$mu[sample(1:nrow(mcmc_samples$mu), 1000, replace=FALSE),]), aes(X1, X2), alpha=.5) +
  # geom_density_2d(aes(X1, X2)) +
  # geom_path(data=hpd, aes(x,y), col='red') +
  # xlab(expression(mu[1])) +
  # ylab(expression(mu[2])) +
  # geom_point(x=mean(mcmc_samples$mu[,1]), y=mean(mcmc_samples$mu[,2]), col='red') +
  # theme_bw()

# #####

# # Trace plot of global parameters

# mu_samples <- data.frame(Iteration=rep(1:10000, 2), Value=c(mcmc_samples$mu[,1], mcmc_samples$mu[,2]), Parameter=rep(c("mu[1]", "mu[2]"), each=10000))

# (g <- ggplot(mu_samples) + 
  # geom_line(aes(Iteration, Value)) + 
  # facet_wrap(~Parameter, scales='free', labeller='label_parsed') + 
  # theme_minimal())
# ggsave('plots/mu_trace.png', g, width=5.5, height=2.5, units='in', scale=1.5)

# ## Computations presented in paper

# mean(mcmc_samples$beta[,1] < mcmc_samples$mu[,1])
# mean(mcmc_samples$beta[,2] < mcmc_samples$mu[,1])
# mean(mcmc_samples$beta[,3] < mcmc_samples$mu[,1])
# mean(mcmc_samples$beta[,4] < mcmc_samples$mu[,1])
# mean(mcmc_samples$beta[,5] < mcmc_samples$mu[,1])
# mean(mcmc_samples$beta[,6] < mcmc_samples$mu[,1])

# mean(log(mcmc_samples$sigma[,1]) < mcmc_samples$mu[,2])
# mean(log(mcmc_samples$sigma[,2]) < mcmc_samples$mu[,2])
# mean(log(mcmc_samples$sigma[,3]) < mcmc_samples$mu[,2])
# mean(log(mcmc_samples$sigma[,4]) < mcmc_samples$mu[,2])
# mean(log(mcmc_samples$sigma[,5]) < mcmc_samples$mu[,2])
# mean(log(mcmc_samples$sigma[,6]) < mcmc_samples$mu[,2])
