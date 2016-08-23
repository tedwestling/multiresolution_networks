# ###############################################################################
# # Multi-resolution blockmodel
# #
# # file: indian_village_postprocess.R
# # 
# # This file does convergence checks, does post processing, and makes plots for 
# #the India village data.  Assuming that we have already run India_village_estimation.R
# ## and have results.  
# #
# #
# #Setup here to use four chains, but can be run with more/fewer.  To use this file, you'll
# #need to be sure that the names of the 
# #
# #This file should be placed in the same directory as  the header.R file.
# #
# # Author: tedwestling
# ###############################################################################
# rm(list = ls())
# # Source header file (should be in top level of working directory)
# source('header.R')

# # We set the seed numerous times throughout the document at points at natural reload points
# set.seed(22311) 
# vilno <- 59
# edgetype <- 'visitcome'

# # Load in household-level network data for the particular edgetype and village
# network <- as.matrix(read.csv(paste0('data/indian_village_raw/1. Network Data/Adjacency Matrices/adj_', edgetype, '_HH_vilno_', vilno, '.csv'), header=FALSE))

# # Remove nodes with degree 0
# fullnetwork <- network
# zros <- which(colSums(network) == 0)
# network <- network[-zros, -zros]
# N <- nrow(network)
# K <- Khat <- 6

# thinfirst=T
# thinfac=4 #keep every fourth
# #note on these runs we've already thinned by 10 within the sampler
# ####for the longer chains, add some thinning

# load('data/results/village_59_mcmc_strongass_v3_1.Rdata')
# chain1a=chain1
# load('data/results/village_59_mcmc_strongass_v3_2.Rdata')
# chain2=chain1
# load('data/results/village_59_mcmc_strongass_v3_3.Rdata')
# chain3=chain1
# load('data/results/village_59_mcmc_strongass_v3_4.Rdata')
# chain4=chain1
# chain1=chain1a

# if(thinfirst==T){
	# #list(beta=beta, sigma=sigma, pi=pi, mu=mu, Sigma=Sigma, B=B, Z=Z, gamma=gamma)
	# keep=which(c(1:length(chain1$beta[,1]))%%thinfac==0)
	# chain1=list(beta=chain1$beta[keep,], sigma=chain1$sigma[keep,], pi=chain1$pi[keep,], mu=chain1$mu[keep,], Sigma=chain1$Sigma[keep,,], B=chain1$B[keep,,], Z=chain1$Z[keep,,], gamma=chain1$gamma[keep,])
	# #
	# chain2=list(beta=chain2$beta[keep,], sigma=chain2$sigma[keep,], pi=chain2$pi[keep,], mu=chain2$mu[keep,], Sigma=chain2$Sigma[keep,,], B=chain2$B[keep,,], Z=chain2$Z[keep,,], gamma=chain2$gamma[keep,])
	# #
	# chain3=list(beta=chain3$beta[keep,], sigma=chain3$sigma[keep,], pi=chain3$pi[keep,], mu=chain3$mu[keep,], Sigma=chain3$Sigma[keep,,], B=chain3$B[keep,,], Z=chain3$Z[keep,,], gamma=chain3$gamma[keep,])
	# #
	# chain4=list(beta=chain4$beta[keep,], sigma=chain4$sigma[keep,], pi=chain4$pi[keep,], mu=chain4$mu[keep,], Sigma=chain4$Sigma[keep,,], B=chain4$B[keep,,], Z=chain4$Z[keep,,], gamma=chain4$gamma[keep,])

# }

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
# # Use asortative spectral clustering estimate as fixed membership vector to rotate towards
# spec <- spectral_cluster(network, Krange = Khat, assortative = TRUE, plot=FALSE, degree_correct = FALSE)
# #uncomment to plot spectral clusters
# #plot_blocked_matrix(network, spec[[1]]$clusters)

# mcmc_samplesc1 <- postprocess_MCMC(chain1, network, fixed_memb = spec[[1]]$clusters)
# mcmc_samplesc2 <- postprocess_MCMC(chain2, network, fixed_memb = spec[[1]]$clusters)
# mcmc_samplesc3 <- postprocess_MCMC(chain3, network, fixed_memb = spec[[1]]$clusters)
# mcmc_samplesc4 <- postprocess_MCMC(chain4, network, fixed_memb = spec[[1]]$clusters)

# ####checkpoint save and load if above is run already
# ####checkpoint save and load if above is run already
# ####checkpoint save and load if above is run already
# #save.image("data/results/village_59_mcmc_strongass_v3_postprocessed.Rdata")
rm(list = ls())
load("data/results/village_59_mcmc_strongass_v1_postprocessed.Rdata")

#assuming here that we've done no burn-in on the sampling, so removing first 1k iterations
nburn=c(1:(length(chain1$mu[,1])*.25))

#all the samples combined together for analysis later
mcmc_samples <- list(beta=rbind(mcmc_samplesc1$beta[-nburn,], mcmc_samplesc2$beta[-nburn,],mcmc_samplesc3$beta[-nburn,], mcmc_samplesc4$beta[-nburn,]),
                sigma=rbind(mcmc_samplesc1$sigma[-nburn,], mcmc_samplesc2$sigma[-nburn,], mcmc_samplesc3$sigma[-nburn,], mcmc_samplesc4$sigma[-nburn,]),
                 pi=rbind(mcmc_samplesc1$pi[-nburn,], mcmc_samplesc2$pi[-nburn,], mcmc_samplesc3$pi[-nburn,], mcmc_samplesc4$pi[-nburn,]),
                 mu=rbind(mcmc_samplesc1$mu[-nburn,], mcmc_samplesc2$mu[-nburn,], mcmc_samplesc3$mu[-nburn,], mcmc_samplesc4$mu[-nburn,]),
                 Sigma=abind(mcmc_samplesc1$Sigma[-nburn,,], mcmc_samplesc2$Sigma[-nburn,,], mcmc_samplesc3$Sigma[-nburn,,], mcmc_samplesc4$Sigma[-nburn,,], along=1),
                 B=abind(mcmc_samplesc1$B[-nburn,,], mcmc_samplesc2$B[-nburn,,], mcmc_samplesc3$B[-nburn,,], mcmc_samplesc4$B[-nburn,,], along=1),
                 Z=abind(mcmc_samplesc1$Z[-nburn,,], mcmc_samplesc2$Z[-nburn,,], mcmc_samplesc3$Z[-nburn,,], mcmc_samplesc4$Z[-nburn,,], along=1),
                 gamma=rbind(mcmc_samplesc1$gamma[-nburn,], mcmc_samplesc2$gamma[-nburn,], mcmc_samplesc3$gamma[-nburn,], mcmc_samplesc4$gamma[-nburn,]))

#now do some mcmc checks with the "monitor" function from stan
#note that this assumes you have at least two chains
#separate arrays for analyzing convergence
#note the first argument takes the form iterations, chains, paramters
require(rstan)
#set this based on the number of chains you've run
howmanychains=4
beta_array=array(dim=c(dim(mcmc_samplesc1$beta)[1],howmanychains,dim(mcmc_samplesc1$beta)[2]))
beta_array[,1,]<-mcmc_samplesc1$beta
beta_array[,2,]<-mcmc_samplesc2$beta
beta_array[,3,]<-mcmc_samplesc3$beta
beta_array[,4,]<-mcmc_samplesc4$beta
print(monitor(beta_array))

sigma_array=array(dim=c(dim(mcmc_samplesc1$sigma)[1],howmanychains,dim(mcmc_samplesc1$sigma)[2]))
sigma_array[,1,]<-mcmc_samplesc1$sigma
sigma_array[,2,]<-mcmc_samplesc2$sigma
sigma_array[,3,]<-mcmc_samplesc3$sigma
sigma_array[,4,]<-mcmc_samplesc4$sigma
monitor(sigma_array)

pi_array=array(dim=c(dim(mcmc_samplesc1$pi)[1],howmanychains,dim(mcmc_samplesc1$pi)[2]))
pi_array[,1,]<-mcmc_samplesc1$pi
pi_array[,2,]<-mcmc_samplesc2$pi
pi_array[,3,]<-mcmc_samplesc3$pi
pi_array[,4,]<-mcmc_samplesc4$pi
monitor(pi_array)


mu_array=array(dim=c(dim(mcmc_samplesc1$mu)[1],howmanychains,dim(mcmc_samplesc1$mu)[2]))
mu_array[,1,]<-mcmc_samplesc1$mu
mu_array[,2,]<-mcmc_samplesc2$mu
mu_array[,3,]<-mcmc_samplesc3$mu
mu_array[,4,]<-mcmc_samplesc4$mu
print(monitor(mu_array))

Sigma_array=array(dim=c(dim(mcmc_samplesc1$Sigma)[1],howmanychains,4))
Sigma_array[,1,]<-matrix(mcmc_samplesc1$Sigma,dim(mcmc_samplesc1$Sigma)[1],4)
Sigma_array[,2,]<-matrix(mcmc_samplesc2$Sigma,dim(mcmc_samplesc2$Sigma)[1],4)
Sigma_array[,3,]<-matrix(mcmc_samplesc3$Sigma,dim(mcmc_samplesc3$Sigma)[1],4)
Sigma_array[,4,]<-matrix(mcmc_samplesc4$Sigma,dim(mcmc_samplesc4$Sigma)[1],4)
monitor(Sigma_array)

#note that 30 here is set up for 6 blocks (6*6-6 on the diagonal)
#need to change if you're using fewer blocks
B_array=array(dim=c(dim(mcmc_samplesc1$B)[1],howmanychains,((Khat*Khat)-Khat)))
B_array[,1,]<-matrix(mcmc_samplesc1$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc1$B)[1],((Khat*Khat)-Khat))
B_array[,2,]<-matrix(mcmc_samplesc2$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc2$B)[1],((Khat*Khat)-Khat))
B_array[,3,]<-matrix(mcmc_samplesc3$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc3$B)[1],((Khat*Khat)-Khat))
B_array[,4,]<-matrix(mcmc_samplesc4$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc4$B)[1],((Khat*Khat)-Khat))
monitor(B_array)

gamma_array=array(dim=c(dim(mcmc_samplesc1$gamma)[1],howmanychains,dim(mcmc_samplesc1$gamma)[2]))
gamma_array[,1,]<-mcmc_samplesc1$gamma
gamma_array[,2,]<-mcmc_samplesc2$gamma
gamma_array[,3,]<-mcmc_samplesc3$gamma
gamma_array[,4,]<-mcmc_samplesc4$gamma
monitor(gamma_array)





prefix="data/results/village_59_mcmc_strongass_k5v1_"
###make trace plots
#mu
pdf(file=paste(prefix,"mutrace.pdf",sep=''),height=10,width=6)
par(mfrow=c(2,1))
matplot(mu_array[-c(1:max(nburn)),,1],type='l',main="Traceplot for mu1",xlab="iteration",ylab="mu1")
matplot(mu_array[-c(1:max(nburn)),,2],type='l',main="Traceplot for mu2",xlab="iteration",ylab="mu2")
dev.off()



pdf(file=paste(prefix,"betatrace.pdf",sep=''),height=10,width=6)
par(mfrow=c(3,2))
matplot(beta_array[-c(1:max(nburn)),,1],type='l',main="Traceplot for beta1",xlab="iteration",ylab="beta1")
matplot(beta_array[-c(1:max(nburn)),,2],type='l',main="Traceplot for beta2",xlab="iteration",ylab="beta2")
matplot(beta_array[-c(1:max(nburn)),,3],type='l',main="Traceplot for beta3",xlab="iteration",ylab="beta3")
matplot(beta_array[-c(1:max(nburn)),,4],type='l',main="Traceplot for beta4",xlab="iteration",ylab="beta4")
matplot(beta_array[-c(1:max(nburn)),,5],type='l',main="Traceplot for beta5",xlab="iteration",ylab="beta5")
if(Khat==6){matplot(beta_array[-c(1:max(nburn)),,6],type='l',main="Traceplot for beta6",xlab="iteration",ylab="beta6")}
dev.off()

pdf(file=paste(prefix,"sigmatrace.pdf",sep=''),height=10,width=6)
par(mfrow=c(3,2))
matplot(sigma_array[-c(1:max(nburn)),,1],type='l',main="Traceplot for sigma1",xlab="iteration",ylab="sigma1")
matplot(sigma_array[-c(1:max(nburn)),,2],type='l',main="Traceplot for sigma2",xlab="iteration",ylab="sigma2")
matplot(sigma_array[-c(1:max(nburn)),,3],type='l',main="Traceplot for sigma3",xlab="iteration",ylab="sigma3")
matplot(sigma_array[-c(1:max(nburn)),,4],type='l',main="Traceplot for sigma4",xlab="iteration",ylab="sigma4")
matplot(sigma_array[-c(1:max(nburn)),,5],type='l',main="Traceplot for sigma5",xlab="iteration",ylab="sigma5")
if(Khat==6){matplot(sigma_array[-c(1:max(nburn)),,6],type='l',main="Traceplot for sigma6",xlab="iteration",ylab="sigma6")}
dev.off()


pdf(file=paste(prefix,"pitrace.pdf",sep=''),height=10,width=6)
par(mfrow=c(3,2))
matplot(pi_array[-c(1:max(nburn)),,1],type='l',main="Traceplot for pi1",xlab="iteration",ylab="pi1")
matplot(pi_array[-c(1:max(nburn)),,2],type='l',main="Traceplot for pi2",xlab="iteration",ylab="pi2")
matplot(pi_array[-c(1:max(nburn)),,3],type='l',main="Traceplot for pi3",xlab="iteration",ylab="pi3")
matplot(pi_array[-c(1:max(nburn)),,4],type='l',main="Traceplot for pi4",xlab="iteration",ylab="pi4")
matplot(pi_array[-c(1:max(nburn)),,5],type='l',main="Traceplot for pi5",xlab="iteration",ylab="pi5")
if(Khat==6){matplot(pi_array[-c(1:max(nburn)),,6],type='l',main="Traceplot for pi6",xlab="iteration",ylab="pi6")}
dev.off()


pdf(file=paste(prefix,"Sigmatrace.pdf",sep=''),height=10,width=6)
par(mfrow=c(2,2))
matplot(Sigma_array[-c(1:max(nburn)),,1],type='l',main="Traceplot for Sigma1",xlab="iteration",ylab="Sigma1")
matplot(Sigma_array[-c(1:max(nburn)),,2],type='l',main="Traceplot for Sigma2",xlab="iteration",ylab="Sigma2")
matplot(Sigma_array[-c(1:max(nburn)),,3],type='l',main="Traceplot for Sigma3",xlab="iteration",ylab="Sigma3")
matplot(Sigma_array[-c(1:max(nburn)),,4],type='l',main="Traceplot for Sigma4",xlab="iteration",ylab="Sigma4")
dev.off()

set.seed(47676728)


# CODE FOR COMPUTING JOINT POSTERIOR MODE:
gamma_vecs <- apply(mcmc_samples$gamma,1,paste0, collapse='.')
gamma_tab <- table(gamma_vecs)
head(sort(gamma_tab))
# If the maximum of gamma_tab is 1 then no membership appeared more than once and the marginal mode should be used instead
gamma_mode <- names(gamma_tab)[which.max(gamma_tab)]
joint_mode <- as.numeric(strsplit(gamma_mode, ".", fixed=TRUE)[[1]])

# CODE FOR COMPUTING MARGINAL POSTERIOR MODES:
marg_mode <- as.numeric(apply(mcmc_samples$gamma, 2, function(memb) names(table(memb))[which.max(table(memb))]))

# CHOOSE ONE FOR DISPLAYING POSITIONS
gamma_mode <- marg_mode

# Plot blocked adjacency matrix - Figure 1
png(paste(prefix,'blocked_adjacency.png',sep=''), width=2, height=2, units='in', res=300)
plot_blocked_matrix(network, gamma_mode, sort=FALSE)
dev.off()

## Make dataset of latent positions
nsamp <- nrow(mcmc_samples$beta)
Z <- ldply(1:N, function(memb) {
  all <- data.frame(mcmc_samples$Z[,memb,])
  names(all) <- paste0("Z", 1:2)
  all$node <- memb
  all$sample <- 1:nsamp
  all$block <- mcmc_samples$gamma[,memb]
  return(all)
})

all_mean_df <- NULL
all_edge_df <- NULL
for(k in 1:K) {
  nodes_in <- which(gamma_mode == k)
  edges <- which(network[nodes_in,nodes_in]==1, arr.ind=TRUE)
  samples_df <- subset(Z, block==k & !is.na(Z1) & node %in% nodes_in)
  mean_df <- ddply(samples_df, .(node), function(subdf) data.frame(Z1=mean(subdf$Z1), Z2=mean(subdf$Z2)))
  mean_df$Block <- k
  all_mean_df <- rbind(all_mean_df, mean_df)
  
  edge_df <- data.frame(x=mean_df$Z1[edges[,1]], xend=mean_df$Z1[edges[,2]], y=mean_df$Z2[edges[,1]], yend=mean_df$Z2[edges[,2]])
  if(nrow(edges) != 0) edge_df$Block <- k
  all_edge_df <- rbind(all_edge_df, edge_df)
}
# 
# ggplot(all_mean_df) + geom_point(aes(Z1, Z2)) + geom_segment(data=all_edge_df, aes(x=x, xend=xend, y=y, yend=yend))+ theme_bw() + facet_wrap(~block, nrow=floor(sqrt(K)), scales='free') +coord_fixed(ratio=1)

hhold <- read.dta('data/indian_village_raw/2. Demographics and Outcomes/household_characteristics.dta')
hhold <- subset(hhold, village==vilno)
hhold$hohreligion <- as.character(hhold$hohreligion)
names(hhold)[3] <- "node"
all_mean_df2 <- merge(hhold, all_mean_df, all.y = TRUE )
all_mean_df2$hohreligion[is.na(all_mean_df2$hohreligion)] <- "Unknown"
all_mean_df2$leader[is.na(all_mean_df2$leader)] <- "Unknown"

# Plot of latent positions - Figure 3
(g <- ggplot(all_mean_df2) + 
  geom_segment(data=all_edge_df, aes(x=x, xend=xend, y=y, yend=yend), color='grey') + 
  geom_point(aes(Z1, Z2, color=as.factor(leader), shape=hohreligion)) + 
  theme_bw() +
  facet_wrap(~Block, nrow=2, scales='free') +
  coord_fixed(ratio=1) + 
  scale_shape_manual(name="HH Religion", breaks=c("HINDUISM", "ISLAM", "CHRISTIANITY", "Unknown"), values=c(1, 2, 3, 4), labels=c("Hindu", "Muslim", "Christian", "Unknown")) +
  scale_color_manual(name="HH Status", values=c("black", "red", "blue"), labels=c("Non-leader", "Leader", "Unknown")) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
ggsave(paste(prefix,'latent_positions.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Posteriors of global parameters
hpd_mu <- HPDregionplot(mcmc(mcmc_samples$mu))
hpd <- data.frame(x=hpd_mu[[1]]$x, y=hpd_mu[[1]]$y)

hpd_blocks <- NULL
for(k in 1:K) {
  hpd_k <- HPDregionplot(mcmc(cbind(mcmc_samples$beta[,k], log(mcmc_samples$sigma[,k]))))
  for(i in 1:length(hpd_k)) hpd_blocks <- rbind(hpd_blocks, data.frame(x=hpd_k[[i]]$x, y=hpd_k[[i]]$y, Block=paste(k), group=paste(k, i)))
}

block_params <- ldply(1:K, function(k) data.frame(beta=mcmc_samples$beta[,k], logsigma=log(mcmc_samples$sigma[,k]), Block=k))
mean_block_params <- ddply(block_params, .(Block), summarise, beta=mean(beta), logsigma=mean(logsigma))

## HPD of block-level parameters and global mean parameters - Figure 2

(g <- ggplot(hpd_blocks) +
  geom_path(aes(x=x, y=y, group=group)) +
  geom_point(data=mean_block_params, aes(beta, logsigma)) +
  xlab(expression(beta)) +
  ylab(expression(log(sigma))) +
  facet_wrap(~Block, nrow=2) +
  geom_path(data=hpd, aes(x,y), linetype=2) +
  geom_point(x=mean(mcmc_samples$beta), y=mean(log(mcmc_samples$sigma)), shape=2) +
  theme_bw() +
  theme(strip.background=element_blank()))
ggsave(paste(prefix,'posterior_parameters.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Density plots of block-specific parameters ## NOT USED IN PAPER

ggplot(block_params) +
  geom_density2d(aes(beta, logsigma)) +
  geom_point(data=mean_block_params, aes(beta, logsigma), color='red') +
  facet_wrap(~Block, nrow=2) +
  xlab(expression(beta)) +
  ylab(expression(log(sigma))) +
  theme_bw() +
  geom_density2d(data=data.frame(beta=mcmc_samples$mu[,1], logsigma=mcmc_samples$mu[,2]), aes(beta, logsigma), color='black') +
  theme(strip.background=element_blank())
ggsave(paste(prefix,'posterior_parameters.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Global mean samples and posterior density ## NOT USED IN PAPER

ggplot(data.frame(mcmc_samples$mu)) +
  geom_point(data=data.frame(mcmc_samples$mu[sample(1:nrow(mcmc_samples$mu), 1000, replace=FALSE),]), aes(X1, X2), alpha=.5) +
  geom_density2d(aes(X1, X2)) +
  geom_path(data=hpd, aes(x,y), col='red') +
  xlab(expression(mu[1])) +
  ylab(expression(mu[2])) +
  geom_point(x=mean(mcmc_samples$mu[,1]), y=mean(mcmc_samples$mu[,2]), col='red') +
  theme_bw()
ggsave(paste(prefix,'posterior_parameters2.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Computations presented in paper

mean(mcmc_samples$beta[,1] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,2] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,3] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,4] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,5] < mcmc_samples$mu[,1])
if(Khat==6){mean(mcmc_samples$beta[,6] < mcmc_samples$mu[,1])}

mean(log(mcmc_samples$sigma[,1]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,2]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,3]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,4]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,5]) < mcmc_samples$mu[,2])
if(Khat==6){mean(log(mcmc_samples$sigma[,6]) < mcmc_samples$mu[,2])}
