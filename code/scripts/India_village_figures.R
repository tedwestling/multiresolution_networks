# ###############################################################################
# # Multi-resolution blockmodel
# #
# # file: indian_village_figures.R
# # 
# # This file makes the plots for the paper.  Assuming that we have already run India_village_estimation.R
# ## and have results and have post-processed them with India_village_postprocess.R.  
# #
# #
# # Author: tedwestling
# ###############################################################################
 rm(list = ls())
#set your working directory to the top level multiresolution_networks folder
setwd("")
source("header.R")
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

sigma_array=array(dim=c(dim(mcmc_samplesc1$sigma)[1],howmanychains,dim(mcmc_samplesc1$sigma)[2]))
sigma_array[,1,]<-mcmc_samplesc1$sigma
sigma_array[,2,]<-mcmc_samplesc2$sigma
sigma_array[,3,]<-mcmc_samplesc3$sigma
sigma_array[,4,]<-mcmc_samplesc4$sigma


pi_array=array(dim=c(dim(mcmc_samplesc1$pi)[1],howmanychains,dim(mcmc_samplesc1$pi)[2]))
pi_array[,1,]<-mcmc_samplesc1$pi
pi_array[,2,]<-mcmc_samplesc2$pi
pi_array[,3,]<-mcmc_samplesc3$pi
pi_array[,4,]<-mcmc_samplesc4$pi

mu_array=array(dim=c(dim(mcmc_samplesc1$mu)[1],howmanychains,dim(mcmc_samplesc1$mu)[2]))
mu_array[,1,]<-mcmc_samplesc1$mu
mu_array[,2,]<-mcmc_samplesc2$mu
mu_array[,3,]<-mcmc_samplesc3$mu
mu_array[,4,]<-mcmc_samplesc4$mu

Sigma_array=array(dim=c(dim(mcmc_samplesc1$Sigma)[1],howmanychains,4))
Sigma_array[,1,]<-matrix(mcmc_samplesc1$Sigma,dim(mcmc_samplesc1$Sigma)[1],4)
Sigma_array[,2,]<-matrix(mcmc_samplesc2$Sigma,dim(mcmc_samplesc2$Sigma)[1],4)
Sigma_array[,3,]<-matrix(mcmc_samplesc3$Sigma,dim(mcmc_samplesc3$Sigma)[1],4)
Sigma_array[,4,]<-matrix(mcmc_samplesc4$Sigma,dim(mcmc_samplesc4$Sigma)[1],4)



#note there are 30 for 6 blocks (6*6-6 on the diagonal)
B_array=array(dim=c(dim(mcmc_samplesc1$B)[1],howmanychains,((Khat*Khat)-Khat)))
B_array[,1,]<-matrix(mcmc_samplesc1$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc1$B)[1],((Khat*Khat)-Khat))
B_array[,2,]<-matrix(mcmc_samplesc2$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc2$B)[1],((Khat*Khat)-Khat))
B_array[,3,]<-matrix(mcmc_samplesc3$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc3$B)[1],((Khat*Khat)-Khat))
B_array[,4,]<-matrix(mcmc_samplesc4$B[is.na(mcmc_samplesc1$B)==F],dim(mcmc_samplesc4$B)[1],((Khat*Khat)-Khat))

gamma_array=array(dim=c(dim(mcmc_samplesc1$gamma)[1],howmanychains,dim(mcmc_samplesc1$gamma)[2]))
gamma_array[,1,]<-mcmc_samplesc1$gamma
gamma_array[,2,]<-mcmc_samplesc2$gamma
gamma_array[,3,]<-mcmc_samplesc3$gamma
gamma_array[,4,]<-mcmc_samplesc4$gamma

#uncomment to see monitor results
#comment to supress if just re-making plots
#monitor(B_array)
#monitor(Sigma_array)
#monitor(mu_array)
#monitor(pi_array)
#monitor(sigma_array)
#monitor(beta_array)

prefix="plots/"
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



hhold <- read.dta('data/indian_village_raw/2. Demographics and Outcomes/household_characteristics.dta')
hhold <- subset(hhold, village==vilno)
######################################################shaded by probabilty within each block
nsamp <- nrow(mcmc_samples$beta)
Z <- ldply(1:N, function(memb) {
  tab <- table(mcmc_samples$gamma[,memb])/nsamp
  this_blocks <- as.numeric(names(tab)[tab > .00017]) ### has appeared in a block at least twice .00017
#print(this_blocks)
ldply(this_blocks, function(block) {
    df <- data.frame(mcmc_samples$Z[mcmc_samples$gamma[,memb] == block,memb,])
    names(df) <- paste0("Z", 1:2)
    df$block <- block
    df$node <- memb
    return(df)
  })
})

all_mean_df <- NULL
all_edge_df <- NULL
for(k in 1:K) {
  nodes_in <- unique(Z$node[Z$block == k])
  edges <- which(network[nodes_in,nodes_in]==1, arr.ind=TRUE)
  samples_df <- subset(Z, block==k & !is.na(Z1) & node %in% nodes_in)
####in the next line ted's code had 10000, changed to nsamp.  NEED TO DOUBLE CHECK THIS!!
  mean_df <- ddply(samples_df, .(node), function(subdf) data.frame(Z1=mean(subdf$Z1), Z2=mean(subdf$Z2), prob=nrow(subdf)/nsamp))
  mean_df$Block <- k
  all_mean_df <- rbind(all_mean_df, mean_df)
  edge_df <- data.frame(x=mean_df$Z1[edges[,1]], xend=mean_df$Z1[edges[,2]], y=mean_df$Z2[edges[,1]], yend=mean_df$Z2[edges[,2]], prob=(mean_df$prob[edges[,1]]+mean_df$prob[edges[,2]])/2,send=edges[,1],rec=edges[,2])
  if(nrow(edges) != 0) edge_df$Block <- k
  all_edge_df <- rbind(all_edge_df, edge_df)
}

###set up for plot by caste
hhold$castesubcaste <- as.character(hhold$castesubcaste)
names(hhold)[3] <- "node"
#merge the covariate information with the block and latent positions
all_mean_df2 <- merge(hhold, all_mean_df, all.y = TRUE )
all_mean_df2$castesubcaste[is.na(all_mean_df2$castesubcaste)] <- "Unknown"
all_mean_df2$leader[is.na(all_mean_df2$leader)] <- "Unknown"

###truncate to only plot nodes with >40pct post prob of being in block
all_mean_df2_tr=all_mean_df2[all_mean_df2$prob>0.4,]
all_edge_df_tr <- NULL
for(k in 1:K) {
  nodes_in <- unique(all_mean_df2_tr$node[all_mean_df2_tr$Block == k])
  edges <- which(network[nodes_in,nodes_in]==1, arr.ind=TRUE)
  samples_df <- subset(Z, block==k & !is.na(Z1) & node %in% nodes_in)
#edge_df_tr <- data.frame(x=all_mean_df2_tr$Z1[edges[,1]], xend=all_mean_df2_tr$Z1[edges[,2]], y=all_mean_df2_tr$Z2[edges[,1]], yend=all_mean_df2_tr$Z2[edges[,2]], prob=(all_mean_df2_tr$prob[edges[,1]]+all_mean_df2_tr$prob[edges[,2]])/2)
	mean_df <- ddply(samples_df, .(node), function(subdf) data.frame(Z1=mean(subdf$Z1), Z2=mean(subdf$Z2), prob=nrow(subdf)/nsamp))
 edge_df_tr <- data.frame(x=mean_df$Z1[edges[,1]], xend=mean_df$Z1[edges[,2]], y=mean_df$Z2[edges[,1]], yend=mean_df$Z2[edges[,2]], prob=(mean_df$prob[edges[,1]]+mean_df$prob[edges[,2]])/2,send=edges[,1],rec=edges[,2])
  if(nrow(edges) != 0) edge_df_tr$Block <- k
  all_edge_df_tr <- rbind(all_edge_df_tr, edge_df_tr)
}

#set up indicators of how many blocks a person is in
bshape=rep(NA,nrow(all_mean_df2_tr))
nshows=table(all_mean_df2_tr$node)
for(nn in 1:length(bshape)){
	ntmp=all_mean_df2_tr$node[nn]
	bshape[nn]<-nshows[names(nshows)==as.character(ntmp)]
}


all_mean_df2_tr<-cbind(all_mean_df2_tr,as.factor(bshape))
names(all_mean_df2_tr)[24]<-"bshape"
rm(bshape)

(g<-ggplot(all_mean_df2_tr) + 
  geom_segment(data=all_edge_df_tr, aes(x=x, xend=xend, y=y, yend=yend, alpha=prob/2)) +
  guides(alpha=FALSE)+
  geom_point(aes(Z1, Z2, alpha=prob,color=as.factor(castesubcaste),shape=bshape)) +
  guides(cex=FALSE)+ 
  #geom_text(aes(label=node, x=Z1, y=Z2, size=prob)) + 
  theme_bw() +
  facet_wrap(~Block, nrow=2, scales='free') +
  #labs(x="First latent dimension",y="Second latent dimension")+
  coord_fixed(ratio=1) + 
 scale_color_manual(name="HH Caste", values=c(1:6), labels=c("General", "Minority", "OBC", "Schedule caste", "Schedule tribe", "Unknown")) +
 scale_shape_manual(name="Block inclusion", values=c(16,17), labels=c("Single block", "Two Blocks"))+
theme(axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
#theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
ggsave(paste(prefix,'latent_positions_shaded.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)
######################################################shaded by probabilty



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
  geom_point(x=mean(mcmc_samples$mu[,1]), y=mean((mcmc_samples$mu[,2])), shape=2) +
  theme_bw() +
  theme(strip.background=element_blank()))

ggsave(paste(prefix,'posterior_parameters.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Density plots of block-specific parameters ## NOT USED IN PAPER

# ggplot(block_params) +
  # geom_density2d(aes(beta, logsigma)) +
  # geom_point(data=mean_block_params, aes(beta, logsigma), color='red') +
  # facet_wrap(~Block, nrow=2) +
  # xlab(expression(beta)) +
  # ylab(expression(log(sigma))) +
  # theme_bw() +
  # geom_density2d(data=data.frame(beta=mcmc_samples$mu[,1], logsigma=mcmc_samples$mu[,2]), aes(beta, logsigma), color='black') +
  # theme(strip.background=element_blank())
# ggsave(paste(prefix,'posterior_parameters2.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

# ## Global mean samples and posterior density ## NOT USED IN PAPER

# ggplot(data.frame(mcmc_samples$mu)) +
  # geom_point(data=data.frame(mcmc_samples$mu[sample(1:nrow(mcmc_samples$mu), 1000, replace=FALSE),]), aes(X1, X2), alpha=.5) +
  # geom_density2d(aes(X1, X2)) +
  # geom_path(data=hpd, aes(x,y), col='red') +
  # xlab(expression(mu[1])) +
  # ylab(expression(mu[2])) +
  # geom_point(x=mean(mcmc_samples$mu[,1]), y=mean(mcmc_samples$mu[,2]), col='red') +
  # theme_bw()
# ggsave(paste(prefix,'posterior_parameters3.png',sep=''), g, width=4.5, height=2.5, units='in', scale=1.5)

## Computations presented in paper

mean(mcmc_samples$beta[,1] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,2] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,3] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,4] < mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,5] < mcmc_samples$mu[,1])
if(Khat==6){mean(mcmc_samples$beta[,6] < mcmc_samples$mu[,1])}

mean(mcmc_samples$beta[,1] > mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,2] > mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,3] > mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,4] > mcmc_samples$mu[,1])
mean(mcmc_samples$beta[,5] > mcmc_samples$mu[,1])
if(Khat==6){mean(mcmc_samples$beta[,6] > mcmc_samples$mu[,1])}

mean(log(mcmc_samples$sigma[,1]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,2]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,3]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,4]) < mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,5]) < mcmc_samples$mu[,2])
if(Khat==6){mean(log(mcmc_samples$sigma[,6]) < mcmc_samples$mu[,2])}

mean(log(mcmc_samples$sigma[,1]) > mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,2]) > mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,3]) > mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,4]) > mcmc_samples$mu[,2])
mean(log(mcmc_samples$sigma[,5]) > mcmc_samples$mu[,2])
if(Khat==6){mean(log(mcmc_samples$sigma[,6]) > mcmc_samples$mu[,2])}

##calculations for the between block and intercepts matrix
#intercepts for each block
round(apply(beta_array,3,quantile,probs=.975),3)
round(apply(beta_array,3,mean),3)

#between block probabilities
#B_array
Bmeans=round(apply(B_array,3,mean),3)
Bsd=round(apply(B_array,3,sd),3)
round(apply(B_array,3,quantile,probs=.975),3)
round(apply(B_array,3,quantile,probs=.025),3)


