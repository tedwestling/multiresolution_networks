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
#setwd("")
setwd("/Users/tylermccormick/Dropbox/git_to_work/archive/multiresolution_networks")
source("header.R")
load("data/results/village_59_mcmc_strongass_v1_postprocessed.Rdata")

require(abind)
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

#block_latent_MCMC <- function(Y, D, K, burn_in, n_samples, thin, v, epsilon=.1, rZ=1, Atheta=matrix(c(2,1,1,1), nrow=2), alpha0=.1, a0=NULL, b0=NULL, m0=c(0,0), s0=0.01, psi0=matrix(c(5.1,0,0,5.1), nrow=2), nu0=5.1, verbose=TRUE, memb_start=NULL, plot_init=FALSE, likelihood=FALSE, true_gamma=NULL, postprocess=TRUE, sample_membs=TRUE, record_acc_probs=FALSE, debug_output=FALSE, perturb_init=TRUE) {
	
	#log_likelihood=function(Y, par_t, alpha0, a0, b0, m0, s0, psi0, nu0) {
	
log_likelihood <- function(Y, par_t, alpha0, a0, b0, m0, s0, psi0, nu0) {
  K <- length(par_t$pi)
  d <- ncol(par_t$Z)
  cD <- 2 * gamma((d+1)/2) / gamma(d/2)
  pi_ll <- sum(par_t$block_n * log(par_t$pi))
  Yoff_ll <- sum(par_t$s[lower.tri(par_t$s, diag=FALSE)] * log(par_t$B[lower.tri(par_t$B, diag=FALSE)]))
  theta <- cbind(as.vector(par_t$beta), as.vector(log(par_t$sigma)))
  theta_ll <- - K * log(2*3.1415)- (K/2) * log(det(par_t$Sigma)) + sum(apply(theta, 1, function(thetak) - .5 * (thetak - par_t$mu) %*% solve(par_t$Sigma) %*% t(thetak - par_t$mu)))
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
  muprior_ll <- - log(2*3.1415) - (1/2) * log(det(par_t$Sigma/ s0)) - (1/2) * (par_t$mu - m0) %*% solve(par_t$Sigma / s0) %*% t(par_t$mu - m0)
  Sigmaprior_ll <- log(diwish(par_t$Sigma, nu0, psi0))
  component_lls=c(Ydiag_ll=Ydiag_ll, Yoff_ll=Yoff_ll, pi_ll=pi_ll, z_ll=z_ll, theta_ll=theta_ll, piprior_ll=piprior_ll, Bprior_ll=Bprior_ll, muprior_ll=muprior_ll, Sigmaprior_ll=Sigmaprior_ll)
  sum(component_lls)
}

	
dens <- sum(network, na.rm=TRUE)/(N^2 - N)
b0=1
a0=(b0*(10*dens))/(1-(10*dens))

postprob=rep(NA,nrow(mcmc_samples$beta))
for(ii in 1:nrow(mcmc_samples$beta)){
	mcmc_samples_tmp <- list(beta=rbind(mcmc_samples$beta[ii,]),
                sigma=rbind(mcmc_samples$sigma[ii,]),
                 pi=rbind(mcmc_samples$pi[ii,]),
                 mu=rbind(mcmc_samples$mu[ii,]),
                 Sigma=abind(mcmc_samples$Sigma[ii,,], along=1),
                 B=abind(mcmc_samples$B[ii,,],along=1),
                 Z=abind(mcmc_samples$Z[ii,,], along=1),
                 gamma=rbind(mcmc_samples$gamma[ii,]))
	par_t<-mcmc_samples_tmp
if(ii%%100==0){print(ii)}
postprob[ii]=log_likelihood(network,mcmc_samples_tmp,alpha0=.1, a0=a0, b0=b0, m0=c(0,0), s0=0.01, psi0=matrix(c(5.1,0,0,5.1), nrow=2), nu0=5.1)
}

postdens=density(postprob)
postmode=which.max(postdens$y)

mcmc_samplespm=list(beta=rbind(mcmc_samples$beta[postmode,]),
                sigma=rbind(mcmc_samples$sigma[postmode,]),
                 pi=rbind(mcmc_samples$pi[postmode,]),
                 mu=rbind(mcmc_samples$mu[postmode,]),
                 Sigma=abind(mcmc_samples$Sigma[postmode,,], along=1),
                 B=abind(mcmc_samples$B[postmode,,],along=1),
                 Z=abind(mcmc_samples$Z[postmode,,], along=1),
                 gamma=rbind(mcmc_samples$gamma[postmode,]))
	
all_mean_df<-cbind(c(1:length(mcmc_samplespm$gamma)),mcmc_samplespm$Z,as.vector(mcmc_samplespm$gamma),rep(1,length(mcmc_samplespm$gamma)))
colnames(all_mean_df)<-c("node","Z1","Z2","Block","prob")
Z=data.frame(all_mean_df)
	
hhold<-read.dta("/Users/tylermccormick/Dropbox/git_to_work/net_comp_thry/data/indian_village_raw/2. Demographics and Outcomes/household_characteristics.dta")
hhold <- subset(hhold, village==vilno)
hhold$castesubcaste <- as.character(hhold$castesubcaste)
names(hhold)[3] <- "node"
#merge the covariate information with the block and latent positions
all_mean_df2 <- merge(hhold, all_mean_df, all.y = TRUE )
all_mean_df2$castesubcaste[is.na(all_mean_df2$castesubcaste)] <- "Unknown"
all_mean_df2$leader[is.na(all_mean_df2$leader)] <- "Unknown"
all_mean_df2_tr=all_mean_df2

all_edge_df_tr <- NULL
for(k in 1:K) {
  nodes_in <- unique(all_mean_df2_tr$node[all_mean_df2_tr$Block == k])
  edges <- which(network[nodes_in,nodes_in]==1, arr.ind=TRUE)
  samples_df <- subset(Z, Block==k & !is.na(Z1) & node %in% nodes_in)
#edge_df_tr <- data.frame(x=all_mean_df2_tr$Z1[edges[,1]], xend=all_mean_df2_tr$Z1[edges[,2]], y=all_mean_df2_tr$Z2[edges[,1]], yend=all_mean_df2_tr$Z2[edges[,2]], prob=(all_mean_df2_tr$prob[edges[,1]]+all_mean_df2_tr$prob[edges[,2]])/2)
	mean_df <- ddply(samples_df, .(node), function(subdf) data.frame(Z1=mean(subdf$Z1), Z2=mean(subdf$Z2), prob=nrow(subdf)))
 edge_df_tr <- data.frame(x=mean_df$Z1[edges[,1]], xend=mean_df$Z1[edges[,2]], y=mean_df$Z2[edges[,1]], yend=mean_df$Z2[edges[,2]], prob=(mean_df$prob[edges[,1]]+mean_df$prob[edges[,2]])/2,send=edges[,1],rec=edges[,2])
  if(nrow(edges) != 0) edge_df_tr$Block <- k
  all_edge_df_tr <- rbind(all_edge_df_tr, edge_df_tr)
}

block_names <- list(
  '1'="Block 1",
  '2'="Block 2",
  '3'="Block 3",
  '4'="Block 4",
    '5'="Block 5",
  '6'="Block 6"
)
block_labeller <- function(variable,value){
  return(block_names[value])
}

prefix="/Users/tylermccormick/Dropbox/git_to_work/archive/multiresolution_networks"

(g<-ggplot(all_mean_df2_tr) + 
  geom_segment(data=all_edge_df_tr, aes(x=x, xend=xend, y=y, yend=yend, alpha=prob/2)) +
  guides(alpha=FALSE)+
  geom_point(aes(Z1, Z2, size=prob,color=as.factor(castesubcaste),shape=as.factor(castesubcaste))) +
	scale_size_continuous(range = c(1,4))+
  #geom_point(aes(Z1, Z2, size=prob,color=as.factor(castesubcaste),shape=bshape)) +
  guides(cex=FALSE)+ 
  geom_vline(xintercept=0,linetype=2,color="grey64")+
  geom_hline(yintercept=0,linetype=2,color="grey64")+
  #geom_text(aes(label=node, x=Z1, y=Z2, size=prob)) + 
  theme_bw() +
  facet_wrap(~Block, nrow=2, scales='free',labeller=block_labeller) +
  #labs(x="First latent dimension",y="Second latent dimension")+
  coord_fixed(ratio=1) + 
 scale_color_manual(name="HH Caste", values=c(1:6), labels=c("General", "Minority", "OBC", "Schedule caste", "Schedule tribe", "Unknown")) +
 scale_shape_manual(name="HH Caste", values=c(15:20), labels=c("General", "Minority", "OBC", "Schedule caste", "Schedule tribe", "Unknown")) +
theme(axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
#theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), strip.background=element_blank()))
ggsave(paste(prefix,'latent_positions_shadedv2_jointpost.png',sep=''), g, width=6.5, height=3, units='in', scale=1.5)
######################################################shaded by probabilty
######################################################end of changes



