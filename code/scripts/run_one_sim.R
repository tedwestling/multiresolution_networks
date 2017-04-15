###############################################################################
# Multi-resolution blockmodel
#
# file: run_one_sim.R

# This file generates a single network with n nodes (command line argument)
# and fits the LS, LPCM, SBM, and LS-SBM to the network with 10% of node pairs
# randomly held out. It then calculates and return the out-of-sample AUC and 
# AUPRC for each of the models overall, within, and betwen blocks.
# 
# See simulations_readme.txt for instructions for running the simulations
# 
# Author: tedwestling
###############################################################################

task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(batch)
parseCommandArgs()

source('header.R')
library(cvAUC)
library(network)
library(latentnet)
#n <- 100 - passed in via command line if run on cluster

# function for calculating AUPRC
auprc <- function(pred, truth) {
  vals <- sort(unique(c(0,pred, 1)))
  pr <- t(sapply(vals, function(val) {
    binpred <- as.numeric(pred >= val)
    tp <- sum(binpred == 1 & truth == 1)
    precision <- tp / sum(binpred == 1)
    recall <- tp / sum(truth == 1)
    c(recall=recall, precision=precision)
  }))
  pr_ord <- order(pr[,1])
  pr <- pr[pr_ord,]
  pr[1,2] <- 1
  rec <- sort(unique(pr[,1]))
  prec <- sapply(rec, function(r) max(pr[pr[,1] == r, 2]))
  auprc <- sum(diff(rec) * prec[-1])
  return(auprc)
}

# function for calculating overall tie probability from LS parameters
logit <- function(x) log(x) - log(1-x)
expit <- function(x) exp(x) / (1 + exp(x))
ls_mean <- function(beta, sigma) {
  f <- function(x) expit(beta - sqrt(2 * sigma^2 * x)) * dchisq(x, df=1)
  integrate(f, lower=0, upper=Inf)$value
}


# Simulation settings
K <- 5

B <- matrix(c(NA, rep(.02, 4),
              .02, NA, rep(.2, 3),
              .02, .2, NA, rep(.02, 2),
              .02, .2, .02, NA, rep(.2, 1),
              .02, .2, .02, .2, NA),
            nrow=K, ncol=K)

beta <- c(.6, 2, 2.1, 4, 4)
sigma <- c(.4, .8, 1.2, 1.6, 2)

# Simulate network
gamma <- rep(1:K, length.out=n) 

seed <- sample(1:10000000, 1)
set.seed(seed)
z <- matrix(NA, nrow=n, ncol=2)
for(k in 1:K) z[gamma == k, ] <- rnorm(sum(gamma ==k) * 2, 0, sigma[k])

y <- matrix(NA, nrow=n, ncol=n)

for(i in 2:n) {
  for(j in 1:(i-1)) {
    if(gamma[i] == gamma[j]) y[i,j] <- y[j,i] <- rbinom(1, 1, expit(beta[gamma[i]] - sqrt(sum( (z[i,] - z[j,])^2 ))))
    else y[i,j] <- y[j,i] <- rbinom(1, 1, B[gamma[i], gamma[j]])
  }
}
#plot_blocked_matrix(y, gamma)

# Hold out edges 
ymiss <- y
inds <- as.matrix(expand.grid(i =1:n, j=1:n))
inds <- inds[inds[,1] < inds[,2],]
inds_gamma <- cbind(gamma[inds[,1]], gamma[inds[,2]])

test <- sample(1:nrow(inds), nrow(inds)/10, replace=FALSE)
train <- setdiff(1:nrow(inds), test)
ymiss[inds[test,]] <- ymiss[inds[test, 2:1]] <-  NA
ytest <- y[inds[test,]]


which_test_within <- which(inds_gamma[test,1] == inds_gamma[test,2])
test_within <- test[which_test_within]
which_test_between <- which(inds_gamma[test,1] != inds_gamma[test,2])
test_between <- test[which_test_between]

ytest_within <- y[inds[test_within,]]
ytest_between <- y[inds[test_between,]]

# Fit LS model
g <- network.initialize(n, directed=FALSE)
net <- network.adjacency(ymiss, g)

control <- ergmm.control(burnin=1000,sample.size= 1000,interval=5)

start <- Sys.time()
ls.fit <- ergmm(net ~ euclidean(d=2), tofit="mcmc", control = control)
ls.time <- as.numeric(Sys.time() - start, units="secs")

ls.pred <- predict(ls.fit)
ls.testpred <- ls.pred[inds[test,2:1]]

ls.auc <- AUC(ls.testpred, ytest)
ls.auc.within <- AUC(ls.pred[inds[test_within, 2:1]], ytest_within)
ls.auc.between <- AUC(ls.pred[inds[test_between, 2:1]], ytest_between)

ls.auprc <- auprc(ls.testpred, ytest)
ls.auprc.within <- auprc(ls.pred[inds[test_within, 2:1]], ytest_within)
ls.auprc.between <- auprc(ls.pred[inds[test_between, 2:1]], ytest_between)

# Fit LPCM - initialize with spectral clustering otherwise it might have trouble
start <- Sys.time()
ls.sbm.ts <- three_stage_est(ymiss, D=2, K=K)
ls.sbm.ts.time <- as.numeric(Sys.time() - start, units='secs')

start <- Sys.time()
lpcm.fit <- ergmm(net ~ euclidean(d=2, G=K), tofit = "mcmc", verbose=TRUE, user.start = list(Z.K=ls.sbm.ts$gamma), control=control)
lpcm.time <- as.numeric(Sys.time() - start, units="secs")

lpcm.pred <- predict(lpcm.fit)
lpcm.pred[upper.tri(lpcm.pred)] <- t(lpcm.pred)[upper.tri(lpcm.pred)] 
lpcm.testpred <- lpcm.pred[inds[test,]]

lpcm.auc <- AUC(lpcm.testpred, ytest)
lpcm.auc.within <- AUC(lpcm.pred[inds[test_within,]], ytest_within)
lpcm.auc.between <- AUC(lpcm.pred[inds[test_between, ]], ytest_between)

lpcm.auprc <- auprc(lpcm.testpred, ytest)
lpcm.auprc.within <- auprc(lpcm.pred[inds[test_within,]], ytest_within)
lpcm.auprc.between <- auprc(lpcm.pred[inds[test_between, ]], ytest_between)

# Fit SBM
sbm.fit <- sbm_mcmc(ymiss, K=K, alpha0=1, print_samples = FALSE, burnin = 1000, samples = 1000, thin=5, a0_within=2, b0_within=2, a0_between = 1, b0_between=3)
sbm.pred.test <- sapply(1:nrow(sbm.fit$memb_samples), function(t) {
  pred <- sbm.fit$B_samples[t,sbm.fit$memb_samples[t,], sbm.fit$memb_samples[t,]]
  pred[inds[test,]]
})
sbm.pred.test <- rowMeans(sbm.pred.test)

sbm.auc <- AUC(sbm.pred.test, ytest)
sbm.auc.within <- AUC(sbm.pred.test[which_test_within], ytest_within)
sbm.auc.between <- AUC(sbm.pred.test[which_test_between], ytest_between)

sbm.auprc <- auprc(sbm.pred.test, ytest)
sbm.auprc.within <- auprc(sbm.pred.test[which_test_within], ytest_within)
sbm.auprc.between <- auprc(sbm.pred.test[which_test_between], ytest_between)

# Fit LS-SBM

dens <- mean(ymiss, na.rm=TRUE)
b0 <- 1
a0 <- (b0*(dens))/(1-dens)

start <- Sys.time()
ls.sbm.fit <- block_latent_MCMC(ymiss, D=2, K=K, burn_in=1000, n_samples=1000, thin=5, rZ=3, v=5, Atheta=matrix(c(1,.4,.4,1),nrow=2), likelihood=FALSE, alpha0 = K, a0=a0, b0=b0, postprocess = FALSE, plot_init=FALSE, debug_output = FALSE)
ls.sbm.time <- as.numeric(Sys.time() - start, units="secs")

predict_our_model <- function(gamma, B, beta, Z) {
  p <- B[gamma, gamma]
  for(k in 1:length(beta)) {
    membs <- which(gamma == k)
    if(length(membs) > 1) p[membs, membs] <- logit.inv(beta[k] - as.matrix(dist(Z[membs,])))
  }
  diag(p) <- 0
  return(p)
}

ls.sbm.preds <- sapply(1:nrow(ls.sbm.fit$Ymiss), function(t) predict_our_model(ls.sbm.fit$gamma[t,], ls.sbm.fit$B[t,,], ls.sbm.fit$beta[t,], ls.sbm.fit$Z[t,,]))

ls.sbm.preds <- rowMeans(ls.sbm.preds) 
ls.sbm.preds <- matrix(ls.sbm.preds, nrow=n, ncol=n)
ls.sbm.testpred <- ls.sbm.preds[inds[test,]]

ls.sbm.auc <- AUC(ls.sbm.testpred, ytest)
ls.sbm.auc.within <- AUC(ls.sbm.testpred[which_test_within], ytest_within)
ls.sbm.auc.between <- AUC(ls.sbm.testpred[which_test_between], ytest_between)

ls.sbm.auprc <- auprc(ls.sbm.testpred, ytest)
ls.sbm.auprc.within <- auprc(ls.sbm.testpred[which_test_within], ytest_within)
ls.sbm.auprc.between <- auprc(ls.sbm.testpred[which_test_between], ytest_between)

# Result 
result <- data.frame(method=rep(c("LS", "LPCM","SBM", "LS-SBM"), each=3),
                     AUC=c(ls.auc, ls.auc.within, ls.auc.between,
                           lpcm.auc,lpcm.auc.within,lpcm.auc.between,
                           sbm.auc, sbm.auc.within, sbm.auc.between,
                           ls.sbm.auc, ls.sbm.auc.within, ls.sbm.auc.between),
                     AUPRC=c(ls.auprc, ls.auprc.within, ls.auprc.between,
                             lpcm.auprc,lpcm.auc.within,lpcm.auprc.between,
                             sbm.auprc, sbm.auprc.within, sbm.auprc.between,
                             ls.sbm.auprc, ls.sbm.auprc.within, ls.sbm.auprc.between),
                     edgetype=rep(c("All", "Within", "Between"), 4),
                     seed=seed,
                     n=n,
                     task_id=task_id)

write.csv(result, file=paste0('results/n', n, 'task', task_id, '.csv'), row.names=FALSE)

