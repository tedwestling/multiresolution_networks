# ###############################################################################
# # Multi-resolution blockmodel
# #
# # file: india_village_lpcm.R
# # 
# # This file fits the latent position cluster model on village 59 for comparison
# #
# #This file should be placed in the same directory as  the header.R file.
# #
# # Author: tedwestling
# ###############################################################################
rm(list = ls())
##set your working directory to the top level multiresolution_networks folder
setwd("")
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

library(VBLPCM)
library(latentnet)

#igraph and sna don't play well together sometimes
detach("package:igraph") 
v.start<-vblpcmstart(network(network,directed=F),G=K)
v.fit<-vblpcmfit(v.start)
##
pdf(file="plots/vlpcm.pdf",height=8,width=8)
plot(v.fit, main="Latent Position Cluster Model",xlab="First latent dimension",ylab="Second latent dimension",cex.axis=1.5,cex.lab=1.3,pty=1:6)
text(x=-5,y=6,labels=expression(paste("Posterior mode of ",beta[0],"=.003",sep='')),cex=1.3)
dev.off()