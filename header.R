###############################################################################
# Multi-resolution blockmodel
#
# file: header.R
# This is the header file for the multi-resolution blockmodel project. It contains
# project meta-data and functions to be included in every program file in the project.
# Author: tedwestling
###############################################################################

# Source function files
source('code/functions/latent_space_vb.R')
source('code/functions/spectral_clustering.R')
source('code/functions/three_stage_estimation.R')
source('code/functions/mcmc_sampler.R')

###############################################################################
# Function: ensure_package
# --------------------------
# This function makes it checks whether a package is installed. If so, it loads it
# and if not it installs it, then loads it.
###############################################################################

ensure_package <- function(package_name) {
  for (i in 1:length(package_name)) {
    name <- package_name[i]
    if(!suppressMessages(require(name, character.only=TRUE, quietly=TRUE, warn.conflicts = FALSE))) {
      install.packages(name, repos='http://cran.cnr.Berkeley.edu')
      suppressMessages(library(name, character.only=TRUE, quietly=TRUE, warn.conflicts = FALSE))
    }
  }
}

ensure_package(c('igraph', 'foreign', 'ggplot2', 'plyr', 'MASS', 'mgcv', 'reshape', 'boot', 'grid', 'sna', 'gridExtra', 'mixer', 'VBLPCM', 'rARPACK', 'EMCluster', 'gtools', 'MCMCpack', 'ellipse', 'clue', 'network', 'latentnet', 'smacof', 'abind', 'emdbook', 'coda'))

logit <- function(x) log(x/(1-x))
expit <- inv.logit <- logit.inv <- function(x) 1/ (1 + exp(-x))