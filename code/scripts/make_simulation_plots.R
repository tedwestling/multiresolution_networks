###############################################################################
# Multi-resolution blockmodel
#
# file: make_simulation_plots.R

# This file takes the simulation results from the cluster and makes plots
# 
# Author: tedwestling
###############################################################################

# The result files were concatenated on the cluster as follows:
# files <- list.files('results', full.names=TRUE)
# library(plyr)
# results <- ldply(files, read.csv)

load('data/simulation_results/results.Rdata')

results$method <- factor(results$method, levels=c("LS", "LPCM", "SBM", "LS-SBM"))

library(plyr)
comparisons <- ddply(results, .(task_id, n, edgetype), function(df) {
  with(df, data.frame(value=c(AUC[method == "LS-SBM"] / AUC[method == "LS"],
                              AUC[method == "LS-SBM"] / AUC[method == "LPCM"],
                              AUC[method == "LS-SBM"] / AUC[method == "SBM"],
                              AUPRC[method == "LS-SBM"] / AUPRC[method == "LS"],
                              AUPRC[method == "LS-SBM"] / AUPRC[method == "LPCM"],
                              AUPRC[method == "LS-SBM"] / AUPRC[method == "SBM"]),
                      comparison=rep(c("LS-SBM / LS", "LS-SBM / LPCM", "LS-SBM / SBM"),2),
                      metric=rep(c("AUC", "AUPRC"), each=3)))
})

comparisons$comparison <- factor(comparisons$comparison, levels=c("LS-SBM / LS", "LS-SBM / LPCM", "LS-SBM / SBM"))

(g <- ggplot(subset(comparisons, metric == "AUC" & n == 300)) +
    geom_boxplot(aes(comparison, value)) +
    facet_grid(~edgetype) + 
    theme_minimal()+
    geom_hline(yintercept=1, linetype=2) +
    xlab(NULL) + 
    ylab("Relative AUC"))
ggsave('plots/simulation_relative_AUC.png', g, width=6.5, height=3, units='in' ,scale=1.5)


(g <- ggplot(subset(comparisons, metric == "AUPRC" & n == 300)) +
    geom_boxplot(aes(comparison, value)) +
    facet_grid(~edgetype) + 
    theme_minimal()+
    geom_hline(yintercept=1, linetype=2) +
    xlab(NULL) + 
    ylab("Relative AUPRC"))
ggsave('plots/simulation_relative_AUPRC.png', g, width=6.5, height=3, units='in' ,scale=1.5)

(g <- ggplot(subset(results, n= 300)) +
    geom_boxplot(aes(method, AUPRC)) +
    facet_grid(~edgetype) +
    theme_minimal() +
    xlab(NULL)) 