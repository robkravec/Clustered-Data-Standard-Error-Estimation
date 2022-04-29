# User defined functions in different file
source("project_helper_methods.R")
library(mvtnorm)
library(dplyr)
library(KRLS)

sim_coverage <- function(n_clust, obs_per_clust, ICC, nrep = 500) {
  
  # Initialize output data frame
  output_df <- data.frame(N_clusters = rep(n_clust, times = 4),
                          Obs_per_cluser = rep(obs_per_clust, times = 4),
                          ICC = rep(ICC, times = 4),
                          Method = c("Individual",
                                     "Cluster",
                                     "Double",
                                     "Weighted"),
                          Coverage = c(NA, NA, NA, NA),
                          Width = c(NA, NA, NA, NA))
  
  # Simulate data
  dat <- create_clusters(cluster_sizes = rep(x = obs_per_clust, times = n_clust),
                         ICC = ICC, nrep = nrep)
  
  # Individual bootstrap
  output_df[1, c(5, 6)] <- apply(do.call(rbind, 
                                         apply(apply(dat[, c(2:(nrep + 1))], 
                                                     2, boot_ind), 
                                               2, coverage)), 2, mean)
  
  # Cluster bootstrap
  output_df[2, c(5, 6)] <- apply(do.call(rbind, 
                                         apply(apply(dat[, c(2:(nrep + 1))], 
                                                     2, boot_clust, dat$Cluster), 
                                               2, coverage)), 2, mean)
  
  # Double bootstrap
  output_df[3, c(5, 6)] <- apply(do.call(rbind, 
                                         apply(apply(dat[, c(2:(nrep + 1))], 
                                                     2, boot_double, dat$Cluster), 
                                               2, coverage)), 2, mean)
  
  # Double weighted bootstrap
  output_df[4, c(5, 6)] <- apply(do.call(rbind, 
                                         apply(apply(dat[, c(2:(nrep + 1))], 
                                                     2, boot_double, dat$Cluster, use_weights=TRUE), 
                                               2, coverage)), 2, mean)
  
  
  # Return results
  return(output_df)
}

args <- commandArgs(trailingOnly = T)
results <- sim_coverage(n_clust = as.numeric(args[1]), obs_per_clust = as.numeric(args[2]),
                        ICC = as.numeric(args[3]))
saveRDS(results, file = paste("./simulation_data/",args[4], sep='')) 