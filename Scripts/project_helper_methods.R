
# helper methods for final project, trying to reduce "bloat" in final markdown

# generate-data-function
create_clusters <- function(cluster_sizes, cluster_means = NA, sigma2 = 1,
                            ICC = 0, nrep = 1) {
  
  # Set cluster means to zero if none are provided
  if(is.na(cluster_means)) {
    cluster_means <- rep(0, length(cluster_sizes))
  }
  
  # Calculate covariance for observations in same cluster
  tau2 <- (ICC * sigma2) / (1 - ICC) 
  
  # Calculate cumulative sum of cluster sizes, and add 0 to beginning
  cluster_cumsum <- c(0, cumsum(cluster_sizes))
  
  # Generate block diagonal covariance matrix
  Sigma <- matrix(data = 0, nrow = sum(cluster_sizes),
                  ncol = sum(cluster_sizes))
  for(i in seq_along(cluster_sizes)) {
    temp_sigma <- sigma2 * diag(cluster_sizes[i]) + 
      tau2 * matrix(data = 1, nrow = cluster_sizes[i], ncol = cluster_sizes[i])
    Sigma[(cluster_cumsum[i] + 1): cluster_cumsum[i + 1], 
          (cluster_cumsum[i] + 1): cluster_cumsum[i + 1]] <- temp_sigma
  }
  
  # Draw from multivariate normal distribution
  means <- rep(cluster_means, cluster_sizes)
  draws <- rmvnorm(n = nrep, mean = means, sigma = Sigma)
  
  # Return clustered data
  cluster_counter <- seq(from = 1, to = length(cluster_sizes), by = 1)
  return(data.frame(Cluster = rep(cluster_counter, cluster_sizes),
                    Draws = t(draws)))
}


# individual-bootstrap-function
# Returns n bootstrap means for a vector, x
boot_ind <- function(x, n = 100) {
  
  # Initialize output
  output <- numeric(n)
  
  # Perform bootstrap
  for (i in 1:n) {
    x_boot <- sample(x, size = length(x), replace = T)
    output[i] <- mean(x_boot)
  }
  
  # Return bootstrap means
  return(output)
}


# cluster-bootstrap-function
# Returns n bootstrap means for a vector, x, with clusters, c
boot_clust <- function(x, c, n = 100) {
  
  # Initialize output
  output <- numeric(n)
  
  # Create original dataframe
  dat <- data.frame(Cluster = c, x = x)
  clusters <- unique(c)
  
  # Perform bootstrap
  for (i in 1:n) {
    
    # Determine counts per cluster
    cluster_boot <- sample(clusters, length(clusters), replace = T) %>% 
      as.data.frame()
    colnames(cluster_boot) <- 'Cluster'
    cluster_boot <- cluster_boot %>% 
      count(Cluster)
    
    # Inner join original dataframe and cluster counts
    joined_dat <- dat %>% 
      inner_join(cluster_boot, by = "Cluster")
    
    # Create bootstrapped dataset
    x_boot <- rep(joined_dat$x, joined_dat$n)
    
    # Calculate mean of bootstrapped data
    output[i] <- mean(x_boot)
  }
  
  # Return bootstrap means
  return(output)
}


# double-bootstrap-function
# Returns n bootstrap means for a vector, x, with clusters, c
# Bootstraps by cluster and then within each cluster
boot_double <- function(x, c, n = 100,
                        use_weights = FALSE, sigma = 1) {
  
  # Initialize output
  output <- numeric(n)
  
  # Perform initial data manipulations 
  clusters <- unique(c)
  dat <- data.frame(Cluster = c, x = x) %>% 
    count(Cluster)
  
  # Perform bootstrap
  for (i in 1:n) {
    # Bootstrap clusters
    cluster_boot <- sample(clusters, length(clusters), replace = T)
    cluster_boot_df <- data.frame(Cluster = cluster_boot) %>% 
      inner_join(dat, by = 'Cluster')
    x_boot <- numeric(sum(cluster_boot_df$n))
    cluster_cumsum <- c(0, cumsum(cluster_boot_df$n))
    
    # For each cluster, perform additional bootstrap
    for (j in seq_along(cluster_boot)) {
      temp_dat <- x[c == cluster_boot[j]]
      
      ## if use weights create them...doesn't make sense to do this every time (could be improved)
      if(use_weights){
        weights = assign_sample_weights(temp_dat, sigma = sigma)
        x_boot[(1 + cluster_cumsum[j]) : (cluster_cumsum[j + 1])] = 
          sample(temp_dat, length(temp_dat), replace = T, prob = weights)
      }
      else{
        x_boot[(1 + cluster_cumsum[j]) : (cluster_cumsum[j + 1])] = 
        sample(temp_dat, length(temp_dat), replace = T)
      }
      
    }
    
    # Calculate mean of bootstrapped data
    output[i] <- mean(x_boot)
  }
  
  # Return bootstrap means
  return(output)
}

# coverage-function
coverage <- function(x, true_mean = 0, method = "normal", conf = 0.95) {
  
  # Reject invalid methods
  if(!(method %in% c("quantile", "normal"))) {
    stop("Method must be \"quantile\" or \"normal\"")
  }
  
  # Define alpha
  alpha = 1 - conf
  
  # Use quantiles of bootstrap sample to create confidence interval
  if (method == "quantile") {
    width <- quantile(x, 1 - alpha / 2) - quantile(x, alpha / 2)
    if ((true_mean >= quantile(x, alpha / 2)) & 
        (true_mean <= quantile(x, 1 - alpha / 2))) {
      return(data.frame(Coverage = 1, Width = width))
    } else {
      return(data.frame(Coverage = 0, Width = width))
    }
  }
  
  # Normal approximation
  se <- sd(x)
  samp_mean <- mean(x)
  width <- 2 * qnorm(conf + alpha / 2) * se
  if ((true_mean >= samp_mean + qnorm(alpha / 2) * se) & 
      (true_mean <= samp_mean + qnorm(conf + alpha / 2) * se)) {
    return(data.frame(Coverage = 1, Width = width))
  } else {
    return(data.frame(Coverage = 0, Width = width))
  }
}


# for a given cluster, return sampling weights
## this could maybe be extended to different kernels/higher dimension space
# returns a vector of weights corresponding to each input point
assign_sample_weights = function(X, sigma = 1){
  
  # define kernel...can be moved to parameter of function call
  mydist_matrix = gausskernel(X = as.matrix(X), sigma = sigma)
  
  # use row sum as total distance metric and normalize
  weights = rowSums(mydist_matrix)
  return(weights/sum(weights))
}





