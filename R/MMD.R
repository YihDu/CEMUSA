# Compute the Maximum Mean Discrepancy (MMD) between two sets of samples
#########################################

GaussianEMDKernelMultiScale <- function(distances, sigmas) {
  kernel_vals <- exp(- (distances^2) / (2 * sigmas^2))
  mean_kernel <- mean(kernel_vals) # get mean of kernel values
  
  return(mean_kernel)
}

# kernel matrix computation
compute_kernel_matrix_multi_scale <- function(samples1, samples2, sigmas = c(0.01)) {
  
  sample_times <- length(samples1)
  K_values <- foreach(i = 1:sample_times, .combine = 'rbind', .export = c("GaussianEMDKernelMultiScale", "swdist")) %dopar% {
    K_row <- numeric(sample_times)
    for (j in 1:sample_times) {
      x_i <- samples1[[i]]
      y_j <- samples2[[j]]
      set.seed(42)
      swd_value <- swdist(x_i, y_j, nproj = 50)
      K_row[j] <- GaussianEMDKernelMultiScale(swd_value$distance, sigmas)
    }
    K_row
  }
  K <- t(K_values)
  return(K)
}

# MMD Computation
compute_mmd_multi_scale_fixed <- function(samples1, samples2, sigmas = c(0.01 , 0.1 , 0.5), nproj = 50) {  
  K_XX <- compute_kernel_matrix_multi_scale(samples1, samples1, sigmas)
  K_YY <- compute_kernel_matrix_multi_scale(samples2, samples2, sigmas)
  K_XY <- compute_kernel_matrix_multi_scale(samples1, samples2, sigmas)

  cat('K_XX:')
  print(K_XX)
  cat('----------------')
  cat('K_YY:')
  print(K_YY)
  cat('----------------')
  cat('K_XY:')
  print(K_XY)
  
  m <- length(samples1)
  n <- length(samples2)
  mmd <- (sum(K_XX) / (m * m)) + (sum(K_YY) / (n * n)) - (2 * sum(K_XY) / (m * n))
  return(mmd)
}