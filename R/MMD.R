# Compute the Maximum Mean Discrepancy (MMD) between two sets of samples
#########################################

GaussianEMDKernelMultiScale <- function(distances, sigmas) {
  kernel_vals <- exp(- (distances^2) / (2 * sigmas^2))
  mean_kernel <- mean(kernel_vals) # get mean of kernel values
  
  return(mean_kernel)
}

# kernel matrix computation
compute_kernel_matrix_multi_scale <- function(samples1, samples2, sigmas = c(0.1)) {
  
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

# simulate 1 c(0.5)
# simulate 2 c(0.1)
# FP FN      c(0.5)

# hBC        0.2
# MERFISH    0.1
# DLPFC      0.05



compute_mmd_multi_scale_fixed <- function(samples1, samples2, sigmas = c(0.5), nproj = 50) {  
  K_XX <- compute_kernel_matrix_multi_scale(samples1, samples1, sigmas)
  K_YY <- compute_kernel_matrix_multi_scale(samples2, samples2, sigmas)
  K_XY <- compute_kernel_matrix_multi_scale(samples1, samples2, sigmas)
  m <- length(samples1)
  n <- length(samples2)
  mmd <- (sum(K_XX) / (m * m)) + (sum(K_YY) / (n * n)) - (2 * sum(K_XY) / (m * n))
  return(mmd)
}

