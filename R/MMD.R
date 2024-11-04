
compute_kernel_matrix <- function(samples1, samples2, sigma = 1.0) {
  n <- dim(samples1)[1]
  m <- dim(samples2)[1]
  
  K_values <- foreach(i = 1:n, .combine = 'rbind', .export = c("sliced_wasserstein_distance", "GaussianEMDKernel", "rand_projections")) %dopar% {
    K_row <- numeric(m)
    for (j in 1:m) {
      x_i <- samples1[i, , ]
      y_j <- samples2[j, , ]
      swd_value <- sliced_wasserstein_distance(x_i, y_j)
      K_row[j] <- GaussianEMDKernel(swd_value, sigma)
    }
    K_row
  }
  
  K <- t(K_values)
  return(K)
}


# 串行的版本
sliced_wasserstein_distance <- function(first_samples, second_samples, num_projections = 16, num_repeat_projections = 32, p = 2, seed = 42) {
  dim <- ncol(first_samples)
  wasserstein_distances <- numeric(num_repeat_projections)
  for (i in 1:num_repeat_projections) {
    projections <- rand_projections(dim, num_projections, seed + i)
    first_projections <- first_samples %*% t(projections)
    second_projections <- second_samples %*% t(projections)
    first_sorted <- apply(first_projections, 2, sort)
    second_sorted <- apply(second_projections, 2, sort)
    
    wasserstein_distance <- (colMeans(abs(first_sorted - second_sorted)^p))^(1/p)
    wasserstein_distances[i] <- mean(wasserstein_distance)
  }
  return(mean(wasserstein_distances))
}


if(FALSE){
# 并行计算
sliced_wasserstein_distance <- function(first_samples, second_samples, num_projections = 16, num_repeat_projections = 32, p = 2, seed = 42) {
  dim <- ncol(first_samples)
  
  # 定义一个内部函数来计算单次重复的 Wasserstein 距离
  compute_wasserstein_distance <- function(i) {
    projections <- rand_projections(dim, num_projections, seed + i)
    first_projections <- first_samples %*% t(projections)
    second_projections <- second_samples %*% t(projections)
    first_sorted <- apply(first_projections, 2, sort)
    second_sorted <- apply(second_projections, 2, sort)
    
    wasserstein_distance <- (colMeans(abs(first_sorted - second_sorted)^p))^(1/p)
    return(mean(wasserstein_distance))
  }
  
  # 注册并行后端
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  # 在并行环境中导入所需的函数和包
  clusterExport(cl, c("rand_projections", "sliced_wasserstein_distance"))
  
  # 使用 foreach 并行计算 Wasserstein 距离
  wasserstein_distances <- foreach(i = 1:num_repeat_projections, .combine = c, .packages = "parallel") %dopar% {
    compute_wasserstein_distance(i)
  }
  
  # 停止并行后端
  stopCluster(cl)
  
  return(mean(wasserstein_distances))
}
}


compute_mmd <- function(samples1, samples2, sigma = 1.0) {
  cat("使用的sigma是: ", sigma, '\n')
  K_XX <- compute_kernel_matrix(samples1, samples1, sigma)
  K_YY <- compute_kernel_matrix(samples2, samples2, sigma)
  K_XY <- compute_kernel_matrix(samples1, samples2, sigma)
  mmd <- mean(K_XX) + mean(K_YY) - 2 * mean(K_XY)
  return(mmd)
}








rand_projections <- function(dim, num_projections = 16, seed = 42) {
  set.seed(seed)
  projections <- matrix(rnorm(num_projections * dim), nrow = num_projections)
  projections <- projections / sqrt(rowSums(projections^2))
  return(projections)
}

# Gaussian + Sliced Wasserstein Distance kernel
GaussianEMDKernel <- function(swd_value, sigma = 1.0) {
  return(exp(-swd_value^2 / (2 * sigma^2)))
}