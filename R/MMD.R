# gaussian_kernel <- function(x, y, sigma) {
#   swd_value <- swd_batch(x, y)
#   return(exp(- (swd_value^2) / (2 * sigma^2)))
# }

# rand_projections <- function(dim, num_projections = 16, random_seed = NULL) {
#   if (!is.null(random_seed)) {
#     set.seed(random_seed)
#   }
#   projections <- matrix(rnorm(num_projections * dim), ncol = dim)
#   projections <- projections / sqrt(rowSums(projections^2))
#   return(projections)
# }

# swd_batch <- function(vectors1, vectors2, n_repeat_projection = 512, proj_per_repeat = 16, random_seed = 42) {
#   if (!all(dim(vectors1) == dim(vectors2))) {
#     stop("The dimensions of vectors1 and vectors2 must match.")
#   }
  
#   batch_size <- dim(vectors1)[1]
#   vector_dim <- dim(vectors1)[2]
  
#   distances <- numeric(n_repeat_projection)
  
#   for (i in seq_len(n_repeat_projection)) {
#     rand <- rand_projections(vector_dim, proj_per_repeat, random_seed = random_seed + i)
    
#     proj1 <- apply(vectors1, 1, function(vec) vec %*% t(rand))
#     proj2 <- apply(vectors2, 1, function(vec) vec %*% t(rand))
    
#     proj1 <- t(apply(proj1, 1, sort))
#     proj2 <- t(apply(proj2, 1, sort))
    
#     d <- rowMeans((proj1 - proj2)^2)
#     distances[i] <- mean(d)
#   }
  
#   result <- mean(distances)
#   return(result)
# }

# compute_kernel_matrix <- function(samples1, samples2, kernel_function, sigma) {
#   n <- length(samples1)
#   m <- length(samples2)

#   cat("m: ", m, "  n: ", n, '\n')

#   K <- matrix(0, n, m)
  
#   for (i in seq_len(n)) {
#     for (j in seq_len(m)) {
#       K[i, j] <- kernel_function(samples1[[i]], samples2[[j]], sigma)  
#     }
#   }
  
#   return(K)
# }

# compute_mmd <- function(samples1, samples2, sigma = 1) {
#   kernel_function <- function(x, y, sigma) {
#     gaussian_kernel(x, y, sigma)
#   }
#   K_XX <- compute_kernel_matrix(samples1, samples1, kernel_function, sigma)
#   K_YY <- compute_kernel_matrix(samples2, samples2, kernel_function, sigma)
#   K_XY <- compute_kernel_matrix(samples1, samples2, kernel_function, sigma)
  
#   return(mean(K_XX) + mean(K_YY) - 2 * mean(K_XY))
# }


# 计算 MMD 的函数
compute_mmd <- function(samples1, samples2, sigma = 1.0) {
  cat("使用的sigma是: ", sigma, '\n')
  K_XX <- compute_kernel_matrix(samples1, samples1, sigma)
  K_YY <- compute_kernel_matrix(samples2, samples2, sigma)
  K_XY <- compute_kernel_matrix(samples1, samples2, sigma)
  mmd <- mean(K_XX) + mean(K_YY) - 2 * mean(K_XY)
  return(mmd)
}

compute_kernel_matrix <- function(samples1, samples2, sigma = 1.0) {
  n <- dim(samples1)[1]
  m <- dim(samples2)[1]
  
  # 并行计算核矩阵元素
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
  
  # 组装核矩阵
  K <- t(K_values)
  return(K)
}

sliced_wasserstein_distance <- function(first_samples, second_samples, num_projections = 512, p = 2) {
  dim <- ncol(first_samples)
  
  # 生成随机投影矩阵
  projections <- rand_projections(dim, num_projections)  # (num_projections, dim)
  
  # 计算投影
  first_projections <- first_samples %*% t(projections)  # (m, num_projections)
  second_projections <- second_samples %*% t(projections)  # (m, num_projections)
  
  # 对每个投影维度排序
  first_sorted <- apply(first_projections, 2, sort)
  second_sorted <- apply(second_projections, 2, sort)
  
  # 计算 Wasserstein 距离
  wasserstein_distance <- (colMeans(abs(first_sorted - second_sorted)^p))^(1/p)
  
  # 返回平均 SWD 值
  return(mean(wasserstein_distance))
}

# 生成随机投影矩阵
rand_projections <- function(dim, num_projections = 512 , seed = 42) {
  
  set.seed(seed)
  
  # 生成形状为 (num_projections, dim) 的随机矩阵
  projections <- matrix(rnorm(num_projections * dim), nrow = num_projections)
  
  # 归一化每个投影向量，使其单位长度
  projections <- projections / sqrt(rowSums(projections^2))
  return(projections)
}

# 高斯 EMD 核函数
GaussianEMDKernel <- function(swd_value, sigma = 1.0) {
  return(exp(-swd_value^2 / (2 * sigma^2)))
}