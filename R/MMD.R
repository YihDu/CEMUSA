# 高斯核函数
gaussian_kernel <- function(x, y, sigma) {
  swd_value <- swd_batch(x, y)
  return(exp(- (swd_value^2) / (2 * sigma^2)))
}

# 生成随机投影矩阵
rand_projections <- function(dim, num_projections = 16, random_seed = NULL) {
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  projections <- matrix(rnorm(num_projections * dim), ncol = dim)
  projections <- projections / sqrt(rowSums(projections^2))
  return(projections)
}

# 计算 Sliced Wasserstein Distance (SWD)
swd_batch <- function(vectors1, vectors2, n_repeat_projection = 512, proj_per_repeat = 16, random_seed = 42) {
  if (!all(dim(vectors1) == dim(vectors2))) {
    stop("The dimensions of vectors1 and vectors2 must match.")
  }
  
  batch_size <- dim(vectors1)[1]
  vector_dim <- dim(vectors1)[2]
  
  distances <- numeric(n_repeat_projection)
  
  for (i in seq_len(n_repeat_projection)) {
    rand <- rand_projections(vector_dim, proj_per_repeat, random_seed = random_seed + i)
    
    # 使用 apply() 处理三维数组中的每个批次
    proj1 <- apply(vectors1, 1, function(vec) vec %*% t(rand))
    proj2 <- apply(vectors2, 1, function(vec) vec %*% t(rand))
    
    proj1 <- t(apply(proj1, 1, sort))
    proj2 <- t(apply(proj2, 1, sort))
    
    d <- rowMeans((proj1 - proj2)^2)
    distances[i] <- mean(d)
  }
  
  result <- mean(distances)
  return(result)
}

# 计算核矩阵
compute_kernel_matrix <- function(samples1, samples2, kernel_function, sigma) {
  n <- length(samples1)
  m <- length(samples2)

  cat("m: ", m, "  n: ", n, '\n')

  K <- matrix(0, n, m)
  
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      K[i, j] <- kernel_function(samples1[[i]], samples2[[j]], sigma)  # 使用每对样本计算核函数
    }
  }
  
  return(K)
}

# 计算 MMD
compute_mmd <- function(samples1, samples2, sigma = 1, is_hist = FALSE) {
  if (is_hist) {
    samples1 <- lapply(samples1, function(x) x / sum(x))
    samples2 <- lapply(samples2, function(x) x / sum(x))
  }
  
  # 使用高斯核计算 MMD
  kernel_function <- function(x, y, sigma) {
    gaussian_kernel(x, y, sigma)
  }
  
  K_XX <- compute_kernel_matrix(samples1, samples1, kernel_function, sigma)
  K_YY <- compute_kernel_matrix(samples2, samples2, kernel_function, sigma)
  K_XY <- compute_kernel_matrix(samples1, samples2, kernel_function, sigma)
  
  return(mean(K_XX) + mean(K_YY) - 2 * mean(K_XY))
}
