get_edge_attributes <- function(graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups) {
  # group to one-hot
#   example:group_to_onehot <- list(
#   A = c(1, 0, 0),
#   B = c(0, 1, 0),
#   C = c(0, 0, 1)
# )

  group_to_onehot <- setNames(lapply(unique_groups, function(group) {
      as.numeric(unique_groups == group)
    }), unique_groups)
  
  # print("group_to_onehot:", group_to_onehot)
  cat("unique_groups:", unique_groups , "\n")

  samples <- list()
  
  for (edge in igraph::E(graph)) {
    u <- igraph::tail_of(graph, edge)
    v <- igraph::head_of(graph, edge)
    group_u <- igraph::V(graph)[u]$label
    group_v <- igraph::V(graph)[v]$label

    # same group non-zero encoding
    if (group_u == group_v) {
      encoding <- group_to_onehot[[group_u]]
      
    } else {
      encoding <- rep(0, length(unique_groups))
    }
    
    if (apply_gene_similarity) {
      gene_similarity_weight <- igraph::edge_attr(graph, "gene_similarity_weight", edge)
      encoding <- encoding * ifelse(is.na(gene_similarity_weight), 1.0, gene_similarity_weight)
    }
    
    if (apply_anomaly_severity_weight) {
      anomaly_severity_weight <- igraph::edge_attr(graph, "anomaly_severity_weight", edge)
      encoding <- encoding * ifelse(is.na(anomaly_severity_weight), 1.0, anomaly_severity_weight)
    }
    
    if (apply_distance_weight) {
      distance_weight <- igraph::edge_attr(graph, "distance_weight", edge)
      encoding <- encoding * ifelse(is.na(distance_weight), 1.0, distance_weight)
    }
    
    samples <- append(samples, list(as.numeric(encoding)))
  }
  
  # to a matrix: (number of edge) * (encoding dim)
  samples <- do.call(rbind, samples)
  return(samples)
}

# fit_kde_and_sample <- function(samples, num_samples, sample_times, bandwidth = NULL, random_seed = NULL) {
#   if (!is.null(random_seed)) set.seed(random_seed)
  
#   # 拟合 KDE
#   cat('开始拟合KDE...\n')
#   flush.console()
#   # 计算时间
#   time_kde <- system.time({
#     kde <- kde(x = samples, h = bandwidth , binned = TRUE)
#   })
#   cat('拟合KDE时间:', time_kde['elapsed'], '秒\n')
#   flush.console()
  
#   # 获取样本的维度
#   dim <- ncol(samples)
  
#   # 初始化多维数组
#   samples_set <- array(NA, dim = c(sample_times, num_samples, dim))
  
#   # 采样
#   for (i in seq_len(sample_times)) {
#     sampled_points <- rkde(n = num_samples, fhat = kde)
#     samples_set[i, , ] <- sampled_points
#   }
  
#   return(samples_set)
# }

fit_kde_and_sample <- function(samples, num_samples, sample_times, bandwidth = NULL, random_seed = NULL) {
  sklearn <- import("sklearn.neighbors")
  np <- import("numpy")
  
  # 将 R 的数据转换为 Python 格式
  samples_py <- r_to_py(as.matrix(samples))
  dim <- ncol(samples)
  # 创建并拟合 KDE 模型
  kde <- sklearn$KernelDensity(kernel = "gaussian", bandwidth = bandwidth)
  kde$fit(samples_py)
  
  # 初始化样本集合
  samples_set <- array(NA, dim = c(sample_times, num_samples, dim))
  
  # 采样
  for (i in seq_len(sample_times)) {
    # 使用 kde$sample() 生成样本
    sampled <- kde$sample(n_samples = as.integer(num_samples), random_state = as.integer(random_seed + i))
    
    # 将样本裁剪在 [0, 1] 范围内
    sampled <- np$clip(sampled, 0, 1)
    
    # 将采样结果转换回 R 格式，并保存
    samples_set[i, , ] <- py_to_r(sampled)
  }
  
  return(samples_set)
}

get_unique_groups <- function(truth_graph, pred_graph) {
  groups_in_truth_G <- unique(igraph::V(truth_graph)$label)
  groups_in_pred_G <- unique(igraph::V(pred_graph)$label)
  unique_groups <- sort(unique(c(groups_in_truth_G, groups_in_pred_G)))
  return(unique_groups)
}

analyze_graph <- function(truth_graph, pred_graph) {
  apply_gene_similarity <- FALSE
  apply_anomaly_severity_weight <- FALSE
  apply_distance_weight <- FALSE
  
  # sample times
  sample_times <- 4
  
  unique_groups <- get_unique_groups(truth_graph, pred_graph)
  
  samples_truth <- get_edge_attributes(truth_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  samples_pred <- get_edge_attributes(pred_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  num_samples <- nrow(samples_truth)


  # 打印形状
  cat("samples_truth shape:", dim(samples_truth), '\n')
  cat("samples_pred shape:", dim(samples_pred), '\n')

  samples_set_truth <- fit_kde_and_sample(samples_truth, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  samples_set_pred <- fit_kde_and_sample(samples_pred, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  
  # 打印形状
  cat("samples_set_truth shape:", dim(samples_set_truth), '\n')
  cat("samples_set_pred shape:", dim(samples_set_pred), '\n')
  return(list(samples_set_truth = samples_set_truth, samples_set_pred = samples_set_pred))
}