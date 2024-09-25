library(MASS)
# 提取边属性并应用权重
get_edge_attributes <- function(graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups) {

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

    if (group_u == group_v) {
      encoding <- group_to_onehot[[group_u]]
      
    } else {
      encoding <- rep(0, length(unique_groups))
    }
    cat("encoding:", encoding, "\n")
    cat("encoding:", as.numeric(encoding), "\n")
    
    # 应用基因相似性权重
    if (apply_gene_similarity) {
      gene_similarity_weight <- igraph::edge_attr(graph, "gene_similarity_weight", edge)
      encoding <- encoding * ifelse(is.na(gene_similarity_weight), 1.0, gene_similarity_weight)
    }
    
    # 应用异常严重度权重
    if (apply_anomaly_severity_weight) {
      anomaly_severity_weight <- igraph::edge_attr(graph, "anomaly_severity_weight", edge)
      encoding <- encoding * ifelse(is.na(anomaly_severity_weight), 1.0, anomaly_severity_weight)
    }
    
    # 应用距离权重
    if (apply_distance_weight) {
      distance_weight <- igraph::edge_attr(graph, "distance_weight", edge)
      encoding <- encoding * ifelse(is.na(distance_weight), 1.0, distance_weight)
    }
    
    samples <- append(samples, list(as.numeric(encoding)))
  }
  
  # 转换为矩阵
  result <- do.call(rbind, samples)
  return(result)
}

# 进行核密度估计并采样
fit_kde_and_sample <- function(samples, num_samples, sample_times, bandwidth = NULL, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)
  
  # 生成二维核密度估计
  kde <- kde2d(samples[, 1], samples[, 2], n = num_samples, h = bandwidth)
  
  # 展开 kde$z 并将其归一化为概率密度函数 (PDF)
  z_flat <- as.vector(kde$z)
  pdf <- z_flat / sum(z_flat)  # 归一化

  # 生成累计分布函数 (CDF)
  cdf <- cumsum(pdf)
  
  # 准备用于存储采样结果的列表
  samples_set <- list()
  
  # 展开网格点的坐标
  x_coords <- rep(kde$x, each = length(kde$y))
  y_coords <- rep(kde$y, times = length(kde$x))
  
  for (i in seq_len(sample_times)) {
    # 从 CDF 中进行采样，生成随机数并找到对应的网格索引
    sampled_indices <- findInterval(runif(num_samples), cdf)
    
    # 根据索引提取相应的 x 和 y 坐标
    sampled_x <- x_coords[sampled_indices]
    sampled_y <- y_coords[sampled_indices]
    
    # 合并为二维坐标对
    sampled <- cbind(sampled_x, sampled_y)
    
    # 将采样值限制在 [0, 1] 范围内，确保采样点合法
    sampled <- pmax(pmin(sampled, 1), 0)
    
    # 保存每次的采样结果
    samples_set <- append(samples_set, list(sampled))
  }
  print("samples_set:")
  print(samples_set)
  return(samples_set)
}

# 获取图中唯一的分组信息
get_unique_groups <- function(truth_graph, pred_graph) {
  groups_in_truth_G <- unique(igraph::V(truth_graph)$label)
  groups_in_pred_G <- unique(igraph::V(pred_graph)$label)
  unique_groups <- sort(unique(c(groups_in_truth_G, groups_in_pred_G)))
  return(unique_groups)
}

# 主函数，分析图并生成样本集
analyze_graph <- function(truth_graph, pred_graph) {
  apply_gene_similarity <- FALSE
  apply_anomaly_severity_weight <- FALSE
  apply_distance_weight <- FALSE
  sample_times <- 1
  
  # 获取唯一分组
  unique_groups <- get_unique_groups(truth_graph, pred_graph)
  
  # 提取真实图和预测图的边属性样本
  samples_truth <- get_edge_attributes(truth_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  samples_pred <- get_edge_attributes(pred_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  
  cat("samples_truth:", samples_truth , '\n')
  cat("samples_pred:", samples_pred, '\n')

  #打印数据类型和形状
  cat("samples_truth 的类型和形状:", typeof(samples_truth), "  ", dim(samples_truth), '\n')
  cat("samples_pred 的类型和形状:", typeof(samples_pred), "  ", dim(samples_pred), '\n')

  cat("samples_truth[,1]:", samples_truth[,1], '\n')
  cat("samples_turth[,2]:", samples_pred[,2], '\n')
  cat("samples_pred[,1]:", samples_pred[,1], '\n')
  cat("samples_pred[,2]:", samples_pred[,2], '\n')

  # 计算核密度估计并采样
  # 样本数量就是samples_truth的行数
  num_samples <- nrow(samples_truth)
  cat('num_samples:', num_samples)
  samples_set_truth <- fit_kde_and_sample(samples_truth, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  samples_set_pred <- fit_kde_and_sample(samples_pred, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  
  print("samples_set_truth_length:")
  print(length(samples_set_truth))
  print("samples_set_pred_length:")
  print(length(samples_set_pred))

  return(list(samples_set_truth = samples_set_truth, samples_set_pred = samples_set_pred))
}