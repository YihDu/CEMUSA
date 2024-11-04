## weight calculation module

# 1.anomaly weight
calculate_anomaly_weight <- function(graph, dict_severity_levels) {
  severity_mapping <- setNames(sapply(dict_severity_levels, function(x) x$severity_level), 
                               sapply(dict_severity_levels, function(x) x$name))
  if (ecount(graph) == 0) {
    stop("Graph contains no edges.")
  }
  
  for (i in seq_along(E(graph))) {
    edge <- E(graph)[i]
    
    u <- ends(graph, edge)[1]
    v <- ends(graph, edge)[2]
    
    group_u <- V(graph)$label[u]
    group_v <- V(graph)$label[v]
    
    if (!group_u %in% names(severity_mapping)) {
      stop(paste("Group", group_u, "not found in severity mapping"))
    }
    if (!group_v %in% names(severity_mapping)) {
      stop(paste("Group", group_v, "not found in severity mapping"))
    }
    
    severity_u <- severity_mapping[group_u]
    severity_v <- severity_mapping[group_v]
    
    anomaly_severity_weight <- (severity_u + severity_v) / 2
    set_edge_attr(graph, "anomaly_severity_weight", index = i, value = anomaly_severity_weight)
  }
}

# 2.gene similarity weight

# For same label assign Similarity to Ground Truth edge and copy
# sim  
# different
# 1-sim

# pearson similarity
# 计算基因相似度
# 计算皮尔逊相似度
# 计算皮尔逊相似度的函数
calculate_pearson_similarity <- function(x, y) {
  cor(x, y, use = "complete.obs")
}

# 计算基因相似度
calculate_gene_similarity <- function(graph, gene_expression_matrix) {
  cat("开始计算Gene similarity\n")
  flush.console()  

  # 使用 split 优化索引构建，并确保 V(graph)$label 中的每个元素都有对应的索引
  time_split <- system.time({
    group_indices <- split(seq_len(nrow(gene_expression_matrix)), V(graph)$label)
    group_means <- lapply(group_indices, function(indices) {
      colMeans(gene_expression_matrix[indices, ], na.rm = TRUE)
    })
  })
  cat('计算组均值时间:', time_split['elapsed'], '秒\n')
  flush.console()

  # 获取所有边的端点
  time_ends <- system.time({
    endpoints <- ends(graph, E(graph), names = FALSE)
  })
  cat('获取边端点时间:', time_ends['elapsed'], '秒\n')
  flush.console()

  # 遍历图中的每条边
  time_loop <- system.time({
    for (edge in igraph::E(graph)) {
      u <- igraph::tail_of(graph, edge)
      v <- igraph::head_of(graph, edge)
      group_u <- igraph::V(graph)[u]$label
      group_v <- igraph::V(graph)[v]$label

      if (group_u == group_v) {
        group_mean <- group_means[[group_u]]
        similarity_u <- calculate_pearson_similarity(gene_expression_matrix[u, ], group_mean)
        similarity_v <- calculate_pearson_similarity(gene_expression_matrix[v, ], group_mean)
        # graph <- set_edge_attr(graph, "gene_similarity_weight", edge, 0.5 * (similarity_u + similarity_v))
        # graph <- set_edge_attr(graph, "gene_similarity_weight", edge, 1)
        graph <- set_edge_attr(graph, "gene_similarity_weight", edge, 1 - calculate_pearson_similarity(gene_expression_matrix[u, ], gene_expression_matrix[v, ]))
      } else {
        graph <- set_edge_attr(graph, "gene_similarity_weight", edge, 1 - calculate_pearson_similarity(gene_expression_matrix[u, ], gene_expression_matrix[v, ]))
      }

      # 打印出来数值
      cat("这条边的gene Similarity weight:", edge_attr(graph, "gene_similarity_weight", edge), "\n")
    }
  })
  cat('遍历边并计算相似度时间:', time_loop['elapsed'], '秒\n')
  flush.console()

  return (graph)
}