# 基于 Jaccard 指数的簇匹配算法

# 计算两个簇之间的 Jaccard 指数
jaccard_coefficient <- function(cluster1, cluster2) {
  # 如果两个簇非空，计算 Jaccard 指数
  if (length(cluster1) > 0 && length(cluster2) > 0) {
    return(length(intersect(cluster1, cluster2)) / length(union(cluster1, cluster2)))
  }
  # 如果任意一个簇为空，返回 0
  return(0)
}

# 根据标签提取簇
get_clusters <- function(labels) {
  # 根据标签将样本索引分组
  clusters <- split(seq_along(labels), labels)
  return(list(values = lapply(clusters, as.integer), keys = names(clusters)))
}

# 使用匈牙利算法匹配预测簇和真实簇
match_clusters <- function(predicted_clusters, true_clusters) {
  # 第一步：计算 Jaccard 系数矩阵
  jaccard_matrix <- outer(
    seq_along(predicted_clusters),  # 遍历预测簇
    seq_along(true_clusters),      # 遍历真实簇
    Vectorize(function(i, j) jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]]))
  )
  
  # 第二步：将 Jaccard 系数矩阵转换为代价矩阵 (1 - Jaccard 指数)
  cost_matrix <- 1 - jaccard_matrix
  
  # 第三步：检查矩阵是否为方阵
  nr <- nrow(cost_matrix)  # 行数（预测簇数）
  nc <- ncol(cost_matrix)  # 列数（真实簇数）
  
  if (nr > nc) {
    # 如果预测簇多于真实簇，添加虚拟列
    cost_matrix <- cbind(cost_matrix, matrix(max(cost_matrix) * 2, nrow = nr, ncol = nr - nc))
  } else if (nc > nr) {
    # 如果真实簇多于预测簇，添加虚拟行
    cost_matrix <- rbind(cost_matrix, matrix(max(cost_matrix) * 2, nrow = nc - nr, ncol = nc))
  }
  
  # 第四步：使用匈牙利算法（solve_LSAP）进行最优匹配
  matching <- solve_LSAP(cost_matrix)
  
  # 返回匹配结果（每个预测簇对应的真实簇索引）
  return(as.list(matching))
}

# 重新分配多余的预测簇
reassign_extra_clusters <- function(unmatched_clusters, true_clusters) {
  # 初始化重新分配结果
  reassigned_clusters <- list()
  
  for (cluster in unmatched_clusters) {
    # 找到与当前预测簇最相似的真实簇
    best_match <- which.max(sapply(true_clusters, function(tc) jaccard_coefficient(cluster, tc)))
    # 将多余的预测簇分配到最相似的真实簇中
    reassigned_clusters[[best_match]] <- c(reassigned_clusters[[best_match]], cluster)
  }
  
  return(reassigned_clusters)
}

# 分裂簇（当真实簇比预测簇多时，基于最近距离）
split_cluster <- function(cluster, unmatched_true_cluster, matched_true_clusters, point_coordinates) {
  # cluster: 预测簇中的点索引
  # unmatched_true_cluster: 未匹配真实簇的点索引
  # matched_true_clusters: 已匹配真实簇的点索引列表
  # point_coordinates: 所有点的坐标矩阵，行号为点的索引，列为坐标维度（如 x, y）

  # 初始化子簇和剩余簇
  sub_cluster <- c()  # 将分配到未匹配真实簇的点
  remaining_cluster <- c()  # 将分配到已匹配真实簇的点

  # 获取未匹配真实簇的所有点坐标
  unmatched_coords <- point_coordinates[unmatched_true_cluster, , drop = FALSE]

  # 获取已匹配真实簇的所有点坐标
  matched_coords <- do.call(rbind, lapply(matched_true_clusters, function(cluster) {
    point_coordinates[cluster, , drop = FALSE]
  }))

  # 遍历预测簇中的每个点
  for (sample in cluster) {
    # 获取当前点的坐标
    sample_coords <- point_coordinates[sample, , drop = FALSE]

    # 计算到未匹配真实簇所有点的距离
    distances_to_unmatched <- apply(unmatched_coords, 1, function(point) {
      sqrt(sum((sample_coords - point)^2))
    })

    # 如果已匹配真实簇存在，计算到其所有点的距离
    if (nrow(matched_coords) > 0) {
      distances_to_matched <- apply(matched_coords, 1, function(point) {
        sqrt(sum((sample_coords - point)^2))
      })
    } else {
      distances_to_matched <- Inf  # 如果没有已匹配真实簇，将距离设置为无限大
    }

    # 找到最近距离
    if (min(distances_to_unmatched) < min(distances_to_matched)) {
      # 如果到未匹配真实簇的距离更近，将点分配到未匹配真实簇
      sub_cluster <- c(sub_cluster, sample)
    } else {
      # 否则分配到已匹配真实簇
      remaining_cluster <- c(remaining_cluster, sample)
    }
  }

  return(list(sub_cluster = sub_cluster, remaining_cluster = remaining_cluster))
}

# 主函数：分配簇并匹配
assign_clusters <- function(pred_clusters, true_clusters, pred_labels, true_keys, point_coordinates) {
  # 第一步：初步匹配预测簇和真实簇
  matching <- match_clusters(pred_clusters, true_clusters)
  
  # 第二步：找出未匹配的预测簇和真实簇
  unmatched_pred <- setdiff(seq_along(pred_clusters), matching)
  unmatched_true <- setdiff(seq_along(true_clusters), unlist(matching))
  
  # 第三步：处理未匹配的预测簇
  if (length(unmatched_pred) > 0) {
    # 将多余的预测簇重新分配到真实簇
    extra_assignments <- reassign_extra_clusters(
      lapply(unmatched_pred, function(i) pred_clusters[[i]]),
      true_clusters
    )
  }
  
  # 第四步：处理未匹配的真实簇
  if (length(unmatched_true) > 0) {
    for (u_true in unmatched_true) {
      # 找到与未匹配真实簇最相似的预测簇
      best_pred <- which.max(sapply(pred_clusters, function(pc) jaccard_coefficient(pc, true_clusters[[u_true]])))
      # 对最相似的预测簇进行分裂
      split_result <- split_cluster(pred_clusters[[best_pred]], true_clusters[[u_true]], 
                                    matched_true_clusters = pred_clusters[unlist(matching)], 
                                    point_coordinates = point_coordinates)
      # 更新原簇并添加子簇
      pred_clusters[[best_pred]] <- split_result$remaining_cluster
      pred_clusters <- c(pred_clusters, list(split_result$sub_cluster))
    }
  }
  
  # 返回调整后的预测标签
  # return(reassign_labels(pred_labels, matching, pred_clusters, true_keys))
}

# 匹配函数：将预测标签与真实标签进行匹配
matching_function <- function(true_labels, pred_labels, point_coordinates) {
  # 获取真实簇和预测簇
  true_clusters <- get_clusters(true_labels)
  pred_clusters <- get_clusters(pred_labels)
  
  # 重新分配标签
  adjusted_labels <- assign_clusters(
    pred_clusters$values, 
    true_clusters$values, 
    pred_labels, 
    true_clusters$keys,
    point_coordinates
  )
  
  return(adjusted_labels)
}

# 测试用例
test_case <- function() {
  # 定义真实标签和预测标签
  true_labels <- c(0, 0, 1, 1, 2, 2, 3, 3)  # 4 个真实簇
  pred_labels <- c(0, 0, 1, 1, 2, 2)  # 3 个预测簇
  
  # 所有点的坐标
  point_coordinates <- matrix(c(
    1, 1,  # 点 1
    2, 1,  # 点 2
    5, 5,  # 点 3
    6, 5,  # 点 4
    8, 8,  # 点 5
    9, 8   # 点 6
  ), ncol = 2, byrow = TRUE)
  
  # 执行匹配函数
  result <- matching_function(true_labels, pred_labels, point_coordinates)
  
  print("匹配后的预测标签：")
  print(result)
}

# 运行测试
test_case()
