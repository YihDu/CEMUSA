# 定义 Jaccard 系数计算函数
jaccard_coefficient <- function(cluster1, cluster2) {
  if (length(cluster1) > 0 && length(cluster2) > 0) {
    return(length(intersect(cluster1, cluster2)) / length(union(cluster1, cluster2)))
  }
  return(0)
}

# 定义获取簇的函数
get_clusters <- function(labels) {
  clusters <- split(seq_along(labels), labels)
  return(list(values = lapply(clusters, as.integer), keys = names(clusters)))
}

# 定义簇匹配函数
match_clusters <- function(predicted_clusters, true_clusters) {
  jaccard_matrix <- outer(
    seq_along(predicted_clusters),
    seq_along(true_clusters),
    Vectorize(function(i, j) jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]]))
  )
  matching <- solve_LSAP(1 - jaccard_matrix)
  return(as.list(matching))
}

# 定义标签重新分配函数
reassign_labels <- function(pred_labels, matching, pred_clusters, true_keys) {
  adjusted_labels <- pred_labels
  for (i in seq_along(matching)) {
    true_label <- true_keys[[matching[[i]]]]
    adjusted_labels[pred_clusters[[i]]] <- true_label
  }
  return(adjusted_labels)
}

# 定义匹配函数
matching_function <- function(true_labels, pred_labels) {
  true_clusters <- get_clusters(true_labels)
  pred_clusters <- get_clusters(pred_labels)
  
  matching <- match_clusters(pred_clusters$values, true_clusters$values)
  adjusted_labels <- reassign_labels(pred_labels, matching, pred_clusters$values, true_clusters$keys)
  
  return(adjusted_labels)
}
