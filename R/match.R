library(clue)

jaccard_coefficient <- function(cluster1, cluster2) {
  intersection <- length(intersect(cluster1, cluster2)) 
  union <- length(unique(c(cluster1, cluster2)))  
  return(intersection / union)  
}

get_clusters <- function(labels) {
  clusters <- split(seq_along(labels), labels)
  return(list(values = lapply(clusters, as.integer), keys = names(clusters)))
}


match_clusters <- function(predicted_clusters, true_clusters) {
  jaccard_matrix <- outer(
    seq_along(predicted_clusters), 
    seq_along(true_clusters),      
    Vectorize(function(i, j) jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]]))
  )
  
  cat("Jaccard 系数矩阵：\n")
  print(jaccard_matrix)

  cost_matrix <- 1 - jaccard_matrix
  
  matching <- solve_LSAP(cost_matrix)

  return(matching)
}

reassign_labels <- function(pred_labels, matching, pred_clusters, true_keys) {
  label_mapping <- setNames(true_keys[matching], names(pred_clusters))
  
  reassigned_labels <- sapply(pred_labels, function(label) {
    matching_index <- which(names(label_mapping) == label)
    if (length(matching_index) > 0) {
      return(label_mapping[matching_index])
    } else {
      return(label) 
    }
  })
  
  return(reassigned_labels)
}


assign_clusters <- function(pred_clusters, true_clusters, pred_labels, true_keys, point_coordinates) {
  matching <- match_clusters(pred_clusters, true_clusters)
  
  cat("匹配结果：\n")
  print(matching)

  return(reassign_labels(pred_labels, matching, pred_clusters, true_keys))
}

matching_function <- function(true_labels, pred_labels) {
  true_clusters <- get_clusters(true_labels)
  pred_clusters <- get_clusters(pred_labels)
  
  adjusted_labels <- assign_clusters(
    pred_clusters$values, 
    true_clusters$values, 
    pred_labels, 
    true_clusters$keys
  )
  
  return(adjusted_labels)
}