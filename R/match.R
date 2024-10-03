jaccard_coefficient <- function(cluster1, cluster2) {
  if (length(cluster1) > 0 & length(cluster2) > 0) {
    intersection_counts <- length(intersect(cluster1, cluster2))
    union_counts <- length(union(cluster1, cluster2))
    return(intersection_counts / union_counts)
  } else {
    return(0)
  }
}

get_clusters_from_labels <- function(labels) {
  clusters <- split(1:length(labels), labels)
  cluster_values <- lapply(clusters, function(x) as.integer(x))
  cluster_keys <- names(clusters)
  return(list(cluster_values, cluster_keys))
}

reassign_clusters <- function(predicted_clusters, true_clusters, matching) {
  reassigned_clusters <- vector("list", length(true_clusters))
  for (pred_index in names(matching)) {
    true_index <- matching[[pred_index]]
    reassigned_clusters[[as.integer(true_index)]] <- union(reassigned_clusters[[as.integer(true_index)]], predicted_clusters[[as.integer(pred_index)]])
  }
  return(reassigned_clusters)
}

match_clusters <- function(predicted_clusters, true_clusters) {
  jaccard_matrix <- matrix(0, nrow = length(predicted_clusters), ncol = length(true_clusters))
  for (i in 1:length(predicted_clusters)) {
    for (j in 1:length(true_clusters)) {
      jaccard_matrix[i, j] <- jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]])
    }
  }
  cost_matrix <- 1 - jaccard_matrix
  matching <- solve_LSAP(cost_matrix)
  return(as.list(matching))
}

reassign_labels <- function(pred_labels, matching, pred_clusters, true_labels_keys) {
  adjusted_labels <- pred_labels
  for (pred_index in names(matching)) {
    true_index <- matching[[pred_index]]
    true_label <- true_labels_keys[[as.integer(true_index)]]
    for (node in pred_clusters[[as.integer(pred_index)]]) {
      adjusted_labels[node] <- true_label
    }
  }
  return(adjusted_labels)
}

label_matching <- function(dataframe, pred_col_name) {
  true_labels <- dataframe$labels
  pred_labels <- dataframe[[pred_col_name]]
  
  clusters_info_true <- get_clusters_from_labels(true_labels)
  true_clusters <- clusters_info_true[[1]]
  true_labels_keys <- clusters_info_true[[2]]
  
  clusters_info_pred <- get_clusters_from_labels(pred_labels)
  pred_clusters <- clusters_info_pred[[1]]
  
  matching <- match_clusters(pred_clusters, true_clusters)
  adjusted_labels <- reassign_labels(pred_labels, matching, pred_clusters, true_labels_keys)
  
  dataframe$matched_label <- adjusted_labels
  return(dataframe[, c("matched_label")])
}
