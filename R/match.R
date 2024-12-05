# Jaccard Index-based clustering matching algorithm

# Jaccard index 
jaccard_coefficient <- function(cluster1, cluster2) {
  if (length(cluster1) > 0 && length(cluster2) > 0) {
    return(length(intersect(cluster1, cluster2)) / length(union(cluster1, cluster2)))
  }
  return(0)
}

get_clusters <- function(labels) {
  clusters <- split(seq_along(labels), labels)
  return(list(values = lapply(clusters, as.integer), keys = names(clusters)))
}

match_clusters <- function(predicted_clusters, true_clusters) {
  # Compute Jaccard coefficient matrix
  jaccard_matrix <- outer(
    seq_along(predicted_clusters),
    seq_along(true_clusters),
    Vectorize(function(i, j) jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]]))
  )
  # To distance matrix(1 - Jaccard coefficient)
  cost_matrix <- 1 - jaccard_matrix
  # Hungarian algorithm to solve
  matching <- solve_LSAP(cost_matrix)
  return(as.list(matching))
}

# Num of predicted clusters > Num of true clusters
reassign_extra_clusters <- function(unmatched_clusters, true_clusters) {
  reassigned_clusters <- list()
  for (cluster in unmatched_clusters) {
    best_match <- which.max(sapply(true_clusters, function(tc) jaccard_coefficient(cluster, tc)))
    reassigned_clusters[[best_match]] <- c(reassigned_clusters[[best_match]], cluster)
  }
  return(reassigned_clusters)
}

# Num of predicted clusters < Num of true clusters
split_cluster <- function(cluster, unmatched_true_cluster) {
  sub_cluster <- list()
  remaining_cluster <- list()
  for (sample in cluster) {
    if (sample %in% unmatched_true_cluster) {
      sub_cluster <- c(sub_cluster, sample)
    } else {
      remaining_cluster <- c(remaining_cluster, sample)
    }
  }
  return(list(sub_cluster = sub_cluster, remaining_cluster = remaining_cluster))
}

assign_clusters <- function(pred_clusters, true_clusters, pred_labels, true_keys) {
  matching <- match_clusters(pred_clusters, true_clusters)
  unmatched_pred <- setdiff(seq_along(pred_clusters), matching)
  unmatched_true <- setdiff(seq_along(true_clusters), unlist(matching))
  
  if (length(unmatched_pred) > 0) {
    extra_assignments <- reassign_extra_clusters(
      lapply(unmatched_pred, function(i) pred_clusters[[i]]),
      true_clusters
    )
  }

  if (length(unmatched_true) > 0) {
    for (u_true in unmatched_true) {
      best_pred <- which.max(sapply(pred_clusters, function(pc) jaccard_coefficient(pc, true_clusters[[u_true]])))
      split_result <- split_cluster(pred_clusters[[best_pred]], true_clusters[[u_true]])
      pred_clusters[[best_pred]] <- split_result$remaining_cluster
      pred_clusters <- c(pred_clusters, list(split_result$sub_cluster))
    }
  }
  return(reassign_labels(pred_labels, matching, pred_clusters, true_keys))
}

matching_function <- function(true_labels, pred_labels) {
  true_clusters <- get_clusters(true_labels)
  pred_clusters <- get_clusters(pred_labels)
  adjusted_labels <- assign_clusters()(pred_clusters$values, true_clusters$values, pred_labels, true_clusters$keys)
  return(adjusted_labels)
}