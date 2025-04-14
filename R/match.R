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
  n_pred <- length(predicted_clusters)
  n_true <- length(true_clusters)
  
  cost_matrix_orig <- 1 - outer(
    seq_along(predicted_clusters), 
    seq_along(true_clusters),
    Vectorize(function(i, j) jaccard_coefficient(predicted_clusters[[i]], true_clusters[[j]]))
  )
  
  if(n_pred == n_true) {
    cost_matrix <- cost_matrix_orig
  } else if(n_pred > n_true) {
    padding <- matrix(1, nrow = n_pred, ncol = n_pred - n_true)
    cost_matrix <- cbind(cost_matrix_orig, padding)
  } else { 
    padding <- matrix(1, nrow = n_true - n_pred, ncol = n_true)
    cost_matrix <- rbind(cost_matrix_orig, padding)
  }
  
  matching_all <- solve_LSAP(cost_matrix)

  if(n_pred >= n_true) {
    assigned <- matching_all[1:n_pred]
    assigned[assigned > n_true] <- NA  
  } else {  
    assigned <- matching_all[1:n_pred]
  }
  
  return(assigned)
}

reassign_labels <- function(pred_labels, matching, pred_clusters, true_keys) {
  
  label_mapping <- sapply(names(pred_clusters), function(key, matching, true_keys) {
    idx <- which(names(pred_clusters) == key)
    mapped <- matching[idx]
    if (is.na(mapped)) {
      return(NA)
    } else {
      return(true_keys[mapped])
    }
  }, matching = matching, true_keys = true_keys, USE.NAMES = TRUE)
  
  reassigned_labels <- sapply(pred_labels, function(label) {
    if(!is.na(label_mapping[label])) {
      return(label_mapping[label])
    } else {
      return(label)
    }
  })
  
  return(reassigned_labels)
}

assign_clusters <- function(pred_clusters, true_clusters, pred_labels, true_keys, point_coordinates = NULL) {
  matching <- match_clusters(pred_clusters, true_clusters)
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
