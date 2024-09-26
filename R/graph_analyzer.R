library(MASS)
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
  
  result <- do.call(rbind, samples)
  return(result)
}

fit_kde_and_sample <- function(samples, num_samples, sample_times, bandwidth = NULL, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)
  
  kde <- kde2d(samples[, 1], samples[, 2], n = num_samples, h = bandwidth)
  
  z_flat <- as.vector(kde$z)
  pdf <- z_flat / sum(z_flat)  

  cdf <- cumsum(pdf)
  
  samples_set <- list()
  
  x_coords <- rep(kde$x, each = length(kde$y))
  y_coords <- rep(kde$y, times = length(kde$x))
  
  for (i in seq_len(sample_times)) {
    sampled_indices <- findInterval(runif(num_samples), cdf)
    
    sampled_x <- x_coords[sampled_indices]
    sampled_y <- y_coords[sampled_indices]
    
    sampled <- cbind(sampled_x, sampled_y)
    
    sampled <- pmax(pmin(sampled, 1), 0)
    
    samples_set <- append(samples_set, list(sampled))
  }
  print("samples_set:")
  print(samples_set)
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
  sample_times <- 1
  
  unique_groups <- get_unique_groups(truth_graph, pred_graph)
  
  samples_truth <- get_edge_attributes(truth_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  samples_pred <- get_edge_attributes(pred_graph, apply_gene_similarity, apply_anomaly_severity_weight, apply_distance_weight, unique_groups)
  
  cat("samples_truth:", samples_truth , '\n')
  cat("samples_pred:", samples_pred, '\n')

  cat("samples_truth :", typeof(samples_truth), "  ", dim(samples_truth), '\n')
  cat("samples_pred :", typeof(samples_pred), "  ", dim(samples_pred), '\n')

  cat("samples_truth[,1]:", samples_truth[,1], '\n')
  cat("samples_turth[,2]:", samples_pred[,2], '\n')
  cat("samples_pred[,1]:", samples_pred[,1], '\n')
  cat("samples_pred[,2]:", samples_pred[,2], '\n')

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