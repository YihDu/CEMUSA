get_edge_attributes <- function(graph, 
                                unique_groups,
                                params) {
  group_to_onehot <- setNames(lapply(unique_groups, function(group) {
      as.numeric(unique_groups == group)
    }), unique_groups)
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
    
    if (params$apply_gene_similarity) {
      gene_similarity_weight <- edge_attr(graph, "gene_similarity_weight", edge)

      if (is.null(gene_similarity_weight) || length(gene_similarity_weight) == 0) {
      } else {
      }
      encoding <- encoding * ifelse(is.na(gene_similarity_weight), 1.0, gene_similarity_weight)
    }
    
    if (params$apply_anomaly_severity_weight) {
      cat("Applying anomaly severity weight when building graphs.\n")
      anomaly_severity_weight <- igraph::edge_attr(graph, "anomaly_severity_weight", edge)
      cat("anomaly_severity_weight:", anomaly_severity_weight , "\n")
      encoding <- encoding * ifelse(is.na(anomaly_severity_weight), 1.0, anomaly_severity_weight)
      cat("encoding:", encoding , "\n")
    }
    samples <- append(samples, list(as.numeric(encoding)))
  }
  # To a matrix: (number_of_edge) * (encoding_dim)
  samples <- do.call(rbind, samples)
  return(samples)
}

fit_kde_and_sample <- function(samples, num_samples, sample_times, bandwidth = NULL, random_seed = NULL) {
  sklearn <- import("sklearn.neighbors")
  np <- import("numpy")
  samples_py <- r_to_py(as.matrix(samples))
  dim <- ncol(samples)
  kde <- sklearn$KernelDensity(kernel = "gaussian", bandwidth = bandwidth)
  kde$fit(samples_py)
  samples_list <- vector("list", sample_times)
  for (i in seq_len(sample_times)) {
    sampled <- kde$sample(n_samples = as.integer(num_samples), random_state = as.integer(random_seed + i))
    sampled <- np$clip(sampled, 0, 1)
    samples_list[[i]] <- py_to_r(sampled)
  }
  return(samples_list)
}

get_unique_groups <- function(truth_graph, pred_graph) {
  groups_in_truth_G <- unique(igraph::V(truth_graph)$label)
  groups_in_pred_G <- unique(igraph::V(pred_graph)$label)
  unique_groups <- sort(unique(c(groups_in_truth_G, groups_in_pred_G)))
  return(unique_groups)
}

analyze_graph <- function(truth_graph, pred_graph , params) {
  sample_times <- 4
  unique_groups <- get_unique_groups(truth_graph, pred_graph)
  samples_truth <- get_edge_attributes(truth_graph, unique_groups , params)
  samples_pred <- get_edge_attributes(pred_graph, unique_groups , params)
  num_samples <- nrow(samples_truth)
  samples_set_truth <- fit_kde_and_sample(samples_truth, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  samples_set_pred <- fit_kde_and_sample(samples_pred, num_samples, sample_times, bandwidth = 0.1, random_seed = 42)
  return(list(samples_set_truth = samples_set_truth, samples_set_pred = samples_set_pred))
}


