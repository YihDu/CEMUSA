get_edge_attributes <- function(graph, unique_groups, params) {
  node_labels <- igraph::V(graph)$label
  edge_list <- igraph::as_edgelist(graph)
  gene_similarity_weights <- igraph::edge_attr(graph, "gene_similarity_weight")
  anomaly_severity_weights <- igraph::edge_attr(graph, "anomaly_severity_weight")
  
  group_to_onehot <- setNames(lapply(unique_groups, function(group) {
    as.numeric(unique_groups == group)
  }), unique_groups)
  
  num_edges <- nrow(edge_list)
  encoding_dim <- length(unique_groups)
  samples <- matrix(0, nrow = num_edges, ncol = encoding_dim)



  for (i in seq_len(num_edges)) {
    u <- edge_list[i, 1]
    v <- edge_list[i, 2]
    group_u <- node_labels[u]
    group_v <- node_labels[v]

    if (group_u == group_v) {
      encoding <- group_to_onehot[[group_u]]
    } else {
      encoding <- rep(0, encoding_dim)
    }
    
    if (params$apply_gene_similarity) {
      gene_similarity_weight <- gene_similarity_weights[i]
      encoding <- encoding * ifelse(is.na(gene_similarity_weight), 1.0, gene_similarity_weight)
    }
    
    if (params$apply_anomaly_severity_weight) {
      # cat("Applying anomaly severity weight when building graphs.\n")
      anomaly_severity_weight <- anomaly_severity_weights[i]
      # cat("anomaly_severity_weight:", anomaly_severity_weight, "\n")
      encoding <- encoding * ifelse(is.na(anomaly_severity_weight), 1.0, anomaly_severity_weight)
      # cat("encoding:", encoding, "\n")
    }
    
    samples[i, ] <- encoding
  }
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
    # sampled <- np$clip(sampled, 0, 1)
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



