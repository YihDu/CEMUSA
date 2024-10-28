build_graph <- function(label_data, coordinate_data, num_neighbors = 6) {
  nbrs <- nn2(coordinate_data, k = num_neighbors + 1)  
  sources <- rep(seq_len(nrow(nbrs$nn.idx)), each = num_neighbors)
  targets <- as.vector(t(nbrs$nn.idx[, -1]))  # 去掉自身的索引并展平
  edge_list <- c(rbind(sources, targets))

  graph <- make_graph(edge_list, directed = FALSE) %>% simplify()

  V(graph)$label <- label_data
  V(graph)$pos <- split(coordinate_data, seq(nrow(coordinate_data)))

  return(graph)
}


copy_weights <- function(truth_graph, pred_graph) {
  truth_gene_weights <- edge_attr(truth_graph, "gene_similarity_weight", E(truth_graph))
  truth_anomaly_weights <- edge_attr(truth_graph, "anomaly_severity_weight", E(truth_graph))
  
  if ("gene_similarity_weight" %in% edge_attr_names(truth_graph)) {
    set_edge_attr(pred_graph, "gene_similarity_weight", value = truth_gene_weights)
  }
  
  if ("anomaly_severity_weight" %in% edge_attr_names(truth_graph)) {
    set_edge_attr(pred_graph, "anomaly_severity_weight", value = truth_anomaly_weights)
  }
  
  return(pred_graph)
}


process_graph <- function(true_labels, cluster_labels , coordinate_data , params) {

  truth_graph <- build_graph(true_labels, coordinate_data)
  pred_graph <- build_graph(cluster_labels, coordinate_data)
  
  if (params$apply_anomaly_severity_weight) {
    message("Applying anomaly severity weight when building graphs.")
    calculate_anomaly_weight(truth_graph, params$severity_weight_dict)
  }

  if (params$apply_gene_similarity) {
    message("Applying gene similarity weight when building graphs.") 
    calculate_gene_similarity(truth_graph, params$gene_exp_matrix)
  }

  if (params$apply_anomaly_severity_weight || params$apply_gene_similarity) {
    copy_weights(truth_graph , pred_graph)
  }

  return(list(truth_graph = truth_graph, pred_graph = pred_graph))
  }


