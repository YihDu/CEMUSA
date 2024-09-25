build_graph <- function(label_data, coordinate_data, num_neighbors = 6) {
  graph <- igraph::make_empty_graph(n = nrow(coordinate_data), directed = FALSE)
  for (i in seq_len(nrow(coordinate_data))) {
    igraph::V(graph)$pos[i] <- list(as.numeric(coordinate_data[i, ]))
  }
  igraph::V(graph)$label <- label_data
  nbrs <- RANN::nn2(coordinate_data, k = num_neighbors + 1)
  for (i in seq_len(nrow(nbrs$nn.idx))) {
    for (n in nbrs$nn.idx[i, -1]) {
      graph <- igraph::add_edges(graph, c(i, n))
    }
  }
  return(graph)
}

copy_weights <- function(truth_graph, pred_graph) {
  for (edge in E(truth_graph)) {
    u <- tail_of(truth_graph, edge)
    v <- head_of(truth_graph, edge)
    
    if (are_adjacent(pred_graph, u, v)) {
      edge_id <- get.edge.ids(pred_graph, c(u, v))
      
      if (!is.null(edge_attr(truth_graph, "gene_similarity_weight", edge))) {
        edge_attr(pred_graph, "gene_similarity_weight", edge_id) <- edge_attr(truth_graph, "gene_similarity_weight", edge)
      }
      
      if (!is.null(edge_attr(truth_graph, "anomaly_severity_weight", edge))) {
        edge_attr(pred_graph, "anomaly_severity_weight", edge_id) <- edge_attr(truth_graph, "anomaly_severity_weight", edge)
      }
      
      if (!is.null(edge_attr(truth_graph, "distance_weight", edge))) {
        edge_attr(pred_graph, "distance_weight", edge_id) <- edge_attr(truth_graph, "distance_weight", edge)
      }
    }
  }
}

process_graph <- function(true_labels, cluster_labels , spatial_coordinates) {

  truth_graph <- build_graph(true_labels, spatial_coordinates, num_neighbors = 6)
  pred_graph <- build_graph(cluster_labels, spatial_coordinates, num_neighbors = 6)
  
  if (FALSE) {

  
  if (config$apply_gene_similarity) {
    message("Applying gene similarity weight to truth graph.")
    
    calculate_gene_similarity(truth_graph, gene_expression_matrix)
  }
  
  if (config$apply_anomaly_severity_weight) {
    message("Applying anomaly severity weight to truth graph.")
    calculate_anomaly_weight(truth_graph, config$severity_levels)
  }

  }
  
  return(list(truth_graph = truth_graph, pred_graph = pred_graph))
}


