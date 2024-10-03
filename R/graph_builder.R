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
  for (i in seq_along(E(truth_graph))) {
    u <- as.numeric(ends(truth_graph, E(truth_graph)[i])[1])
    v <- as.numeric(ends(truth_graph, E(truth_graph)[i])[2])
    
    if (are_adjacent(pred_graph, u, v)) {
      edge_id <- get.edge.ids(pred_graph, c(u, v))
      
      if ("gene_similarity_weight" %in% edge_attr_names(truth_graph)) {
        truth_gene_weight <- edge_attr(truth_graph, "gene_similarity_weight", i)
        if (!is.na(truth_gene_weight)) {
          set_edge_attr(pred_graph, "gene_similarity_weight", edge_id, truth_gene_weight)
        }
      }
      
      if ("anomaly_severity_weight" %in% edge_attr_names(truth_graph)) {
        truth_anomaly_weight <- edge_attr(truth_graph, "anomaly_severity_weight", i)
        if (!is.na(truth_anomaly_weight)) {
          set_edge_attr(pred_graph, "anomaly_severity_weight", edge_id, truth_anomaly_weight)
        }
      }
      
      if ("distance_weight" %in% edge_attr_names(truth_graph)) {
        truth_distance_weight <- edge_attr(truth_graph, "distance_weight", i)
        if (!is.na(truth_distance_weight)) {
          set_edge_attr(pred_graph, "distance_weight", edge_id, truth_distance_weight)
        }
      }
    }
  }
}



process_graph <- function(true_labels, cluster_labels , spatial_coordinates , params) {

  truth_graph <- build_graph(true_labels, spatial_coordinates, num_neighbors = 6)
  pred_graph <- build_graph(cluster_labels, spatial_coordinates, num_neighbors = 6)
  
  if (params$apply_anomaly_severity_weight) {
    message("Applying anomaly severity weight when building graphs.")
    calculate_anomaly_weight(truth_graph, params$severity_levels)
  }

  if (params$apply_gene_similarity) {
    message("Applying gene similarity weight when building graphs.") 
    calculate_gene_similarity(truth_graph, params$gene_expression_matrix)
  }
  copy_weights(truth_graph , pred_graph)
  return(list(truth_graph = truth_graph, pred_graph = pred_graph))
  }


