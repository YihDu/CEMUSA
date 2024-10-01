calculate_gene_similarity <- function(graph, gene_expression_matrix) {
  for (edge in igraph::E(graph)) {
    u <- igraph::tail_of(graph, edge)$name
    v <- igraph::head_of(graph, edge)$name
    similarity <- cor(gene_expression_matrix[u, ], gene_expression_matrix[v, ])
    igraph::E(graph)[u %--% v]$gene_similarity_weight <- similarity
  }
}

calculate_anomaly_weight <- function(graph, severity_levels) {
  for (edge in igraph::E(graph)) {
    u <- igraph::tail_of(graph, edge)$name
    v <- igraph::head_of(graph, edge)$name
    severity <- mean(severity_levels[c(u, v)])
    igraph::E(graph)[u %--% v]$anomaly_severity_weight <- severity
  }
}

# match function
