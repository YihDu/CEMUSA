## weight calculation module

# 1.anomaly weight
calculate_anomaly_weight <- function(graph, dict_severity_levels) {
  severity_mapping <- setNames(sapply(dict_severity_levels, function(x) x$severity_level), 
                               sapply(dict_severity_levels, function(x) x$name))
  if (ecount(graph) == 0) {
    stop("Graph contains no edges.")
  }
  
  for (i in seq_along(E(graph))) {
    edge <- E(graph)[i]
    
    u <- ends(graph, edge)[1]
    v <- ends(graph, edge)[2]
    
    group_u <- V(graph)$label[u]
    group_v <- V(graph)$label[v]
    
    if (!group_u %in% names(severity_mapping)) {
      stop(paste("Group", group_u, "not found in severity mapping"))
    }
    if (!group_v %in% names(severity_mapping)) {
      stop(paste("Group", group_v, "not found in severity mapping"))
    }
    
    severity_u <- severity_mapping[group_u]
    severity_v <- severity_mapping[group_v]
    
    anomaly_severity_weight <- (severity_u + severity_v) / 2
    print(paste('Anomaly Severity Weight:', anomaly_severity_weight))
    
    set_edge_attr(graph, "anomaly_severity_weight", index = i, value = anomaly_severity_weight)
  }
}

# 2.gene similarity weight

# For same label assign Similarity to Ground Truth edge and copy
# sim  
# different
# 1-sim

# pearson similarity
calculate_pearson_similarity <- function(x, y) {
  vx <- x - mean(x)
  vy <- y - mean(y)
  return(sum(vx * vy) / (sqrt(sum(vx ^ 2)) * sqrt(sum(vy ^ 2))))
}

# gene similarity
calculate_gene_similarity <- function(graph, anndata, is_preprocessed = TRUE) {
  if (!is_preprocessed) {
    preprocessed_data <- preprocess_anndata(anndata)
    gene_expression_matrix <- preprocessed_data[, preprocessed_data@var$highly_variable]
  } else {
    gene_expression_matrix <- anndata@X
  }
  
  gene_expression_matrix <- as.matrix(gene_expression_matrix)
  
  group_means <- list()
  group_indices <- list()
  
  for (i in seq_len(vcount(graph))) {
    group <- V(graph)$group[i]
    if (!group %in% names(group_indices)) {
      group_indices[[group]] <- c()
    }
    group_indices[[group]] <- c(group_indices[[group]], i)
  }
  
  for (group in names(group_indices)) {
    group_expression_matrix <- gene_expression_matrix[group_indices[[group]], ]
    group_mean_vector <- colMeans(group_expression_matrix)
    group_means[[group]] <- group_mean_vector
  }
  
  pearson_matrix <- cor(gene_expression_matrix)
  
  for (edge in E(graph)) {
    endpoints <- ends(graph, edge)
    u <- endpoints[1]
    v <- endpoints[2]
    
    group_u <- V(graph)$group[u]
    group_v <- V(graph)$group[v]

    if (group_u == group_v) {
      group_mean <- group_means[[group_u]]
      similarity_u <- calculate_pearson_similarity(gene_expression_matrix[u, ], group_mean)
      similarity_v <- calculate_pearson_similarity(gene_expression_matrix[v, ], group_mean)
      E(graph)[edge]$gene_similarity_weight <- 0.5 * (similarity_u + similarity_v)
    } 
    else {
      E(graph)[edge]$gene_similarity_weight <- 1 - pearson_matrix[u, v]
    }
    
    print(paste('Gene Similarity Weight:', E(graph)[edge]$gene_similarity_weight))
  }
}