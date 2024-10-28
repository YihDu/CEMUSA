# build KNN graph and add coordinates/label
# igraph package?
library(data.table)

build_graph3 <- function(label_data, coordinate_data, num_neighbors = 6) {

  
  # 计算K近邻
  nbrs <- nn2(coordinate_data, k = num_neighbors + 1)
  
  # 使用data.table构建边列表
  dt <- data.table(
    source = rep(seq_len(nrow(nbrs$nn.idx)), each = num_neighbors),
    target = as.vector(t(nbrs$nn.idx[, -1]))
  )
  
  # 创建无向图并去重边
  graph <- graph_from_data_frame(dt, directed = FALSE) %>% simplify()
  
  # 分配顶点属性
  V(graph)$label <- label_data
  V(graph)$pos <- split(coordinate_data, seq(nrow(coordinate_data)))
  
  return(graph)
}

build_graph1 <- function(label_data, coordinate_data, num_neighbors = 6) {
  graph <- igraph::make_empty_graph(n = nrow(coordinate_data), directed = FALSE)
  for (i in seq_len(nrow(coordinate_data))) {
    igraph::V(graph)$pos[i] <- list(as.numeric(coordinate_data[i, ]))
  }
  igraph::V(graph)$label <- label_data  # vector
  nbrs <- RANN::nn2(coordinate_data, k = num_neighbors + 1)
  for (i in seq_len(nrow(nbrs$nn.idx))) {
    for (n in nbrs$nn.idx[i, -1]) {
      graph <- igraph::add_edges(graph, c(i, n))
    }
  }
  return(graph)
}

build_graph2 <- function(label_data, coordinate_data, num_neighbors = 6) {
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
  # 获取truth_graph的所有边属性
  truth_gene_weights <- edge_attr(truth_graph, "gene_similarity_weight", E(truth_graph))
  truth_anomaly_weights <- edge_attr(truth_graph, "anomaly_severity_weight", E(truth_graph))
  
  # 直接复制边属性到pred_graph
  if ("gene_similarity_weight" %in% edge_attr_names(truth_graph)) {
    set_edge_attr(pred_graph, "gene_similarity_weight", value = truth_gene_weights)
  }
  
  if ("anomaly_severity_weight" %in% edge_attr_names(truth_graph)) {
    set_edge_attr(pred_graph, "anomaly_severity_weight", value = truth_anomaly_weights)
  }
  
  return(pred_graph)
}


process_graph <- function(true_labels, cluster_labels , coordinate_data , params) {

  truth_graph <- build_graph2(true_labels, coordinate_data)
  pred_graph <- build_graph2(cluster_labels, coordinate_data)
  
  if (params$apply_anomaly_severity_weight) {
    message("Applying anomaly severity weight when building graphs.")
    calculate_anomaly_weight(truth_graph, params$severity_weight_dict)
  }

  if (params$apply_gene_similarity) {
    message("Applying gene similarity weight when building graphs.") 
    calculate_gene_similarity(truth_graph, params$gene_exp_matrix)
  }

  # 如果 上述有一个为真，则复制权重
  if (params$apply_anomaly_severity_weight || params$apply_gene_similarity) {
    copy_weights(truth_graph , pred_graph)
  }

  return(list(truth_graph = truth_graph, pred_graph = pred_graph))
  }


