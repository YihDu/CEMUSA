SGD <- function(
  true_labels , 
  cluster_labels , 
  spatial_coordinates,
  match_cluster_labels = FALSE,
  params = list()
){
  params <- modifyList(list(
      apply_gene_similarity = FALSE,
      apply_anomaly_severity_weight = FALSE,
      gene_exp_matrix = NULL,
      severity_weight_dict = NULL,
      AD_weight = 0.7,
      sigma = 1e-2), params)
  
  # using match function
  if (match_cluster_labels) {
    time_match <- system.time({
      cluster_labels <- matching_function(true_labels , cluster_labels)
    })
    cat('匹配标签时间:', time_match['elapsed'], '秒\n')
    flush.console()
  }
  
  time_build_graph <- system.time({
    graph_list <- process_graph(true_labels, cluster_labels, spatial_coordinates , params)
  })
  cat('构建图时间:', time_build_graph['elapsed'], '秒\n')
  flush.console()

  
  time_analysis <- system.time({
    result <- analyze_graph(graph_list$truth_graph, graph_list$pred_graph , params)
  })
  cat('分析图时间:', time_analysis['elapsed'], '秒\n')
  
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)  
  time_mmd <- system.time({
    mmd_value <- compute_mmd(result$samples_set_truth, result$samples_set_pred, sigma = params$sigma)
  })
  cat('计算MMD时间:', time_mmd['elapsed'], '秒\n')
  
  stopCluster(cl)

  # print
  cat('SGD:', mmd_value, '\n')
}