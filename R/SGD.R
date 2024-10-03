SGD <- function(
  true_labels , 
  cluster_labels , 
  spatial_coordinates,
  params = list()
){
  params <- modifyList(list(
      apply_gene_similarity = FALSE,
      apply_anomaly_severity_weight = FALSE,
      AD_weight = 0.7,
      sigma = 1e-2), params)
 graph_list <- process_graph(true_labels, cluster_labels, spatial_coordinates , params)
 result <- analyze_graph(graph_list$truth_graph, graph_list$pred_graph)
 cat("result$sample_set_truth type and shape:", typeof(result$samples_set_truth), "  ", length(result$samples_set_truth), '\n')
 cat("result$sample_set_pred type and shape:", typeof(result$samples_set_pred), "  ", length(result$samples_set_pred), '\n')
 mmd_value <- compute_mmd(result$samples_set_truth, result$samples_set_pred, sigma = 1)
 return(mmd_value)
}