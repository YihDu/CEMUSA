############################################
##                CEMUSA                  ##
############################################

CEMUSA <- function(
  true_labels , 
  cluster_labels , 
  spatial_coordinates,
  match_cluster_labels = FALSE,
  params = list()
){

  if (is.null(params) || any(is.na(params))) {
    params <- list()
  }
  if (!is.list(params)) {
    stop("params must be a list")
  }
  default_params <- list(
    apply_gene_similarity = FALSE,
    apply_anomaly_severity_weight = FALSE,
    gene_exp_matrix = NULL,
    severity_weight_dict = NULL
  )
  params <- modifyList(default_params, params)
  
  if (match_cluster_labels) {
    cluster_labels <- matching_function(true_labels, cluster_labels)
  }

  graph_list <- process_graph(true_labels, cluster_labels, spatial_coordinates , params)
  result <- analyze_graph(graph_list$truth_graph, graph_list$pred_graph , params) 
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)  
  mmd_value <- compute_mmd_multi_scale_fixed(result$samples_set_truth, result$samples_set_pred)
  stopCluster(cl)
  cat("CEMUSA Value: ", mmd_value, "\n")
  return(mmd_value)
}