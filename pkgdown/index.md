---
title: "Documentation"
description: "User Guide of SAS package"
template: home
---

## Introduction

`SAS` is a evaluation metric for Spatial Transcriptomics data. It addresses the limitations of existing clustering evaluation metrics by accounting for label agreement, spatial
locations, and error severity simultaneously.

## Installation

To install the package from GitHub, use the following command:

```r
# Install the package from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("YihDu/SAS")
```

## Spatial Alignment Score (SAS)

```r
library(SAS)
SAS(true_labels, cluster_labels, spatial_coordinates , match_cluster_labels = TRUE , params = list)
```

## Parameters
- **true_labels** :  
  Ground truth class labels to be used as a reference.

- **cluster_labels** : 
  Cluster labels to evaluate.

- **spatial_coordinates** : 
  A matrix or data frame containing the spatial coordinates (e.g., `x`, `y`) of each data point.

- **match_cluster_labels** : (default = TRUE)  
  If `TRUE`, the function will attempt to match the   `cluster_labels`  with the `true_labels` using an internal matching function. This is useful when the labels are not already matched (e.g., matched by external information like marker genes).

- **params** : list (optional)  
  A list of additional parameters that control more detailed aspects of the evaluation.

## Example Usage

### Example 1: Basic Usage (Reproduce the Case I in the paper)

```r
load('Simulate_Case_CenterEdge.RData')

metadata <- seurat_object@meta.data
coordinates <- metadata[, c("spatial_x", "spatial_y")]
truth_labels <- metadata$truth_label
pred_labels_edge <- metadata$edge_error

SAS = SAS(true_labels = truth_labels, 
           cluster_labels = pred1_labels , 
           spatial_coordinates = coordinates , 
           match_cluster_labels = FALSE)
```r

Click [here]() to download the data metioned above.

### Example 2: When considering the Error Severity

```r


```r


## Cite `SAS`
Jiaying Hu<sup>†</sup>, Yihang Du<sup>†</sup>, Suyang Hou, Yueyang Ding, Hao Wu and Xiaobo Sun&#35;.*SAS:A clustering evaluation metric for spatial transcriptomics.*,2024







