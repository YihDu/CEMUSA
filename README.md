# SAS: **S**patial **A**lignment **S**core

`SAS` is a evaluation metric for Spatial Transcriptomics data. It addresses the limitations of existing clustering evaluation metrics by accounting for label agreement, spatial
locations, and error severity simultaneously.

## Installation
### 1. Install the package from GitHub
```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("YihDu/SAS")
```

### 2. Load the package
```r
library(SAS)
```

## How to use `SAS`
Please refer to the [documentation](https://yihdu.github.io/SGD/) for details.
