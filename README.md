# SGD: **S**patial **G**rouping **D**iscrepancy

`SGD` is a evaluation metric for Spatial Transcriptomics data. It addresses the limitations of existing clustering evaluation metrics by accounting for label agreement, spatial
locations, and error severity simultaneously.

## Installation
### 1. Install the package from GitHub
```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("YihDu/SGD")
```

### 2. Load the package
```r
library(SGD)
```

## How to use `SGD`
Please refer to the [documentation](https://yihdu.github.io/SGD/) for details.
