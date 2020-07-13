## install the dependencies (required packages)
### CRAN
required.packages <- c("ggplot2","reshape2","gtools","matrixStats","stringr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)

install.packages("ggplot2")
install.packages("reshape2")
install.packages("gtools")
install.packages("matrixStats")
install.packages("stringr")

### biocLite
source("http://bioconductor.org/biocLite.R")
##Install package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")



