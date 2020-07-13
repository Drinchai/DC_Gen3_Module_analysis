## install the dependencies (required packages)
### CRAN
required.packages <- c("ggplot2","reshape2","gtools","matrixstats","stringr")
missing.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(missing.packages)) install.packages(missing.packages)
### biocLite
source("http://bioconductor.org/biocLite.R")
required.biocLite.packages <- c("GEOquery", "ComplexHeatmap")
missing.biocLite.packages <- required.biocLite.packages[!(required.biocLite.packages %in% installed.packages()[,"Package"])]
if(length(missing.biocLite.packages)) biocLite(missing.biocLite.packages)

## load packages into the memory
lapply(required.packages, library, character.only = TRUE)
lapply(required.biocLite.packages, library, character.only = TRUE)


##Install package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")


install.packages("ggplot2")
install.packages("reshape2")
install.packages("gtools")
