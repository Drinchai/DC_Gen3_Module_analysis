# Setup environment
rm(list=ls())

# Set working directory (location on your computer)
setwd("~/DC_Gen3_Module_analysis")

# load dependencies
library(matrixStats)
library(gtools)
library(stringr)
library(ComplexHeatmap)

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106"

##load module matrix
load("./R data/DC_ModuleGen3_2019.Rdata")

##Load module matrix
load(paste0("./R data/", GSE_ID, "_", platform, "_mod_matrix.Rdata"))

#prepare data for analysis
###########
df_raw = as.matrix(dat.mod.Gen3)          # replace "dat.mod.Gen3" with data_matrix in raw expression data
sample_info = sample.info                 # replace "sample_info" with sample annotation 
mod_func = dat.mod.func.Gen3              # repleace "mod_func" with Gene module annotation table  

#### make sure that expression matrix and sample information are the same order
df_raw = df_raw[,rownames(sample_info)]
colnames(df_raw) == rownames(sample_info)

# Difference
Diff.mod.ind.sin <- df_raw[,]
Diff.mod.ind.sin [,] <- NA

k=1
for (k in 1:nrow(df_raw)) {
  signature = rownames(df_raw)[k]
  test.table <- sample_info 
  test.table$scores <- df_raw[k,]
  T4 <- test.table
  T3 <- test.table[test.table$Group %in% c("Control"),]
  Diff.mod.ind.sin[k,] <- (T4$scores-(mean(T3$scores)))
}

Diff.mod.ind.sin <- as.data.frame(Diff.mod.ind.sin)

## fold change 
FC.mod.ind.sin <- df_raw[,]
FC.mod.ind.sin [,] <- NA

for (k in 1:nrow(df_raw)) {
  signature = rownames(df_raw)[k]
  test.table <- sample_info 
  test.table$scores <- df_raw[k,]
  T4 <- test.table
  T3 <- test.table[test.table$Group %in% c("Control"),]
  FC.mod.ind.sin[k,] <- foldchange(T4$scores,mean(T3$scores))
}

FC.mod.ind.sin <- as.data.frame(FC.mod.ind.sin)


#############################################
# Calculate percentage of response ##
############################################ 
#logical check ##
test.up <- (FC.mod.ind.sin > 1.5) + (Diff.mod.ind.sin > 100) == 2          # TRUE Up gene, Both TRUE
test.down <- (FC.mod.ind.sin < -1.5) + (Diff.mod.ind.sin < -100) == 2      # TRUE down gene, Both TRUE

### prepare gene annotation table
Gene.matrix = dat.mod.func.Gen3[rownames(test.up),] 
#####UP GENE#######
up.mods.sin <- data.frame(Module = row.names(test.up), test.up,gene=0)                    # create a new blank table
up.mods.sin [,] <- NA
up.mods.sin <- up.mods.sin [-c(2:nrow(up.mods.sin)),]

for (i in 1:length(unique(Gene.matrix$Module))){                                         # length of module
  module <- unique(as.character(Gene.matrix$Module))[i]                                  # look for only unique module
  sums <- colSums(test.up[Gene.matrix$Module==module,])                                  # sum upgene of each column by module 
  genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])                                # sum number of gene in each module
  up.mods.sin <- rbind(up.mods.sin,c(module,sums,genes))                                 # paste result into a new fake table
}

up.mods.sin <- up.mods.sin [-1,]

rownames(up.mods.sin) <- up.mods.sin$Module
up.mods.sin$Module <- NULL
up.mods.sin.cal <- up.mods.sin
up.mods.sin <- as.data.frame(lapply(up.mods.sin, as.numeric))                                       # convert data frame to be numberic
up.mods.sin <- (up.mods.sin/up.mods.sin$gene)*100 
rownames(up.mods.sin) <-rownames(up.mods.sin.cal)
up.mods.sin <- up.mods.sin[,-ncol(up.mods.sin)]


#####DOWN GENE#######
down.mods.sin <- data.frame(Module = row.names(test.down), test.down,gene=0)                         # create a new blank table
down.mods.sin [,] <- NA                                                                              # create a new blank table
down.mods.sin <- down.mods.sin [-c(2:nrow(down.mods.sin)),]

for (i in 1:length(unique(Gene.matrix$Module))){
  module <- unique(Gene.matrix$Module)[i]
  sums <- colSums (test.down[Gene.matrix$Module==module,])
  genes <- nrow(Gene.matrix[Gene.matrix$Module==module,])
  down.mods.sin <- rbind(down.mods.sin,c(module,sums,genes))
}

down.mods.sin <- down.mods.sin [-1,]

rownames(down.mods.sin) <- down.mods.sin$Module
down.mods.sin$Module <- NULL
down.mods.sin.cal <- down.mods.sin
down.mods.sin <- as.data.frame(lapply(down.mods.sin, as.numeric))                                   # convert data frame to be numberic
down.mods.sin <- (down.mods.sin/down.mods.sin$gene)*100 
rownames(down.mods.sin) <-rownames(down.mods.sin.cal)
down.mods.sin <- down.mods.sin[,-ncol(down.mods.sin)]


## Prepare data for ploting ## 
Sum.mod.sin <- down.mods.sin                                                                       ## prepare a new matrix for new data
Sum.mod.sin[,] <- NA                                                                               ## Empty matrix

for (i in 1: nrow(up.mods.sin)){
  for (j in 1:ncol(down.mods.sin)){
    up = up.mods.sin[i,j]
    down = down.mods.sin[i,j]
    if (up > down) {
      res = up
    }
    if (down > up){
      res = -down
    }
    if (up == down){
      res = 0
    }
    Sum.mod.sin[i,j] = res
  }
}

################################################

save(Sum.mod.sin, sample_info ,mod_func, file = paste0("./R data/", GSE_ID, "_", platform, "_sum_mod_sin.Rdata"))
