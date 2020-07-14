# Setup environment
rm(list = ls())
setwd("~/DC_Gen3_Module_analysis")

library(preprocessCore)
# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" # "GPL6106" "GPL6947" # Available platforms for each GSE dataset can be checked by running code below

# Load Module list ##
load("./R data/DC_ModuleGen3_2019.Rdata")

# Option 1 (from Rdata file downloaded and processed with GEOquery (see script 1.1))
load(paste0("./R data/", GSE_ID, "_", platform, "_data.matrix.Rdata"))
data.matrix = as.data.frame(data.matrix)
data.matrix$ID = rownames(data.matrix)

##quantile normalization
data.matrix.nor <- normalize.quantiles(as.matrix(data.matrix))
colnames(data.matrix.nor) = colnames(data.matrix)
rownames(data.matrix.nor) = rownames(data.matrix)

## prepare data at the gene level by aggregate value of each gene 
data.matrix.nor$Symbol = Probe.annotation.table$Symbol[match(data.matrix.nor$ID, Probe.annotation.table$ID)]
data.data.matrix.nor = data.matrix.nor[-which(data.matrix.nor$Symbol == ""),]                                          # remove probes without annotated gene Symbol
data.matrix.nor$ID = NULL
data.matrix.nor = aggregate(data.matrix.nor,FUN = mean,by=list(data.matrix.nor$Symbol))                           # calculate average of each gene
data.matrix.nor$Symbol = NULL
rownames(data.matrix.nor) = data.matrix.nor$Group.1
data.matrix.nor$Group.1 = NULL

sample.info$Type = gsub(sample.info$Type,pattern = "type 2 diabetes",replacement = "T2D")             # Modify text
sample.info$Type = gsub(sample.info$Type,pattern = "Other infections",replacement = "Other")          # Modify text

# preparing data ##
head(data.matrix.nor)                   # The data need to be chceked whether they are log2 transformed or raw data.
rownames(data.matrix.nor)               # check gene symbol
data.matrix.nor[data.matrix.nor<10]=10  # fillter gene that has expression value < 10 = 10

### Prepare expression matrix with module list
df1=Module_listGen3
df2=data.frame(data.matrix.nor)
df2$Gene = rownames(df2)

#Annotate gene module to expression matrix
df.mod = merge(df1,df2,by="Gene",all=F)

nrow(df.mod)                     # Number of available module genes

rownames(df.mod) = df.mod$Module_gene
dat.mod.func.Gen3 = df.mod[,c(1:8)]
dat.mod.Gen3 = df.mod[,-c(1:8)]

table(dat.mod.func.Gen3$Module) # Check number of available genes per module

save(dat.mod.Gen3,dat.mod.func.Gen3,sample.info,file = paste0("./R data/", GSE_ID, "_", platform, "_mod_matrix.Rdata"))
