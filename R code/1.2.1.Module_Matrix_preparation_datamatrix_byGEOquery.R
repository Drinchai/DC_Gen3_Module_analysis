
# Setup environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DC_Gen3_Module_analysis")

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" # "GPL6106" "GPL6947"                                                 # Available platforms for each GSE dataset can be checked by running code below

# Load Module list ##
load("./R data/DC_ModuleGen3_2019.Rdata")

# Option 1 (from Rdata file downloaded and processed with GEOquery (see script 1.1))
load(paste0("./R data/", GSE_ID, "_", platform, "_data.matrix.Rdata"))
data.matrix = as.data.frame(data.matrix)
data.matrix$ID = rownames(data.matrix)

## prepare data at the gene level by aggregate value of each gene 
data.matrix$Symbol = Probe.annotation.table$Symbol[match(data.matrix$ID, Probe.annotation.table$ID)]
data.matrix = data.matrix[-which(data.matrix$Symbol == ""),]                                          # remove probes without annotated gene Symbol
data.matrix$ID = NULL
data.matrix = aggregate(data.matrix,FUN = mean,by=list(data.matrix$Symbol))                           # calculate average of each gene
data.matrix$Symbol = NULL
rownames(data.matrix) = data.matrix$Group.1
data.matrix$Group.1 = NULL

sample.info$Type = gsub(sample.info$Type,pattern = "type 2 diabetes",replacement = "T2D")             # Modify text
sample.info$Type = gsub(sample.info$Type,pattern = "Other infections",replacement = "Other")          # Modify text

# preparing data ##
head(data.matrix)               # The data need to be chceked whether they are log2 transformed or raw data.
rownames(data.matrix)           # check gene symbol
data.matrix[data.matrix<10]=10  # fillter gene that has expression value < 10 = 10

### Prepare expression matrix with module list
df1=Module_listGen3
df2=data.frame(data.matrix)
df2$Gene = rownames(df2)

#Annotate gene module to expression matrix
df.mod = merge(df1,df2,by="Gene",all=F)

nrow(df.mod)                     # Number of available module genes

rownames(df.mod) = df.mod$Module_gene
dat.mod.func.Gen3 = df.mod[,c(1:8)]
dat.mod.Gen3 = df.mod[,-c(1:8)]

table(dat.mod.func.Gen3$Module) # Check number of available genes per module

save(dat.mod.Gen3,dat.mod.func.Gen3,sample.info,file = paste0("./R data/", GSE_ID, "_", platform, "_mod_matrix.Rdata"))
