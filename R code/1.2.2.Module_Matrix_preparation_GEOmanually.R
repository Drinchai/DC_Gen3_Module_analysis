#set enviroment
rm(list = ls())
setwd("~/DC_Gen3_Module_analysis")

# Load Module list ##
load("./R data/DC_ModuleGen3_2019.Rdata")

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" # "GPL6106" "GPL6947" # Available platforms for each GSE dataset can be checked by running code below

##load expression data
# Option 2 (from mamnully downloaded txt file; Supplementary table 1: "GSE13015-GPL6106_series_matrix.txt")
data.matrix = read.table(file = paste0("./Dataset/", GSE_ID, "-", platform, "_series_matrix.txt"),sep="\t",header = TRUE,stringsAsFactors = FALSE)
rownames(data.matrix)= data.matrix$ID_REF
data.matrix$ID_REF = NULL

## Probe annotation could be manully dowload from GEO 
# "GPL6106":  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6106
# "GPL6947": https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6947 

Probe.annotation.table <- read.delim(paste0("./Dataset/Probe_annotation/", platform, ".txt"))  # "GPL6106": Supplementary table 2
rownames(Probe.annotation.table)=Probe.annotation.table$ID
Probe.annotation.table = Probe.annotation.table[rownames(data.matrix),]
rownames(Probe.annotation.table)==rownames(data.matrix)
data.matrix$Symbol = Probe.annotation.table$Symbol

## Sample annotation could be downloaded from series matrix; Supplementary table 2: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13015/matrix/
sample.info = read.csv(file = "./Dataset/Clinical Information/GSE13015_GPL6106sample_info.csv",stringsAsFactors = F)
rownames(sample.info)= sample.info$GEOID
rownames(sample.info)==colnames(data.matrix)
sample.info$Group = sapply(strsplit(sample.info$illness,"/"),"[",2)
sample.info$illness_group = sapply(strsplit(sample.info$illness,"/"),"[",1)

sample.info$Group = gsub(sample.info$Group,pattern = "type 2 diabetes",replacement = "T2D")
sample.info$Group = gsub(sample.info$Group,pattern = "Other infections",replacement = "Other")

## prepare data at the gene level by aggregate value of each gene 
data.matrix = aggregate(data.matrix,FUN = mean,by=list(data.matrix$Symbol))
data.matrix$Symbol = NULL
rownames(data.matrix) = data.matrix$Group.1
data.matrix$Group.1 = NULL

# check data ##
head(data.matrix)                          # The data need to be chceked whether they are log2 transformed or raw data.
rownames(data.matrix)                      # check gene symbol
data.matrix[data.matrix<10]=10             # fillter gene that has expression value < 10 = 10

### Prepare expression matrix with module list
df1=Module_listGen3
df2=data.frame(data.matrix)
df2$Gene = rownames(df2)

#Annotate gene module to expression matrix
df.mod = merge(df1,df2,by="Gene",all=F)

rownames(df.mod) = df.mod$Module_gene
dat.mod.func.Gen3 = df.mod[,c(1:8)]
dat.mod.Gen3 = df.mod[,-c(1:8)]

save(dat.mod.Gen3,dat.mod.func.Gen3,sample.info,file = "./R data/GSE13015_GPL6106_mod_matrix.Rdata")
