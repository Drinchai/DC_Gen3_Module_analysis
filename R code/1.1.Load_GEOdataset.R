# Setup environment
rm(list=ls())

# Set working directory (location on your computer)
setwd("~/DC_Gen3_Module_analysis")

#Load package
library(GEOquery)

# GSE13015 will be used as an example datasets. In this dataset was included by two different platform "GPL6106" and "GPL6947" in a GSE file.
# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" # "GPL6106" "GPL6947" # Available platforms for each GSE dataset can be checked by running code below

#GET GSE soft file and matrix
Dat <- getGEO(GSE_ID, GSEMatrix=FALSE)

#check sample names
names(GSMList(Dat))

#platforms used in this GSE
names(GPLList(Dat))

#show platforms by sample
GSM.platforms <- lapply(GSMList(Dat ),function(x) {Meta(x)$platform}) 
df = data.frame(GSM.platforms)

# get GSM ID for specific platform
GSM_IDs = colnames(df)[which(df == platform)]

#example of an GSM expression vector 
Table(GSMList(Dat)[[1]])[1:100,]

# Get gene annotation data from specified platform
Probe.annotation.table <- Table(GPLList(Dat)[[platform]])[,1:10]

#Probeset extrated from GPL of GSM 1 
probesets <- Table(GPLList(Dat)[[platform]])$ID                                          # This can be different in different platform (illumina or Affy)
rownames(Probe.annotation.table) <- make.names(Probe.annotation.table$ID, unique=TRUE)


#creating the expression matrix ordered by the GPL order of probes
data.matrix <- do.call('cbind',lapply(GSMList(Dat),function(x) {
  tab <- Table(x)
  mymatch <- match(probesets,tab$ID_REF)
  return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})


rownames(data.matrix) <- probesets                                                     # give rowname = probsets 
data.matrix <- data.matrix[,colSums(is.na(data.matrix))<nrow(data.matrix)]             # get rid of Na column, this will remove samples of other GPL platform
data.matrix <- data.matrix[complete.cases(data.matrix),]                               # remove probes without signal
data.matrix[1:5,]

### Match probe and gene names
rownames(Probe.annotation.table) = gsub(rownames(Probe.annotation.table),pattern = "X",replacement = "")
ProbeID.Dat <- Probe.annotation.table[which(rownames(Probe.annotation.table)%in%rownames(data.matrix)),]
rownames(ProbeID.Dat)==rownames(data.matrix)
rownames(data.matrix) <- ProbeID.Dat$ID

##########################
##Loading Annotation data
##########################

file_name = paste0(GSE_ID, "-", platform, "_series_matrix.txt.gz")

dat_info  <- getGEO(GSE_ID,GSEMatrix=TRUE)
Phenotipic.data <- pData(dat_info[[file_name]])   # file_name could be replaced by "1", in case the GSE# has only one platform 
Phenotipic.characteristics <- Phenotipic.data[,grep(x = colnames(Phenotipic.data), pattern = "characteristics")]
Phenotipic.characteristics <- data.frame(lapply(Phenotipic.characteristics, as.character), stringsAsFactors=FALSE)

#Fix colnames
fix.col.names <- t(as.data.frame(strsplit(x = as.character(Phenotipic.characteristics[1,]),split = ":")))
colnames(Phenotipic.characteristics) <- fix.col.names[,1]

#Fix data
for (i in 1:ncol(Phenotipic.characteristics)) {
  x<-colnames(Phenotipic.characteristics[i])
  x<- paste0(x,": ")
  print (x)
  x<-gsub(x=x ,pattern = "\\(",replacement = "\\\\(")
  x<-gsub(x=x ,pattern = "\\)",replacement = "\\\\)")
  Phenotipic.characteristics[,i] <- gsub(x = Phenotipic.characteristics[,i],pattern = x,replacement = "")
}

rownames(Phenotipic.characteristics) <- rownames(Phenotipic.data)

sample.info <- Phenotipic.characteristics

## Clean up some column in this datasets for future analysis
sample.info$Illness = gsub(sample.info$Illness,pattern = "/",replacement = "_")
sample.info$Group = sapply(strsplit(sample.info$Illness,"_",fixed = TRUE),"[",1)
sample.info$Type = sapply(strsplit(sample.info$Illness,"_",fixed = TRUE),"[",2)

dir.create("./R data", showWarnings = FALSE)

save(data.matrix,sample.info, Probe.annotation.table, file = paste0("./R data/", GSE_ID, "_", platform, "_data.matrix.Rdata"))

