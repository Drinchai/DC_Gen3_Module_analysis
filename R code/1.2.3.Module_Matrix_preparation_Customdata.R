#set enviroment
rm(list = ls())
setwd("")

# Load Module list ##
load("DC_ModuleGen3_2019.Rdata")

### Load expression data Data
load("./GSEXXXXX_data.matrix.Rdata")

### Prepare expression matrix with module list
df1=Module_listGen3                       # This is module list annotation table
df2=data.frame(data.matrix)               # expression data (from your own datasets or from step 1)
df2$Gene = rownames(df2)

#Annotate gene module to expression matrix
df.mod = merge(df1,df2,by="Gene",all=F)   # match df1 and df2 by Gene symbol

rownames(df.mod) = df.mod$Module_gene
dat.mod.func.Gen3 = df.mod[,c(1:8)]
dat.mod.Gen3 = df.mod[,-c(1:8)]

save(dat.mod.Gen3,dat.mod.func.Gen3,sample_info,file = "./GSEXXXXX_mod_matrix.Rdata")

