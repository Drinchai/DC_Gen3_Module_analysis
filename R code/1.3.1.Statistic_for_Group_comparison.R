# Setup environment
rm(list=ls())
# Set working directory (location on your computer)
setwd("~/DC_Gen3_Module_analysis")
## dependency
library("gtools")

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106" 

# Load Module list ##
load("./R data/DC_ModuleGen3_2019.Rdata")

# Load module gene annotation matrix 
load("./R data/GSE13015_GPL6106_mod_matrix.Rdata")

## dependecy
library("gtools")

#prepare data for analysis
###########
df_raw = as.matrix(dat.mod.Gen3)          # replace "dat.mod.Gen3" with data_matrix in raw expression data
sample_info = sample.info                 # replace "sample_info" with sample annotation 
mod_func = dat.mod.func.Gen3              # repleace "mod_func" with Gene module annotation table  

#### make sure that expression matrix and sample information are the same order
df_raw = df_raw[,rownames(sample_info)]
colnames(df_raw) == rownames(sample_info)

#############################################
# Statistic analysis ##
############################################ 
dat_log2 <- as.matrix(log(df_raw,2))      # tranform data to log 2

## prepare entry table 
group.test = c("Melioidosis", "Other")

########################
##### T test 
########################

tt_pval = data.frame(matrix(ncol = length(group.test), nrow = nrow(dat_log2)))
colnames(tt_pval) = group.test
rownames(tt_pval) = rownames(dat_log2)

# Check if rownames of sample.info and colnames of dat_log2 are in the same order before running loop below
rownames(sample.info) == colnames(dat_log2)

k=1
for (k in 1:nrow(dat_log2)) {
  signature = rownames(dat_log2)[k]
  test.table <- sample_info
  test.table$scores <- dat_log2[k,]
  i=1
  for (i in 1:length(group.test)) {
    group = group.test[i]
    T2 <- test.table[test.table$Type == group,]               # "Type"; the selected column could be changed to your interested group comparison
    T1 <- test.table[test.table$Group==c("Control"),]         # "Group"; the selected column could be changed to your interested group comparison
    if(all(T1$scores == T2$scores)){
      tt_pval[signature,group] = 1
    }else{
      tt_pval[signature,group] <- t.test(x =T1$scores,y=T2$scores,paired = FALSE)$p.value
    }
  }
}

pvalue_Group <- data.frame(tt_pval)

pvalue_Group.FDR <- apply(pvalue_Group,2,function(x) p.adjust(x,method = "fdr")) ## Apply multiple correction testing
pvalue_Group.FDR <- as.data.frame(pvalue_Group.FDR)

####################################
####calculate fold change ##
####################################

FC.group = data.frame(matrix(ncol = length(group.test), nrow = nrow(df_raw)))
colnames(FC.group) = group.test
rownames(FC.group) = rownames(df_raw)

k=1
for (k in 1:nrow(df_raw)) {
  signature = rownames(df_raw)[k]
  test.table <- sample_info
  test.table$scores <- df_raw[k,]
  for (i in 1:length(group.test)) {
    group = group.test[i]
    T2 <- test.table[test.table$Type==group,]              # "Type"; the selected column could be changed to your interested group comparison
    T1 <- test.table[test.table$Group==c("Control"),]      # "Group"; the selected column could be changed to your interested group comparison
    FC.group[signature,group] <- foldchange(mean(T2$scores),mean(T1$scores))
  }
}   

FCgroup <- data.frame(FC.group)


#############################################
# Calculate percentage of response ##
############################################ 
#logical check ##
Group.up <- (FCgroup > 0)+(pvalue_Group.FDR < 0.1) == 2          # TRUE Up gene, Both TRUE

Group.down <- (FCgroup < 0) + (pvalue_Group.FDR < 0.1) == 2      # TRUE down gene, Both TRUE


################################################
Gene.matrix <- dat.mod.func.Gen3[rownames(Group.up),]  
Gene.matrix$Module <- as.character(Gene.matrix$Module)

up.mods.group <- data.frame(matrix(ncol =2+length(group.test), nrow = 1))        # create a new blank table
colnames(up.mods.group) = c("Module", group.test, "genes")

i=1
for (i in 1:length(unique(Gene.matrix$Module))){                                    # length of module
  module <- unique(Gene.matrix$Module)[i]                                           # look for only unique module
  sums <- colSums(Group.up[Gene.matrix$Module==module,])                            # sum upgene of each column by module 
  genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,])               # sum number of gene in each module
  up.mods.group <- rbind(up.mods.group,c(module,sums,genes))                        # paste result into a new fake table
}

# Calculate percentage of genes that are upregulated in each module 
up.mods.group <-up.mods.group[-1,]
rownames(up.mods.group) <- up.mods.group$Module
up.mods.group$Module <- NULL
up.mods.group.cal <- up.mods.group
up.mods.group <- as.data.frame(lapply(up.mods.group, as.numeric))                    # convert data frame to be numeric
up.mods.group <- (up.mods.group/up.mods.group$genes)*100 
rownames(up.mods.group) <-rownames(up.mods.group.cal)
up.mods.group <- up.mods.group[,-ncol(up.mods.group)]



#####DOWN GENE#######
down.mods.group <- data.frame(matrix(ncol =2+length(group.test), nrow = 1))        # create a new blank table
colnames(down.mods.group) = c("Module", group.test, "genes")                                                                            

for (i in 1:length(unique(Gene.matrix$Module))){
  module <- unique(Gene.matrix$Module)[i]
  sums <- colSums (Group.down[Gene.matrix$Module==module,])
  genes <- nrow(dat.mod.func.Gen3[dat.mod.func.Gen3$Module==module,])
  down.mods.group <- rbind(down.mods.group,c(module,sums,genes))
}
down.mods.group<-down.mods.group[-1,]

# Calculate percentage of genes that are downregulated in each module 
rownames(down.mods.group) <- down.mods.group$Module
down.mods.group$Module <- NULL
down.mods.group.cal <- down.mods.group
down.mods.group <- as.data.frame(lapply(down.mods.group, as.numeric))
down.mods.group <- (down.mods.group/down.mods.group$genes)*100
rownames(down.mods.group) <- rownames(down.mods.group.cal)
down.mods.group <- down.mods.group[,-ncol(down.mods.group)]

## Prepare data for ploting ## 
res.mods.group <- up.mods.group[,]                # prepare a new matrix for new data
res.mods.group[,1:ncol(res.mods.group)] <- NA                          # Empty matrix

i=1
for (i in 1: nrow(up.mods.group)){
  for (j in 1:ncol(up.mods.group)){
    up = up.mods.group[i,j]
    down = down.mods.group[i,j]
    if (up > down) {
      res = up
    }
    if (down > up){
      res = -down
    }
    if (up == down){
      res = 0
    }
    res.mods.group[i,j] = res
  }
}

################################################

save(res.mods.group, file = paste0("./R data/", GSE_ID,"_", platform, "_res_mods_group.Rdata"))

