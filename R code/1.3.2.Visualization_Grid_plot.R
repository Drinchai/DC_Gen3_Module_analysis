# Setup environment
rm(list=ls())

# Set working directory (location on your computer)
setwd("~/Dropbox (TBI-Lab)/DC_Gen3_Module_analysis")

library(reshape2)
library(ggplot2)

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106"

# Load data
load(paste0("./R data/", GSE_ID, "_", platform, "_res_mods_group.Rdata"))
load("./R data/DC_ModuleGen3_2019.Rdata")

## prepared cluter position
Group_plot = res.mods.group
Group_plot <-Group_plot[rownames(Gen3_ann),]
rownames(Group_plot)==rownames(Gen3_ann)                         # check if rownames is the same 
rownames(Group_plot) <- Gen3_ann$position
Group_plot <- as.data.frame(Group_plot)

head(Group_plot)

# creat new grid with all filtered cluster##
mod.group1 <- matrix(nrow=38,ncol=42)       
rownames (mod.group1) <- paste0("A",c(1:38))
colnames (mod.group1) <- paste0("",c(1:42))
##

diseases = colnames(Group_plot)
N.disease = length(diseases)

i=1
for (i in 1:N.disease){
  disease = diseases[i]
  for (i in 1 : nrow(Group_plot)){
    Mx <- as.numeric(gsub(x = strsplit (rownames(Group_plot)[i],"\\.")[[1]][[1]],pattern = "A",replacement = ""))
    My <- as.numeric(strsplit (rownames(Group_plot)[i],"\\.")[[1]][[2]])
    mod.group1[Mx,My] <- Group_plot[,disease][i] 
  }
  mod.group <- mod.group1[-c(9:14,19:23),]
  melt_test <- melt(mod.group,id.var=c("row.names"))
  colnames(melt_test) = c("Aggregate","Sub_aggregate","%Response")
  pdf(paste0("./Figure/Group comparison/", GSE_ID, "_", platform, "_Group_comparison_to_healthy", disease, "_Grid.pdf"), height = 5.5, width = 8.5)
  plot = ggplot(melt_test, aes(Aggregate, as.factor(Sub_aggregate))) +
    geom_tile(color="#E6E6E6" , size = 0.2, fill=color )+
    geom_point(aes(colour=`%Response`),size=4.5)+ 
    ylab("") +
    xlab("") +
    labs(title= disease)+
    theme(axis.text.x = element_text(angle = -90, hjust = 0))+
    scale_color_gradient2(low = "blue", mid="white", high = "red",limits=c(-100,100), na.value = "#E6E6E6", guide = "colourbar")+
    theme_light() +
    theme(panel.grid.minor = element_line(colour="black", size=0.9))+
    coord_flip() + 
    scale_x_discrete(limits = rev(levels(melt_test$Aggregate))) +
    theme(panel.border = element_rect(color = "black",size = 0.5),
          axis.text.x = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=2,face="plain"),
          axis.text.y = element_text(colour="black",size=9,angle=0,hjust=0.5,vjust=0.5,face="plain"))
  
  plot(plot)
  dev.off()
}

