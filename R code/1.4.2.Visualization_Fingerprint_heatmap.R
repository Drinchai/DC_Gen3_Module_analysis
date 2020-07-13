## Setup environment
rm(list=ls())

# Set working directory (location on your computer)
setwd("~/DC_Gen3_Module_analysis")

# Load dependencies
library(stringr)
library(ComplexHeatmap)

# Set parameters
GSE_ID = "GSE13015"
platform = "GPL6106"


# Load data
load("./R data/DC_ModuleGen3_2019.Rdata")
load(paste0("./R data/", GSE_ID, "_", platform, "_sum_mod_sin.Rdata"))

########### PREPARING DATA FOR PLOTTTING STEP ############################################
##############################################################################################################
Sum.mod.sin = Sum.mod.sin[rownames(Gen3_ann),]
rownames(Sum.mod.sin) == rownames(Gen3_ann)

rownames(Sum.mod.sin) <- paste(Gen3_ann$Module, Gen3_ann$Function_New, sep = ".")

Sum.mod.sin.comp <- Sum.mod.sin[apply(Sum.mod.sin[,], 1, function(x) !all(x==0)),]                 # Rows sum=0

############################################################
##################### MODULES GEN3 and MODULE WITH FUNCTION DEFINED #######################################
#modules with function deffined

Module.list <- unique(mod_func[,c("Module","Function")])                                                             # creat new dataframe from Module
Module.list$Modules <- paste(Module.list$Module, Module.list$Function, sep = ".")
rownames(Module.list) <- Module.list$Modules

mod.with.function <- Module.list$Modules[which(Module.list$Function!="TBD")]                         # select module that have only function
Sum.mod.sin.comp.withF <- Sum.mod.sin.comp[rownames(Sum.mod.sin.comp) %in% mod.with.function,]       # selected only modules that have function in this dataset

####################################################################################
####### DOT Heatmap by complexHeatmap ####

df_plot = Sum.mod.sin.comp.withF
df_plot = df_plot[rowSums(df_plot != 0) >20,]   # keep only rows that have values for more than 20 individuals

########## An example of DISPLAY DATA > 15 %
df_plot[abs(df_plot) < 15] <- 0

sample_info = sample_info[order(sample_info$Type),]
sample_info = sample_info[order(sample_info$Group),]

df_plot = df_plot[,rownames(sample_info)]

colnames(df_plot)==rownames(sample_info)


col_fun = circlize::colorRamp2(c(-100,0,100), c("blue", "white", "red"))

ha_column = HeatmapAnnotation(df = data.frame(Illness = sample_info$Group,
                                              Diseases = sample_info$Type), 
                              show_annotation_name = TRUE,
                              col = list(Illness = c("Control" ="#64C0FE", "Sepsis" = "#FE5386"),
                                         Diseases=c("healthy" ="#64C0FE","Melioidosis"="#FF7F00","Other"="#FEFB01","Recovery"="#113fc7","T2D"="#D6D6D6")))


#DOT HEATMAP

pdf(paste0("./Figure/Individual comparison/", GSE_ID, "_", platform, "_Module_Gen3_individual_FC1.5diff100_20perct_nocluster.pdf"), height = 27, width = 18)
ht=Heatmap(df_plot,
           cluster_rows = TRUE,
           cluster_columns = T,
           height = unit(4, "mm")*nrow(df_plot), 
           width  = unit(4, "mm")*ncol(df_plot), 
           rect_gp = gpar(type = "none"),
           top_annotation = ha_column,
           name = "% Response",
           row_names_max_width = unit(10,"in"),
           row_title_gp = gpar(fontsize = 0.1),
           column_names_gp = gpar(fontsize = 12),
           row_names_gp = gpar(fontsize = 13),
           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.circle(x = x, y = y, r = unit(1.85, "mm") ,gp = gpar(fill = col_fun(df_plot[i, j]), col = NA))
           }
)
draw(ht,heatmap_legend_side = "left", annotation_legend_side = "left", padding = unit(c(2, 20, 2, 2), "mm"))

dev.off()


