# #####Load GSE118829
# library(limma)
# library(ggplot2)
# library(edgeR)
# library(reshape)
# #library(plotly)
# library(plyr)
# library(Biobase)
# library(gplots)
# library(preprocessCore)
# 
# i = list.files("Data/GSE118829/",pattern = "GSM*")[1]
# DF2 <- read.table(paste0( "Data/GSE118829/",i), header =TRUE, sep = "\t")
#   colnames(DF2)[2] <- i
# for (i in list.files("Data/GSE118829/",pattern = "GSM*")) {print (i)
# 
# DF1 <- read.table(paste0( "Data/GSE118829/",i), header =TRUE, sep = "\t")
# DF2 <- data.frame(DF2, DF1[,2])
#   colnames(DF2)[length(DF2)] <- i
# }
# 
# head(DF2)
# SampleAnnotation <- colnames(DF2)[-1]
#     SampleAnnotation1 <-  data.frame(t(data.frame(strsplit(SampleAnnotation[1:70], split = "_"))))
#       rownames(SampleAnnotation1) <- colnames(DF2)[2:71]
#       colnames(SampleAnnotation1) <- c("SampleID", "Disease", "Cell")
#       SampleAnnotation1$Treatment <- ""
#       SampleAnnotation1 <- SampleAnnotation1[,c(1,2,4,3)]
#       SampleAnnotation2 <-  data.frame(t(data.frame(strsplit(SampleAnnotation[71:337], split = "_"))))
#       rownames(SampleAnnotation2) <- colnames(DF2)[72:338]
#       colnames(SampleAnnotation2) <- c("SampleID", "Disease", "Treatment", "Cell")
#       SampleAnnotation <- rbind(SampleAnnotation1, SampleAnnotation2)
#       SampleAnnotation$Donor <- paste0(SampleAnnotation$Disease, SampleAnnotation$Treatment)
#       SampleAnnotation$Donor <- gsub('[A-Z]','', SampleAnnotation$Donor)
#       SampleAnnotation$Disease <- gsub('[0-9]','', SampleAnnotation$Disease)
#       SampleAnnotation$Treatment <- gsub('[0-9]','', SampleAnnotation$Treatment)
# 
# ##################
# #Create Expression set
# ##################
# gset <- new("ExpressionSet", exprs = ((as.matrix(DF2[,2:ncol(DF2)]))))
# fData(gset) <- DF2[,1:2]
# #Add max value per gene
# fData(gset)$maxExpr <- apply(exprs(gset), 1, max)
# pData(gset) <- data.frame(SampleAnnotation)
# gset$Cell <- gsub(".txt|.1", "", gset$Cell)
# gset$Group <- paste(gset$Disease, gset$Treatment, gset$Cell, sep=".")
# 
# #Test Row names and column names match
# table(row.names(pData(gset)) == colnames(exprs(gset)))
# 
# 
# 
# ##################
# #Normalize Data
# ##################
# Factor <- factor(gset$Group)
# design <- model.matrix(~0+Factor)
# 
# NormalizeData <- function(data, design=design){
#   dge <- DGEList(counts=data)
#   dge <- calcNormFactors(dge)
#   v <- voom(dge, design, plot=TRUE)
#   return(v)
# }
# 
# DF_Norm <- NormalizeData((exprs(gset)), design)
# DFS<- DF_Norm$E
# DFSMAX <- DFS[apply(DFS, 1, max) >4,]
# 
# 
# ###########################
# ########Update GSET with Normalized Values
# ###########################
# #assayDataElement(gset, "raw_data")<- exprs(gset)
# 
# exprs(gset) <-(DFS)
# # For Cytoreason analysis no need for cutoffs
# gsetGSE118829 <- gset
# save(gsetGSE118829, file="Data/GSE118829.rds")
# 
# Cutoff <- 4
# gsetMax <- gset[fData(gset)$maxExpr>Cutoff]
# 
# 
# boxplot(exprs(gsetMax))
# #Histograms
# HistFunc(exprs(gset))
# HistFunc(exprs(gset)[fData(gset)$maxExpr>Cutoff,])
# 
# #PCA PLOT
# PCAFunction(gset[fData(gsetMax)$maxExpr>Cutoff,],labelCol="none",  shape1 ="Treatment", color1="Cell")
# PCAFunction(gset[fData(gset)$"Limmaanova" <0.001],labelCol="none",  shape1 ="Group", color1="Group")

PlotGSE118829 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("gsetGSE118829")) {
    load("Data/GSE118829.rds")
  }
  PlotData <- data.frame(Value=exprs(gsetGSE118829 )[fData(gsetGSE118829)$Gene==gene,])
  if(dim(PlotData)[1]==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  PlotData_t <-data.frame(Value=PlotData, pData(gsetGSE118829))
  #PlotData_t$Group = PlotData_t$condition
  
  my_comparisons <- list( c("healthy", "osteoarthritis"), 
                          c("healthy", "rheumatoid arthritis"),
                          c ("osteoarthritis" , "rheumatoid arthritis"))
  maxval <- max(PlotData_t$Value)
  
  a <-  ggplot(PlotData_t,aes(x=Cell,y=Value,label = Disease,fill = Treatment)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE55335 Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
  return(a)
  
}
