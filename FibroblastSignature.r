# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Aug 9 13:11:56 EDT 2019


#Server: www.ncbi.nlm.nih.gov
#Query: acc=GSE63626&platform=GPL570&type=txt&groups=All&colors=dfeaf4&selection=000000000000000000000000000000000000000000000000000000000000000&padj=fdr&logtransform=auto&columns=ID&columns=adj.P.Val&columns=P.Value&columns=F&columns=Gene+symbol&columns=Gene+title&num=250&annot=ncbi

# Unable to generate script analyzing differential expression.
#      Invalid input: at least two groups of samples should be selected.

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE63626", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

library(hgu133plus2.db)

select(hgu133plus2.db, c("1007_s_at","1053_at"), c("SYMBOL", "GENENAME")) ##  This is just a trying example
PROBES<- as.character(row.names(exprs(gset)))
OUT <- select(hgu133plus2.db,keys= PROBES, columns=c("SYMBOL", "GENENAME"))

table(table(OUT$PROBEID))
OUT <- OUT[!(duplicated(OUT$PROBEID)),]

table(row.names(exprs(gset))==OUT$PROBEID)
fData(gset) <- OUT

 fibrgset <- gset
save(fibrgset, file="Data/Fibroblastgset.RData")


PlotFIBData <- function(gene="CDH11"){
  if(!exists("fibrgset")) {
    load("~/app/RAExpressionDashboard/Data/Fibroblastgset.RData")
    fData(fibrgset)$gene_name <- fData(fibrgset)$SYMBOL
  }
  gene <- as.character(gene)
  
  PlotData1 <- exprs(fibrgset)[fData(fibrgset)$gene_name == gene,]
  PlotData1 <-   PlotData1 [!is.na(row.names(PlotData1)),]
       Index=1
      if(!is.null(nrow(PlotData1))) {MAX <- apply(PlotData1, 1, max)
                              Index <- which(MAX==max(MAX))}
      
      if(!Index==1){PlotData1 <- PlotData1 [Index,]}
  PlotData_t <-data.frame(Value=PlotData1, pData(fibrgset))
  
  PlotData_t$group <- factor(PlotData_t$description)
  PlotData_t$group <-gsub("Normal human | fibroblasts","", PlotData_t$group)
  # Add groups for FAB+PDPN+ thy1 +/-
  PlotData1 <- exprs(fibrgset)[fData(fibrgset)$gene_name %in% c("PDPN", "FAP", "THY1"),]
  
  PlotData_t1 <-data.frame(t(PlotData1))
  
  colnames(PlotData_t1) <- fData(fibrgset)[fData(fibrgset)$gene_name %in% c("PDPN", "FAP", "THY1"),"gene_name"]
  
  PlotData_t <- data.frame(PlotData_t, PlotData_t1)

  a <-  ggplot(PlotData_t,aes(x=group,y=log(Value),fill = tissue.ch1)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter( position = position_jitterdodge(), size=3)+
    labs(x='', y = 'Expression', title=paste0("Single Cell Fibroblast Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())

  return(a)
} 

PlotFIBData("FAP")


VAR <- apply(exprs(fibrgset),1,var)

which(row.names(fibrgset)=="210809_s_at")

