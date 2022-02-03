

# SSRA <- read.table("Data/GSE109449_singlecell_rnaseq_gene_tpm.tsv", sep="\t", header=T)
# annotation <- read.table("Data/GSE109449_singlecell_rnaseq_metadata.tsv", sep="\t", header=T)
# 
# annotation$sample_name==colnames(SRA.m)[3:386]
# 



# library("org.Hs.eg.db") # remember to install it if you don't have it already
# symbols <- mapIds(org.Hs.eg.db, keys = as.character(SSRA$ID_REF), keytype = "ENSEMBL", column="SYMBOL")
# SYMbols <- data.frame(gene.name = symbols)
# 
# SSRA.m <- merge(SSRA, SYMbols, by.x="ID_REF", by.y=0, all.x=0)
# 
# SRA.m <- SSRA.m[,c(1,386, 2:385)]


##################
#Create Expression set
# ##################
# SSRA <- new("ExpressionSet", exprs = ((as.matrix(SRA.m [,3:ncol(SRA.m)]))))
# fData(SSRA) <- SRA.m[,1:2]
# #Add max value per gene
# fData(SSRA)$maxExpr <- apply(exprs(SSRA), 1, max)
# pData(SSRA) <- data.frame(annotation)


PlotSSRA <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("SSRA")) {
    load("~/projects/CDH11_membrance/CDH11/Data/BrennerFib.RData")
    SSRA <- gset_RNAseq
    fData(SSRA)$gene_name <- fData(SSRA)$external_gene_name
  
  }
  gene <- as.character(gene)

  PlotData <- exprs(SSRA)[which(fData(SSRA)$gene_name==gene),]
    
  PlotData_t <-data.frame(Value=PlotData, pData(SSRA))
     #PlotData_t$CD34 <- PlotData_t$CD34_protein
     ##PlotData_t$THY1 <- PlotData_t$THY1_protein
     #PlotData_t$CD11 <- PlotData_t$CDH11_protein
       PlotData_t$CD34 <-gsub("+", "CD34H",PlotData_t$CD34, fixed = T )
       PlotData_t$CD34 <-gsub("-", "CD34L",PlotData_t$CD34, fixed = T )
       PlotData_t$THY1 <-gsub("+", "THY1H",PlotData_t$THY1, fixed = T )
       PlotData_t$THY1 <-gsub("-", "THY1L",PlotData_t$THY1, fixed = T )
       PlotData_t$CDH11 <-gsub("+", "CDH11H",PlotData_t$CDH11, fixed = T )
       PlotData_t$CDH11 <-gsub("-", "CDH11L",PlotData_t$CDH11, fixed = T )
     PlotData_t$group <- paste(PlotData_t$CD34, PlotData_t$THY1, PlotData_t$CDH11, sep="_")

    PlotData1 <- exprs(SSRA)[fData(SSRA)$gene_name %in% c("PDPN", "FAP", "THY1"),]

  PlotData_t1 <-data.frame(t(PlotData1))
  
  colnames(PlotData_t1) <- fData(SSRA)[fData(SSRA)$gene_name %in% c("PDPN", "FAP", "THY1"),2]

  PlotData_t <- data.frame(PlotData_t, PlotData_t1)
  PlotData_t$FABGroup [PlotData_t$FAP>200 & PlotData_t$PDPN>75]<- as.character(
    PlotData_t[PlotData_t$FAP>200 & PlotData_t$PDPN>75,"group"])
 
  PlotData_t$FABGroup <- gsub("CD34[HL]_THY1L_CDH11[HL]", "FAPH_PDPNH_THY1L",PlotData_t$FABGroup )
  PlotData_t$FABGroup <- gsub("CD34[HL]_THY1H_CDH11[HL]", "FAPH_PDPNH_THY1H",PlotData_t$FABGroup )
  
  my_comparisons <- list( c("FAPH_PDPNH_THY1H" , "FAPH_PDPNH_THY1L"))
  maxval <- max(PlotData_t$Value)+1
  a <-  ggplot(PlotData_t,aes(x=group,y=Value)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter(  size=3)+
    labs(x='', y = 'Expression', title=paste0("Single Cell Fibroblast Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
 a <- a + geom_boxplot(data= PlotData_t[PlotData_t$FABGroup %in% c("FAPH_PDPNH_THY1L" , "FAPH_PDPNH_THY1H"),],aes(x=FABGroup, y=Value),  outlier.shape = NA)+ 
   geom_jitter(data= PlotData_t[PlotData_t$FABGroup %in% c("FAPH_PDPNH_THY1L" , "FAPH_PDPNH_THY1H"),],aes(x=FABGroup, y=Value),size=3)+
  geom_vline(xintercept = 7.5)
 a <- a + 
 guides(shape = FALSE)+
   theme(legend.position=c(0,1))+guides(fill = guide_legend(override.aes = list(size = 0)))
   return(a)
} 


#PlotSSRA(gene="FAP")
