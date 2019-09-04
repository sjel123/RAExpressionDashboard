#####Load GSE124731 T cell innateness profiling

DF <- read.table("Data/GSE124731_low_input_rnaseq_gene_normalized.txt", header=TRUE)
ANN <- read.table("Data/GSE124731_low_input_rnaseq_meta_data.txt", header=TRUE)


##################
#Create Expression set
##################
gset <- new("ExpressionSet", exprs = ((as.matrix(DF[,2:ncol(DF)]))))
fData(gset) <- data.frame(DF[,1])
#Add max value per gene
fData(gset)$maxExpr <- apply(exprs(gset), 1, max)
pData(gset) <- data.frame(ANN)
row.names(pData(gset)) <- pData(gset)$sampleID

#Test Row names and column names match  
table(row.names(pData(gset)) == colnames(exprs(gset)))


library("org.Hs.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Hs.eg.db, keys = as.character(SSRA$ID_REF), keytype = "ENSEMBL", column="SYMBOL")
SYMbols <- data.frame(gene.name = symbols)

Merg <- merge(fData(gset), SYMbols, by.x=1, by.y=0, all.x=TRUE)

table(fData(gset)[,1]==Merg$DF...1.)
fData(gset) <- Merg


##################
#Study Design
##################
gset$Group <- gset$cell_type
  gset$Group <- factor(gset$Group)
   levels(gset) <- (c("CD4", "CD8",  "MAIT", "iNKT", "Vd1", "Vd2", "NK"))

  InnateGradient <- gset
  save (InnateGradient, file = "Data/InnateCellGradient.RData")
 
   PlotTcellInnate <- function(gene="CDH11"){
    if(!exists("InnateGradient")) {
      load("~/app/RAExpressionDashboard/Data/InnateCellGradient.RData")
      
    }
    gene <- as.character(gene)
    fData(InnateGradient)$gene_name <- fData(InnateGradient)$gene.name
    PlotData1 <- exprs(InnateGradient)[fData(InnateGradient)$gene_name == gene,]
    PlotData1 <- PlotData1[!is.na(row.names(PlotData1)),]
    PlotData_t <-data.frame(Value=PlotData1, pData(InnateGradient))
    PlotData_t$Group <- factor(PlotData_t$Group, levels(PlotData_t) <- 
                                 (c("CD4", "CD8",  "MAIT", "iNKT", "Vd1", "Vd2", "NK")))
    my_comparisons <- list( c("CD4", "CD8"),
                            c("CD4", "MAIT"),
                            c("CD4", "iNKT"), 
                            c("CD4", "Vd1"),
                            c("CD4", "Vd2"),
                            c("CD4", "NK")
    )
    maxval <- max(PlotData_t$Value)
    a <-  ggplot(PlotData_t,aes(x=Group,y=Value, fill=Group)) +
      geom_boxplot(outlier.shape = NA)+geom_jitter( position = position_jitterdodge(), size=3)+
      labs(x='', y = 'Expression', title=paste0("T cell Innate Gradient Expression Data ",gene))+
      theme_bw()+ 
      theme(legend.position="none")+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  
      a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = maxval+3)     # Add global p-value
    
    a
    return(a)
  } 
  
   PlotTcellInnate("DTHD1")
gene="TRDJ1"
