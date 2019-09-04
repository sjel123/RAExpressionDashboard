

GSM <- list.files("Data/GSE118829/",pattern = "GSM*")

multmerge = function(mypath){
  filenames=list.files(path=mypath,pattern = "GSM*",full.names=TRUE)
  datalist = lapply(filenames, function(x){read.table(file=x,header=T, sep="\t")})
  Reduce(function(x,y) {merge(x,y, by="Gene")}, datalist)
}  


mymergeddata = multmerge("Data/GSE118829/")
names(mymergeddata) <- c("Gene", GSM)


##################
#Create Expression set
##################
gset <- new("ExpressionSet", exprs = ((as.matrix(mymergeddata[,2:ncol(mymergeddata)]))))
fData(gset) <- mymergeddata[,1:2]
#Add max value per gene
fData(gset)$maxExpr <- apply(exprs(gset), 1, max)

Annotation <- names(mymergeddata[2:ncol(mymergeddata)])
Annotation[1:69] <- gsub("_CD" ,"__CD",Annotation[1:69])
Annotation <- data.frame(t(as.data.frame(sapply(Annotation, function(x) strsplit(x, "_")))))
Annotation$SubjectID <-  as.character(c(as.character(Annotation$X2[1:69]), as.character( Annotation$X3[70:336])))
Annotation$SubjectID <- as.numeric(gsub("MTX|NT|HC|IFX|TCZ|SF", "", Annotation$SubjectID))
Annotation$X2 <- gsub("[0-9]", "", Annotation$X2)
Annotation$X3 <- gsub("[0-9]", "", Annotation$X3)
  names(Annotation) <- c("ID", "Disease","Treatment", "Cell", "SubjectID")
Annotation$Cell <- gsub(".txt", "", Annotation$Cell)
  
  pData(gset) <- data.frame(Annotation)
  #Test Row names and column names match  
  table(row.names(pData(gset)) == colnames(exprs(gset)))
  names(fData(gset))[1] <- "gene_name"
  fData(gset) <- fData(gset)[,-2]
  ####If gene_name doesnt exist that add row_name to gene_name variable
  Val <- grep("gene_name", names(fData(gset)))
  if(identical(Val, integer(0))) {fData(gset)$gene_name <- row.names(fData(gset))}
  ##################
  #Study Design
  ##################
  textplot(table(pData(gset)$Group))
  textplot(table(pData(gset)$Group, pData(gset)$TMT))
  title("Study Design",xlab = "Time_point", ylab = "Treatment")
  
  ##############################
  ####ToDO ADD ENTREZ IDS#######
  ##############################
  
  
  ##################
  #Normalize Data
  ##################
  gset$Group <- paste(gset$Disease, gset$Treatment, gset$Cell, sep="_")
  Factor <- factor(gset$Group)
  design <- model.matrix(~0+Factor)
  
  NormalizeData <- function(data, design=design){
    dge <- DGEList(counts=data)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot=TRUE)
    return(v)
  }
  
  DF_Norm <- NormalizeData((exprs(gset)), design)
  DFS<- DF_Norm$E
  DFSMAX <- DFS[apply(DFS, 1, max) >4,]
  
  
  ###########################
  ########Update GSET with Normalized Values
  ###########################
  ##################
  #Create Expression set
  ##################
  gset1 <- new("ExpressionSet", exprs = ((as.matrix(DFS))))
  fData(gset1) <-   fData(gset)
  pData(gset1) <- pData(gset)
  gset <- gset1
  gsetMax <- gset[fData(gset)$maxExpr>Cutoff]
  
  
  
  ###########################
  ########QC
  ###########################
  source ("RNASeqQC2/QCFunction.R")
  source ("RNASeqQC2/QCFunction.R")
  source ("RNASeqQC2/ClusteringFunctions.R")
  source ("RNASeqQC2/RegulatedGenesFunction.R")
  
  #########SampleAnnotation
  library(knitr)
  kable(pData(gset))
  
  ###########Scatterplot
  Sample1 <- colnames(exprs(gset))[1]
  Sample2 <- colnames(exprs(gset))[2]
  ScatterFunction (exprs(gset),Sample1,Sample2)
  ScatterFunction (exprs(gset)[fData(gset)$maxExpr>Cutoff,],Sample1,Sample2)
  ScatterFunctionMatrix (exprs(gset)[fData(gset)$maxExpr>Cutoff,])
  
  ##########BoxPlot
  par(mfrow=c(1,1))
  boxplot(exprs(gset))
  boxplot(exprs(gset)[fData(gset)$maxExpr>Cutoff,])
  
  #Histograms
  HistFunc(exprs(gset))
  HistFunc(exprs(gset)[fData(gset)$maxExpr>Cutoff,])
  
  #PCA PLOT
  PCAFunction(gset[fData(gset)$maxExpr>Cutoff,],labelCol="none",  shape1 ="Cell", color1="Disease")
  PCAFunction(gset[fData(gset)$"Limmaanova" <0.001],labelCol="none",  shape1 ="Group", color1="Group")
  
  #Correlation heatmap
  
  CorFunction(gset[fData(gset)$maxExpr>Cutoff,], orderby = "Group", Labels="Group")
  CorFunction(gset[fData(gset)$"Limmaanova" <0.001&fData(gset)$"gene_type"=="protein_coding",], orderby = "Methyl", Label="Group")
  CorFunction(gset[fData(gset)$"Limmaanova" <0.001,], orderby = "Group", Label="Group")
  
  ###Density Plots
  #Histograms
  DensityPlotFunction (gset)
  DensityPlotFunction (gsetMax)
  DensityPlotFunction (DFSMAX)
  
  ## Identification of differential genes
  #####Limma Aanlysis to identify differential expressed genes
  CONtrasts <- c("rheumatoid.arthritis - healthy",
                 "osteoarthritis - healthy",
                 "rheumatoid.arthritis - osteoarthritis")
  
  
  
  efit1 <- LimmaFunction()[[1]] 
  cont.wt1 <- LimmaFunction()[[2]]
  save(efit1, cont.wt1, file="Data/efit1.RData")    
  ############
  
  topTable(efit1, coef=1)
  topTable(efit1, coef=2)
  topTable(efit1, coef=3)
  topTable(efit1, coef=4)
  maxExpr<- fData(gset)$maxExpr>0
  plot(x=topTable(efit1, coef=1, n="inf")[maxExpr,10],
       y=topTable(efit1, coef=2, n="inf")[maxExpr,10])
  
  
  #Volanoplot
  volcanoplot(efit1, coef=1,  highlight=20, names=efit1$genes$gene_name)
  volcanoplot(efit1, coef=1,  highlight=20, names=efit1$genes$parent)
  #MvsA plot
  limma::plotMA (efit1,coef=1, xlab="Average log-Expression", ylab= "Log Fold Change")
  o <- order(efit1$p.value[,1])
  x <- efit1$Amean
  y <- efit1$coefficients[,1]
  G <- efit1$genes$parent
  
  text(x[o[1:20]], y[o[1:20]], labels=G[o[1:20]])
  
  
  ###Combine limma output pvalue into single table  
  
  test1 <- CombineLimmaOutput(efit = efit1, CONTRASTs = cont.wt1)
  
  #### Add output to gset
  
  Index <- which(names(fData(gset))=="maxExpr")
  fData(gset) <- data.frame(fData(gset)[1:Index], as.data.frame(test1))
  #### Save gset
  TcellTretgset <- gset

  
  ###
  #Regulated Genes Table
  ###
  source("RNASeqQC2/RegulatedGenesFunction.R")
  RegGeneTable(maxExpr = 3.0, pVal= 0.0001, FC=1.4)
  
  #Kmeans clustering
  
  KmeansOutput <- KmeansClusterFun(gset,n=5, shade=50, ESET=gset, orderBy="Group",PCutoff=0.00005) 
  KmeansOutput <- KmeansClusterFun(gsetPC,n=5, shade=50, ESET=gset, orderBy="TMT",PCutoff=0.0005101) 
  KmeansOutput <- KmeansClusterFun(gset,n=5, shade=50, ESET=gset, orderBy="Group",PCutoff=0.05) 
  
  ###################################################

  gset$Treatment[gset$Treatment==""] <- "none"
  gset$Treatment<- factor(gset$Treatment, levels = c("none", "NT", "IFX", "MTX", "SF", "TCZ"))
  TcellTretgset <- gset
  
  save (TcellTretgset, file="~/app/RAExpressionDashboard/Data/TCellgset")
  
PlotTcellTretData <- function(gene="CDH11"){
    if(!exists("TcellTretgset")) {
      load("~/app/RAExpressionDashboard/Data/TCellgset")
  
    }
    gene <- as.character(gene)
    
    PlotData1 <- exprs(TcellTretgset)[fData(TcellTretgset)$gene_name == gene,]
    PlotData_t <-data.frame(Value=PlotData1, pData(TcellTretgset))
    
    my_comparisons <- list( c("none", "NT"),
                            c("NT", "IFX"), 
                            c("NT", "MTX"),
                         
                            c("NT", "TCZ")
    )
    maxval <- max(PlotData_t$Value)
    a <-  ggplot(PlotData_t,aes(x=Treatment,y=Value,fill = Disease)) +
      geom_boxplot(outlier.shape = NA)+geom_jitter( position = position_jitterdodge(), size=3)+
      labs(x='', y = 'Expression', title=paste0("Treated T cell Expression Data ",gene))+
      theme_bw()+ facet_wrap(~Cell)+
      theme(legend.position="none")+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
    a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = maxval+6)     # Add global p-value


    return(a)
  } 
  
  VART <- apply(exprs(TcellTretgset), 1,function(x) var(x))
  
  fData(TcellTretgset)[names(VART[order(VART,decreasing = TRUE)][1:10]),]
  
  PlotTcellTretData("SOCS3")
  
  
  ## Identification of differential genes
  #####Limma Aanlysis to identify differential expressed genes
  CONtrasts <- c("RA_IFX_CD4TN - RA_NT_CD4TN",
                 "RA_MTX_CD4TN - RA_NT_CD4TN",
                 "RA_TCZ_CD4TN - RA_NT_CD4TN")
  
  
  
  efit1 <- LimmaFunction()[[1]] 
  cont.wt1 <- LimmaFunction()[[2]]
  save(efit1, cont.wt1, file="Data/efit1.RData")    
  ############
  
  topTable(efit1, coef=1)
  topTable(efit1, coef=2)
  topTable(efit1, coef=3)
  