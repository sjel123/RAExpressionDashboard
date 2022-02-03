

if (!file.exists("Meta")) {load("/app/Shiny/IndexTable2/Data/MetaAnalysis2/Meta2_0.rds")}
  Meta <- Meta1
#source("/app/Shiny/IndexTable2/R/Factor.R")
source("MetaDataPlots.R")
generow1 <- which(Meta$gene==Gene)
  p <- (PlotFunctionMeta(num = generow1, meta = Meta)[[3]])
  p1 <- (PlotFunctionMeta(num = generow1, meta = Meta)[[4]])
  
  
  
  PPTCreate( list("BTK", "BMX", "ITK", "TXK", "TEC"), Index = TRUE, RA = TRUE)
 
  Gene= list("BTK", "BMX", "ITK", "TXK", "TEC")
  Gene= list( "BMX",  "TEC")
  Gene= list("BTK", "ITK", "TXK")
  
  print(plotGeneList())
  
 plotGeneList <- function(){
   require(ggplot2)
   require(reshape2)
  Outcome <- list()
  OutcomeGtex <- list()
  OutcomeFantam <- list()
  OutcomeProtein<- list()
    for(i in Gene){
      Outcome[[i]]<-
          print(PrintFunction(i))
      OutcomeGtex[[i]]<-
        print(PrintGtexFunction(i))
      OutcomeFantam[[i]]<-
        print(PrintFantomFunction(i))
      OutcomeProtein[[i]]<-
        print( PrintProteinFunction(i))
    }

  Outcome1 <- Outcome
  InDEXList <- list(Outcome1, OutcomeGtex, OutcomeFantam, OutcomeProtein)
  for(i in 1:4 ){
    Outcome <- do.call(cbind,  InDEXList[[i]])
      Outcome <- Outcome[,c(1,2,seq(4,ncol(Outcome),2))]
  print(head(Outcome))
   Outcome.m <-  melt(Outcome)
    Outcome.m$variable <- gsub(".mean_run", "", Outcome.m$variable) 
    colnames(Outcome.m)[1] <- "BTK.Cell"

   S1<- ggplot( Outcome.m, aes(y=  BTK.Cell, x=variable, size=value, color=log(value+1))) + geom_point(alpha = 0.8) + 
     theme_classic()
   S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(0,7))
   if (i==4){ 
     S1<- ggplot( Outcome.m, aes(y=  BTK.Cell, x=variable, size=log10(value+1), color=log10(value+1))) + geom_point(alpha = 0.8) + 
     theme_classic()
      S1 = S1 + scale_color_gradient2(low = "red2",  mid="grey", high = "mediumblue", midpoint=3.5, space = "Lab", limit = c(0,7.2))
        }
   print(S1)
  }
 return(S1)
  }
 
 
 # Print Function
 # Plot Dice Data for genes of interest
 
  PrintFunction <- function(Gene=i, index=1){
    if(!exists("B_CELL_NAIVE")) {
      load("~/projects/DICE2/DICE.RData")
    }
    
    DICE_Cell_Names <- c("B_CELL_NAIVE","MONOCYTES","M2","NK",
                         "TREG_MEM", "CD4_NAIVE", "CD4_STIM", 
                         "TREG_NAIVE", "TFH", "TH1", "THSTAR", 
                         "TH17", "TH2", "CD8_NAIVE", "CD8_STIM")
    Data <- data.frame(); k=0; Data1=data.frame()
    Index <- grep(paste0("^", Gene, ";"), M2$Additional_annotations)
    if(identical(Index, integer(0))){Index=1}

    print(sprintf("Index %s", Index))
    for (j in DICE_Cell_Names){
      k=k+1
      Data[1:length(get(j)[Index,]),1]<- j
      Data[1:length(get(j)[Index,]),2]<- t(get(j)[Index,])
      Data <- Data[-1:-3,1:2 ]
      Data1<- data.frame(rbind(Data,Data1))
    }
    
    Data1 <- data.frame(Cell=Data1[,1], Value=as.numeric(Data1[,2]))
    #colnames(Data1) <- c("Cell", "Value")
    
    if(Index==1 & Gene != "TSPAN6"){Data1$Value=rep(0,nrow(Data1))}
    
    #Data1$Value <- as.numeric(as.character(Data1$Value))
    require(dplyr)
   require(tidyr)
    Data1 <- Data1 %>%
      group_by(Cell) %>%
      dplyr::summarise(mean_run = mean(as.numeric(as.character(Value))))
    
    
   # Data1 %>% group_by(Cell) %>% summarise(sum(Value))
    
    p <- ggplot(Data1, aes(x=Cell, y=mean_run))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2, height = 0)
    p <- p + ggtitle(Gene)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsize))
    p <- p + ylab("Expression Value")
    p
    return(Data1)
  }
  
  PrintGtexFunction <- function(Gene=i, index=1){
    if(!exists("Gtex")) {
      load("./Data/MembraneDataAll.RData")
    }
    index <- grep(paste0("^", Gene, "$"), Gtex$Description) 
    plotData <-data.frame(t(data.frame(Gtex[index[1],3:55])))
      plotData$Tissue <- row.names(plotData)
      colnames(plotData) <- c("Value", "Tissue")
       plotData$Tissue <- gsub("Cells...", "", plotData$Tissue)
       plotData$Tissue <- gsub("...", "", plotData$Tissue, fixed = T)
       plotData$Tissue2 <- strtrim(plotData$Tissue, 20)

    
    Data1 <- data.frame(Cell=plotData$Tissue, Value=plotData$Value)
    #colnames(Data1) <- c("Cell", "Value")
    
    if(Index==1 & Gene != "TSPAN6"){Data1$Value=rep(0,nrow(Data1))}
    
    #Data1$Value <- as.numeric(as.character(Data1$Value))
    require(dplyr)
    require(tidyr)
    Data1 <- Data1 %>%
      group_by(Cell) %>%
      dplyr::summarise(mean_run = mean(as.numeric(as.character(Value))))
    
    
    # Data1 %>% group_by(Cell) %>% summarise(sum(Value))
    
    p <- ggplot(Data1, aes(x=Cell, y=mean_run))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2, height = 0)
    p <- p + ggtitle(Gene)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsize))
    p <- p + ylab("Expression Value")
    p
    return(Data1)
  } 

  PrintFantomFunction <- function(Gene=i, index=1){
    gene=Gene
    if(!exists("MeanCluster_Fantom_t")) {
      load("./Data/MembraneDataAll.RData")
    }
   
    if (length(grep(gene, Fantom$X)>0)){
      index <- grepl(paste0(gene,"$"), Fantom$X) 
      plotData <- data.frame(t(MeanCluster_Fantom_t[index,]))
        plotData[] <- lapply(plotData, as.character)
        colnames(plotData) <- (plotData[52,])
          plotData$Cell <- row.names(plotData)
          colnames(plotData) <- make.names(colnames(plotData))
          plotData$Cell <- (row.names(plotData))
        index2 <- grep("p1\\.|Cell", colnames(plotData))
        plotData <- plotData[,index2]
      if(is.null(nrow(plotData))){plotData=data.frame(Value=rep(0,length(plotData)), Cell=plotData)}
      if(!is.null(colnames(plotData))){
        colnames(plotData) <- gsub("p1\\.", "", colnames(plotData))
        index3 <- grep(paste0("^", gene, "$"), colnames(plotData) )
        index3 <- ifelse(identical(index3, integer(0)),1, index3)
        colnames(plotData)[index3] <- "Value"
      }
    }
    
    Data1 <- data.frame(Cell=plotData$Cell, Value=plotData$Value)
    
      Data1 <- Data1[-which(as.character(Data1$Cell)=="GeneSymbol"),]
    #colnames(Data1) <- c("Cell", "Value")
    
    if(Index==1 & Gene != "TSPAN6"){Data1$Value=rep(0,nrow(Data1))}
    
    #Data1$Value <- as.numeric(as.character(Data1$Value))
    require(dplyr)
    require(tidyr)
    Data1 <- Data1 %>%
      group_by(Cell) %>%
      dplyr::summarise(mean_run = mean(as.numeric(as.character(Value))))
    
    
    # Data1 %>% group_by(Cell) %>% summarise(sum(Value))
    
    p <- ggplot(Data1, aes(x=Cell, y=mean_run))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2, height = 0)
    p <- p + ggtitle(Gene)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsize))
    p <- p + ylab("Expression Value")
    p
    return(Data1)
  } 
 
  PrintProteinFunction <- function(Gene=i, index=1){
    gene=Gene
    if(!exists("ImmuneProtgset")) {
      load("~/CSI/projects/proteomics.dev/Data/Proteomics.RData")
      ImmuneProtgset <- gset
      rm(gset)  
    }
    gene <- as.character(gene)
    index <- grep(paste0("^", gene, "$"), fData(ImmuneProtgset)$Gene.names.1) 
    plotData <-(data.frame(exprs(ImmuneProtgset)[index,]))
    colnames(plotData) <- fData(ImmuneProtgset)$Gene.names.1[index]
     plotData$Cell <- pData(ImmuneProtgset)$Cell
     plotData$State <- pData(ImmuneProtgset)$State
     colnames(plotData) <- c("Value", "Cell", "State")
    #Data1$Value <- as.numeric(as.character(Data1$Value))
     
     index <- grep( gene, fData(ImmuneProtgset)$Gene.names.1) 
     plotData <-(data.frame(exprs(ImmuneProtgset)[index[1],]))
     plotData$Cell <- pData(ImmuneProtgset)$Cell
     plotData$State <- pData(ImmuneProtgset)$State
     colnames(plotData) <- c("Value", "Cell", "State")
     
    require(dplyr)
    require(tidyr)
    Data1 <- plotData %>%
      group_by(Cell) %>%
      dplyr::summarise(mean_run = mean(as.numeric(as.character(Value))))
    
    
    # Data1 %>% group_by(Cell) %>% summarise(sum(Value))
    
    p <- ggplot(Data1, aes(x=Cell, y=mean_run))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2, height = 0)
    p <- p + ggtitle(Gene)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsize))
    p <- p + ylab("Expression Value")
    p
    return(Data1)
  }
  