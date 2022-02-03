#####Melanocyte selective genes
#library
library(dplyr)
library(plyr)
### Load Data


proteintogene <- read.table("uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab", header=T, sep="\t", quote="", comment.char = "")
  head(proteintogene)

skinatlas <- read.table("skinatlas_data_data.tsv", header=T, sep="\t", stringsAsFactors = F)
  skinatlas <- as.data.frame(t(skinatlas),stringsAsFactors = F)
  colnames(skinatlas)  <- as.character(skinatlas[1,])
  head(skinatlas)
    skinatlas <- skinatlas[-1,]
    skinatlas$ids <- row.names(skinatlas)
    proteintoGene <- merge(data.frame(ids=skinatlas$ids), proteintogene[,c(1:2,5)],by.x="ids", by.y="Entry", all.x=T )
    skinatlas <- join(skinatlas, proteintoGene)
      skinatlas$Entry.name <- gsub("_HUMAN", "", skinatlas$Entry.name)
      skinatlas[,1:13] <- as.data.frame(lapply(skinatlas[,1:13], as.numeric))

### Functions
#From:Fantom5 database
MelanocyteSelective <- function(num=2){
  if(!exists("MeanCluster_Fantom_t")) {
    load("./Data/MembraneDataAll.RData")
  }
  df <- Fantom
  names(df)
  MelanocyteIndex <- grep("Melan", names(df))
    MelanocyteMax <- apply(df[,MelanocyteIndex], 1, max)
    OthersMac   <- apply(df[,c(-1, -MelanocyteIndex)], 1, max)
    OthersMac   <- apply(df[,c(-1, -MelanocyteIndex)], 1, function(x)(sort(x,TRUE)[num]))
    MelaSelFC <- ((MelanocyteMax+1)/(OthersMac+1))
    DF <- data.frame(df$X, MelanocyteMax, OthersMac,  MelaSelFC)
      DF <- DF[order(DF$MelaSelFC, decreasing = T),]
  return(DF)
}

    DF[grep("DCT$", DF$df.X),]
    DF[grep("TYRP1$", DF$df.X),]
    DF[grep("TYR$", DF$df.X),]
    DF[grep("B4DNE6$", DF$df.X),]
    
    DF1 <- DF[1:100,]
      unique(gsub("^p[0-9]+@", "", DF1$df.X))
  
###From GTEX
SkinSelective <- function(){
  if(!exists("Gtex")) {
    load("./Data/MembraneDataAll.RData")
  }
 
  plotData <-data.frame(Gtex[,3:55])
    dfSkin <- plotData
    names(dfSkin)
    SkinIndex <- grep("Skin", names(dfSkin))
    SkinMax <- apply(dfSkin[,SkinIndex], 1, max)
    OthersMax   <- apply(dfSkin[,c(-1, -SkinIndex)], 1, max)
    SkinSelFC <- ((SkinMax+1)/(OthersMax+1))
    DFSkin <- data.frame(Gtex$Description, SkinMax, OthersMax,  SkinSelFC)
    DFSkin <- DFSkin[order(DFSkin$SkinSelFC, decreasing = T),]
   return(DFSkin)
}
#Merge(DFSKIN and MelanocyteSelective )

MelaSkinMerge <- function(){
  Temp1 <- DFMela
    Temp1$genename<- gsub("p[0-9]+@","", Temp1$df.X)
    Temp2 <- merge(Temp1, DFSkin, by.x="genename", by.y="Gtex.Description")
      Temp2 <- Temp2[order(Temp2$MelaSelFC, decreasing=T),]
      return(Temp2)
}

SelectMelocyteSkinSpecificGenes <- function(Melfc=2, SkinFC=1){
  MelaSkinSpecific <- MelaSkinMerge()
  MelaSkinSpecific <- MelaSkinSpecific[c(MelaSkinSpecific$MelaSelFC>Melfc & MelaSkinSpecific$SkinSelFC>SkinFC),]

  require(data.table) ## 1.9.2
  group <- as.data.table(MelaSkinSpecific)
  #If you want to keep all the entries corresponding to max values of pt within each group:
  
  group <- unique(group[group[, .I[MelaSelFC == max(MelaSelFC)], by=genename]$V1])
return(group)
}

DFSkin <- SkinSelective() #Gtex
DFMela <- MelanocyteSelective(num=1) #Fantom%
SelectMelocyteSkinSpecificGenes(Melfc=2, SkinFC=1.5)[1:100,1]




uniprot_mapping <- function(ids) {
  uri <- 'http://www.uniprot.org/uniprot/?query='
  idStr <- paste(ids, collapse="+or+")
  format <- '&format=tab'
  fullUri <- paste0(uri,idStr,format)
  dat <- read.delim(fullUri)
  dat
}


      
SkinAtlasExpress <- function(gene="TYRP1"){
  index<- which(skinatlas$Entry.name%in%gene)
  index2<- which(skinatlas$Gene.names%in%gene)
  index3<- unlist(lapply(gene, function(x) grep(x, skinatlas$Gene.names)))
    output <- skinatlas[c(index,index2, index3),]

    return(output)
}      

genelist <- SelectMelocyteSkinSpecificGenes(Melfc=2, SkinFC=0)[1:100,1]
O1 <- SkinAtlasExpress(gene=as.data.frame(genelist)[,1])
O1[order(O1$Melanocyte, decreasing = T),]

O1[unlist(lapply(lapply(O1$Gene.names, function(x)(strsplit(as.character(x), " "))), function(y)y%in%as.data.frame(genelist)[,1])),]


data.frame(lapply(O1$Gene.names, function(x)(strsplit(as.character(x), " "))))
       
library(tidyr)
O1 <- O1 %>%
  separate(Gene.names, c("one", "two", "three", "four", "five", "six")," ", convert = T)

INDEX <- NULL
for (i in 1:nrow(genelist)){
  for(j in 16:21){
    INDEX <- c(INDEX, grep(paste0("^", as.data.frame(genelist)[i,1], "$"), O1[,j]))
  }
}  

O1 <- unique(O1[INDEX,])
O1[order(O1$Melanocyte, decreasing = T),]





####PLot of expression

SkinPlot <- function(gene="MLANA"){
  require(ggplot2)
  require(cowplot)
  if(!exists("MeanCluster_Fantom_t")) {
    load("./Data/MembraneDataAll.RData")
  }
  df <- Fantom
  plotdata1 <- df[grep(gene, df$X),]
  MelanocyteIndex <- grep("Melan", names( plotdata1))
  MelanocyteMax <- apply( plotdata1[,MelanocyteIndex], 1, max)
  index <- which(MelanocyteMax==max(MelanocyteMax))
  plotdata1 <- data.frame(t(plotdata1[index,]))
    colnames(plotdata1) <- "Gene"
    plotdata1 <- data.frame(plotdata1[-1,])
    colnames(plotdata1) <- "Gene"
    plotdata1$Gene <- as.numeric(as.character(plotdata1$Gene))
    plotdata1$Tissue <- row.names(plotdata1)
  p <- ggplot(plotdata1, aes(y=Gene, x=Tissue))+geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  #Gtex
  if(!exists("Gtex")) {
    load("./Data/MembraneDataAll.RData")
  }
  
  PlotDataGtex <-data.frame(Gtex[,c(1,3:55)])
  PlotDataGtex <- PlotDataGtex[PlotDataGtex$Description==gene,]
  PlotDataGtex <- PlotDataGtex[-1]
  PlotDataGtex <- data.frame(t(PlotDataGtex))
  PlotDataGtex$Tissue <- row.names(PlotDataGtex)
  colnames(PlotDataGtex)[1]<- "Gene"
  q <- ggplot(  PlotDataGtex, aes(y=Gene, x=Tissue))+geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(gene)
  q
  
  SkinAtlasPlotData <- skinatlas[grep(gene, skinatlas$Gene.names),]
    if(nrow(SkinAtlasPlotData)==0){SkinAtlasPlotData = skinatlas[1,]
                                   SkinAtlasPlotData[1:13] <- 0
                                   SkinAtlasPlotData[14:16] <- "O"
    }
  SkinAtlasPlotData1 <- SkinAtlasPlotData[1:13] 
  SkinAtlasPlotData1 <- data.frame(t(  SkinAtlasPlotData1 ))
  colnames(SkinAtlasPlotData1) <- "Exprs"
  SkinAtlasPlotData1$Tissue <- row.names(SkinAtlasPlotData1)
  pp <- ggplot(  SkinAtlasPlotData1, aes(y=Exprs, x=Tissue))+geom_point() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(SkinAtlasPlotData$Gene.names)
  pp

  print(plot_grid(p,plot_grid(q,pp),nrow=2))
  #return(list(p,q,pp))
  }

Output <- SkinPlot(gene="MLANA")
plot_grid(Output[[1]],plot_grid(plotlist = c(Output[[2]],Output[[3]])),nrow=2)

SkinPlot(gene="S100A1")

MarkerList <- c("HIBCH","ASPH","TYR",
"CDH3","PLP1","KIAA1598","SOX10",
"DAB2","SORBS1",#"TTYH3","TPCN2",
"FYCO1","CHCHD6","LIMS2")

MarkerList2 <- c("MBP", "TNS3",
"LYST",#"PLEKHA5",
"TIMP2","GOLGA7B",
"SLC6A17",#"MED15",
#"MCF2L",
"SLC6A8",
#"HPS4",
"TSPAN10","CCL18","ACP5",
"PRKD3",#"RGS12",
"PDE3A","TNS1")

MarkerList <- SelectMelocyteSkinSpecificGenes(Melfc=10, SkinFC=0)[,1]
c(MarkerList)

pdf("MelanocyteMarkers3.pdf", width=15, height = 15)
  for (i in c(MarkerList$genename)){print(i)
    SkinPlot(gene=i)
  }

dev.off()


pdf("MelanocyteMarkers2.pdf", width=15, height = 15)
for (i in MarkerList2){
  SkinPlot(gene=i)
}
dev.off()

SkinPlot(gene="GMPR")
#######Selection of Melanocyte Markers


#Step1: 
Markers1 <- SkinSelective()
  Markers1 <- dim(Markers1[Markers1$SkinSelFC>2,]) #302
  
#Step2:
  DFSkin <- SkinSelective() #Gtex
  DFMela <- MelanocyteSelective(num=1) #Fantom%
  SelectMelocyteSkinSpecificGenes(Melfc=2, SkinFC=0)[1:100,1]
  
  SelectMelocyteSkinSpecificGenes
  

