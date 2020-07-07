

library(ggplot2)
library(DBI)
library(Biobase)
library(ggpubr)
library(cowplot)

textsize=16
#load("./Data/MembraneDataAll.RData")

#' Get values from shinyTree object
#'
#' @param tree_input_value Input values selected from shinyTree objects
#' @param crop_selected the name of the crop
#' @description This function gets all the selected values from shinyTree checkboxes
#' @export
#'

#if(is.null(input$tree)) {print("omar")}
# get_selected(tree)
get_tree_value <- function(tree_input_value,crop_selected){
  
  trait_selected <- unlist(get_selected(tree_input_value))
  trait_selected <- stringr::str_replace_all(string = trait_selected,pattern = ":.*" ,replacement = "")
  
  #trial_headers <- fbmodule::list_modules(crop=crop_selected)
  
  if(crop_selected == "potato")     {tbl <- table_module_potato     } #using internal RDA
  if(crop_selected == "sweetpotato"){tbl <- table_module_sweetpotato} #using internal RDA
  
  #mdl <- tbl[tbl$crop_selected == crop_selected, c("module", "module_name")]
  mdl <- tbl[tbl$crop_selected == crop_selected, c("TRIAL_ABBR", "TRIAL")]
  mdl <- paste0(mdl[,2], " (", mdl[, 1],")")
  trial_headers <- sort(unique(mdl))
  
  trial_headers <- str_trim(gsub("\\(.*","", trial_headers ),side = "both")
  trait_selected <- trait_selected[!is.element(el = trait_selected, set = trial_headers)]
  #print(trait_selected)
}


PlotGTEXFunction <- function(gene="HAPLN3") {  
  if(!exists("Gtex")) {
    load("./Data/MembraneDataAll.RData")
  }
  index <- grep(paste0("^", gene, "$"), Gtex$Description) 
  plotData <-data.frame(t(data.frame(Gtex[index[1],3:55])))
  plotData$Tissue <- row.names(plotData)
  colnames(plotData) <- c("Value", "Tissue")
  
  a <- ggplot(plotData[,], aes(x=Tissue, y=as.numeric(Value), fill = Tissue))+geom_col()+
    theme(legend.position="none")+ ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
  
  b <- ggplot(plotData[-52,],aes(x=reorder(Tissue,-as.numeric(Value)),y=as.numeric(Value),label = Tissue,fill = Tissue))+
    geom_bar(stat = 'identity')+
    labs(x='', y = 'Expression') + ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  return(list(a,b))      
}      

PlotFantonFunction <- function(gene="HAPLN3") { 
  if(!exists("MeanCluster_Fantom_t")) {
    load("./Data/MembraneDataAll.RData")
  }
  df <- data.frame()
  a <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
  b <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
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
     
      a <- ggplot(plotData[-52,], aes(x=Cell, y=as.numeric(Value), fill = Cell))+geom_col()+
        theme(legend.position="none")+ ggtitle(paste0("Fantom5 Expression Data ", gene))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
      
      b <- ggplot(plotData[-52,],aes(x=reorder(Cell,-as.numeric(Value)),y=as.numeric(Value),label = Cell,fill = Cell))+
        geom_bar(stat = 'identity')+
        labs(x='', y = 'Expression', title=paste0("Fantom5 Expression Data ", gene))+
        theme_bw()+
        theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  }#end if
      
       return(list(a,b))      
}

#Load Data.  run only once when the server.R is called
ImmNavPlot <-function( Symbol = "IL6"){
  ptm <- proc.time()
  con <- dbConnect(RSQLite::SQLite(), "~/projects/ImmuneNavigator/Data/immunnav")
  #1 Identify probes for genes of interest
  PID <- dbGetQuery(con, paste0("SELECT * FROM GeneAnnotation WHERE SYMBOL = '", Symbol,"'"))
  #2 Query Expression Data for individual tables for PID 
  Query <- NULL
  for (i in c("immnav1","immnav2","immnav3","immnav4","immnav5")){print (i)
    # Pass multiple sets of values with dbBind():
    rs <- dbSendQuery(con, paste0("SELECT * FROM ", i,  " WHERE row_names =?"))
    rs <- dbBind(rs, list(PID$PROBEID))
    Query <- c(Query, dbFetch(rs))
  }
  
  Query <- as.data.frame(Query)
  Query_t=as.data.frame(t(Query[,-1]))
  names(Query_t) <- Query[,1]
  Query <- Query_t
    IndexRow <- grep("row", rownames(Query),invert = T)
    Query <-Query[IndexRow,]
  #3 QueryMeta data
  Cells <- dbGetQuery(con, paste0("SELECT * FROM Meta"))
  Cells$Sample <- make.names(Cells$Sample)
  Query_m <- merge(Query, Cells, by.x=0, by.y="Sample", all=T)
  names(Query_m)[2:(ncol(Query_m)-2)] <- "Symbol"
  index <- which (names(Query_m)=="Symbol")
  ifelse (length(index)==1, Query_m[,index] <- Query_m[,index],
          Query_m[,index]<- apply(Query_m[,index],2,as.numeric))
 
  
   #4 Plot Data
  #Determine highest expressed Row
  #Conver factors to numeric
  #indx <- sapply(Query_m, is.factor)
  #Query_m[indx] <- lapply(Query_m[indx], function(x) as.numeric(as.character(x)))
  
  MaxExp <- apply(Query_m[,index],2,max) 
  index3 <- which(MaxExp==max(MaxExp))
  index4 <-(2:(ncol(Query_m)-2))[index3]
  index5 <- which (colnames(Query_m)%in% c("Row.names","row_names","Cell"))
  Query_m <- Query_m[,c(index5, index4)]
    colnames(Query_m)[4] <- "value"
    Query_m$value <- as.numeric(as.character((Query_m$value)))
  p <- ggplot(Query_m, aes(x=Cell, y=value))+ geom_boxplot(outlier.shape = NA)+
    geom_jitter(width=0.2)+ggtitle(Symbol) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = textsize))
  
  dbDisconnect(con) 
  print(proc.time() - ptm)
  #     #Merge with GSE
  #Query_mm <-  merge(Query_m, GSE, by.x="Row.names", by.y="names", all.x=TRUE)
  # q <- ggplot(Query_mm, aes(x=Cell, y=log2(as.numeric(Symbol))))+ geom_boxplot(outlier.shape = NA)+
  #   geom_jitter(width=0.2)+ggtitle(Symbol) + facet_wrap(~series_id)+
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(p)
}

###DICE2 Data set
plotDataFunction <- function(Gene=Gene, index=1){
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
  
  #Index <- which(M2$Additional_annotations==Gene)
  #print(M2$Additional_annotations[Index])
  #Index <- Gene
  #print(M2$Additional_annotations[Index])
  print(sprintf("Index %s", Index))
  for (j in DICE_Cell_Names){
    k=k+1
    Data[1:length(get(j)[Index,]),1]<- j
    Data[1:length(get(j)[Index,]),2]<- t(get(j)[Index,])
    Data <- Data[-1:-3,1:2 ]
    Data1<- data.frame(rbind(Data,Data1))
  }
  
  Data1 <- data.frame(as.matrix(Data1)[,1:2])
  colnames(Data1) <- c("Cell", "Value")
  
  if(Index==1 & Gene != "TSPAN6"){Data1$Value=rep(0,nrow(Data1))}
 
  
  p <- ggplot(Data1, aes(x=Cell, y=as.numeric(as.character(Value))))+geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2, height = 0)
  p <- p + ggtitle(Gene)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size = textsize))
  p <- p + ylab("Expression Value")
  return(p)
}  

###Immune Cell Proteomics
plotImmuneProteomics <- function(gene="HAPLN3") {  
  if(!exists("ImmuneProtgset")) {
    load("~/CSI/projects/proteomics.dev/Data/Proteomics.RData")
    ImmuneProtgset <- gset
    rm(gset)  
  }

index <- grep( gene, fData(ImmuneProtgset)$Gene.names.1) 
plotData <-(data.frame(exprs(ImmuneProtgset)[index[1],]))
plotData$Cell <- pData(ImmuneProtgset)$Cell
plotData$State <- pData(ImmuneProtgset)$State
colnames(plotData) <- c("Value", "Cell", "State")

a <- ggplot(plotData[,], aes(x=Cell, y=as.numeric(Value),fill = factor(State)))+geom_boxplot()+
   ggtitle(label=paste0("Immune Cell Protein Expression Data ", gene)) +
  geom_point(position=position_dodge(width=.75))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
return(a)      
}      

#AMP Low Input
PlotLowInputAMP <- function(gene="CDH11"){
  if(!exists("lowinputGSET")) {
    load("~/projects/CDH11_membrance/CDH11/Data/ProteinClass/Data/LowInput.RData")
  }
  
  PlotData <- exprs(lowinputGSET)[fData(lowinputGSET)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(lowinputGSET))
  a <-  ggplot(PlotData_t,aes(x=DiseaseTissue,y=Value,label = DiseaseTissue,fill = DiseaseTissue)) + 
    geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_jitter(width=0.2)+
    labs(x='', y = 'Expression', title=paste0("AMP Low Input Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  
  return(a)
} 

#AMP Low Input Phase II
PlotLowInputAMPII <- function(gene="CDH11"){
  if(!exists("lowinputAMPPhaseII")) {
    load("~/Projects/AMP_Phase2_RA/Data/AMPPhase2.RData")
    lowinputAMPPhaseII <- gset
    rm(gset)  
  }
  
  PlotData <- exprs(lowinputAMPPhaseII)[fData(lowinputAMPPhaseII)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(lowinputAMPPhaseII))
  a <-  ggplot(PlotData_t,aes(x=interaction(Group, Diagnosis),y=Value,label = Group,fill = Group)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_jitter(width=0.2)+
    labs(x='', y = 'Expression')+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a
  PlotData_t$Cohort <- gsub(" APPROVED", "", PlotData_t$Cohort)
  #PlotData_t$Cohort <- gsub("", "OA", PlotData_t$Cohort)
   Index <- grep("Group", PlotData_t$Cohort, invert = T)
   PlotData_t$Cohort[Index] <- "OA"
  b <-  ggplot(PlotData_t,aes(x=Group,y=Value,fill = Cohort)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_point(aes(fill=Cohort),shape=21,size=3,color="black",stroke=1,position = position_jitterdodge(jitter.width=0.1))+
    labs(x='', y = 'Expression')+
    theme_bw()+
   # theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  b

  # now add the title
  title <- ggdraw() + 
    draw_label(
      paste0("AMP Phase II Low Input Expression Data ",gene),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  Output <- plot_grid(
    title, plot_grid(a,b, align = 'h', rel_widths =   c(.7, 1)),
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  

  return(Output)
} 

#AMP Low Input
PlotAMP <- function(gene="CDH11"){
  if(!exists("LeucoRichGSET")) {
    load("~/projects/AMPViewer/amp_phase1_ra_viewer/data/LeucoRichAmp.Rdata")
    LeucoRichGSET <- gset
    fData(LeucoRichGSET)$gene_name <- fData(LeucoRichGSET)$gene.name
  }
  gene <- as.character(gene)
  #print(gene)
  PlotData <- exprs(LeucoRichGSET)[fData(LeucoRichGSET)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(LeucoRichGSET))
    PlotData_t <- PlotData_t[,c(1,5,91)]
    PlotData_t$Group1 <- factor(PlotData_t$Group1)
    print(PlotData_t[1:10,])
  my_comparisons <- list( c("leukocyte.poor.RA" , "leukocyte.rich.RA"))
  maxval <- max(PlotData_t$Value)+1
  a <-  ggplot(PlotData_t,aes(x=Group1,y=Value,fill = Group1)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2)+ facet_wrap(~cell_type, nrow=1)+
    labs(x='', y = 'Expression', title=paste0("AMP Low Input Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  return(a)
} 

#AMP Low Input
PlotAMPSingleCell <- function(gene="CDH11"){
  if(!exists("LeucoRichGSET")) {
    load("~/projects/CDH11_membrance/CDH11/Data/AmpSingleCell.RData")
  }
  PlotData <- as.data.frame(t(MeanCluster_t[row.names(MeanCluster_t)==gene,]))
    PlotData$Cell <- row.names(PlotData)
    names(PlotData)[1] <- "Value"
    PlotData$Name <- c("NaiveB", "MemB", "ABC", "PlasmaBlasts", "CD34+Sub", "HLA+Sub", "DKK3+Sub", "CD55+Lining",
                       "IL1B+ pro-infl", "NUPR1+", "C1QA+", "IFN-activated", "CCR7+ CD4+", "FOXP3+ Tregs",
                       "PD-1+ Tph/Tfh", "GZMK+ CD8+", "GNLY+ GZMB+ CTLs", "GZMK+/GZMB+")

  a <-ggplot(PlotData,aes(x=Cell,y=Value,fill = Cell, label=Name)) +
    geom_bar(stat = "identity")+
    labs(x='', y = 'Expression', title=paste0("AMP Single Cell Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
    # Add global p-value
  #Add text annotation
  a <- a + geom_text(y=0, angle = 90, vjust = 1, hjust = 0)
  return(a)
} 

#Brenner Fibriblast
PlotBrennerFibro <- function(gene="CDH11"){
  if(!exists("BrennerFibroGSET")) {
    load("~/projects/CDH11_membrance/CDH11/Data/BrennerFib.RData")
    BrennerFibroGSET <- gset_RNAseq
    fData(BrennerFibroGSET)$gene_name <- fData(BrennerFibroGSET)$external_gene_name
  }
  gene <- as.character(gene)
  
  PlotData <- exprs(BrennerFibroGSET)[fData(BrennerFibroGSET)$gene_name==gene,]
  PlotData <- PlotData [!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(BrennerFibroGSET))
  
  PlotData_t$group <- factor(PlotData_t$group)
  # Add groups for FAB+PDPN+ thy1 +/-
  PlotData1 <- exprs(BrennerFibroGSET)[fData(BrennerFibroGSET)$gene_name %in% c("PDPN", "FAP", "THY1"),]
  
  PlotData_t1 <-data.frame(t(PlotData1))
  
  colnames(PlotData_t1) <- fData(BrennerFibroGSET)[fData(BrennerFibroGSET)$gene_name %in% c("PDPN", "FAP", "THY1"),2]
  
  PlotData_t <- data.frame(PlotData_t, PlotData_t1)
  PlotData_t$FABGroup [PlotData_t$FAP>150 & PlotData_t$PDPN>200]<- as.character(
    PlotData_t[PlotData_t$FAP>150 & PlotData_t$PDPN>200,"group"])
  
  PlotData_t$FABGroup <- gsub("CD34[HL]_THY1L_CDH11[HL]", "FABH_PDPNH_THY1L",PlotData_t$FABGroup )
  PlotData_t$FABGroup <- gsub("CD34[HL]_THY1H_CDH11[HL]", "FABH_PDPNH_THY1H",PlotData_t$FABGroup )
  
  my_comparisons <- list( c("FABH_PDPNH_THY1H" , "FABH_PDPNH_THY1L"))
  maxval <- max(PlotData_t$Value)+1
  a <-  ggplot(PlotData_t,aes(x=group,y=Value,fill = group)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter(width=0.2)+
    labs(x='', y = 'Expression', title=paste0("Brenner Fibroblasts ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a <- a + geom_boxplot(aes(x=FABGroup, fill=FABGroup), outlier.shape = NA)+ 
    geom_jitter(aes(x=FABGroup, fill=FABGroup),width=0.2)+
    geom_vline(xintercept = 7.5)
  return(a)
} 

ggplot(PlotData_t, aes(x=FABGroup, y=donor))+geom_jitter(height=0.2, width=0.2)

#Brenner Fibriblast-Micro
PlotBrennerFibroMicro <- function(gene="CDH11"){
  if(!exists("microgset")) {
    load("~/projects/CDH11_membrance/CDH11/Data/fibroblastmicrosubset.RData")
    fData(microgset)$gene_name <- fData(microgset)$y
  }
  gene <- as.character(gene)
  
  PlotData <- exprs(microgset)[fData(microgset)$gene_name==gene,]
 # PlotData <- PlotData [!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(microgset))
  
  PlotData_t$group <- factor(PlotData_t$Group)
  # Add groups for FAB+PDPN+ thy1 +/-
  PlotData1 <- exprs(microgset)[fData(microgset)$gene_name %in% c("PDPN", "FAP", "THY1"),]
  
  PlotData_t1 <-data.frame(t(PlotData1))
  
  colnames(PlotData_t1) <- fData(microgset)[fData(microgset)$gene_name %in% c("PDPN", "FAP", "THY1"),"gene_name"]
  
  PlotData_t <- data.frame(PlotData_t, PlotData_t1)

  maxval <- max(PlotData_t$Value)+1
  a <-  ggplot(PlotData_t,aes(x=group,y=Value,fill = disease.ch1, shape=disease.ch1)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter( position = position_jitterdodge(), size=3)+
    labs(x='', y = 'Expression', title=paste0("Brenner Fibroblast Micro ",gene))+
    theme_bw()+
    guides(fill = FALSE)+
    theme(legend.position=c(0,1))+
    #guides(shape = guide_legend(override.aes = list(size = 1)))+
    theme(legend.key.size = unit(2,"line"))+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a <- a + 
    guides(shape = FALSE)+
    theme(legend.position=c(0,1))+guides(fill = guide_legend(override.aes = list(size = 0)))
  return(a)
} 

#Brenner Single cell Fibriblast
PlotSSRA <- function(gene="CDH11"){
  if(!exists("SSRA")) {
    load("~/projects/CDH11_membrance/CDH11/Data/BrennerFib.RData")
    
    fData(SSRA)$gene_name <- fData(SSRA)$gene.name
  }
  gene <- as.character(gene)
  
  PlotData <- exprs(SSRA)[which(fData(SSRA)$gene_name==gene),]
  PlotData_t <-data.frame(Value=PlotData, pData(SSRA))
    PlotData_t$CD34 <- PlotData_t$CD34_protein
    PlotData_t$THY1 <- PlotData_t$THY1_protein
      PlotData_t$CD11 <- PlotData_t$CDH11_protein
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
  a <-  ggplot(PlotData_t,aes(x=group,y=Value,fill = disease, shape=disease)) +
    geom_boxplot(outlier.shape = NA)+geom_jitter( position = position_jitterdodge(), size=3)+
    labs(x='', y = 'Expression', title=paste0("Single Cell Fibroblast Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a <- a + geom_boxplot(aes(x=FABGroup, y=Value), outlier.shape = NA)+ 
    geom_jitter( aes(x=FABGroup, y=Value), position = position_jitterdodge(), size=3)+
    geom_vline(xintercept = 7.5)
  a <- a + 
    guides(shape = FALSE)+
    theme(legend.position=c(0,1))+guides(fill = guide_legend(override.aes = list(size = 0)))
  return(a)
} 

#Fibroblast signature
PlotFIBData <- function(gene="CDH11"){
  if(!exists("fibrgset")) {
    load("~/app/RAExpressionDashboard/Data/Fibroblastgset.RData")
    fData(fibrgset)$gene_name <- fData(fibrgset)$SYMBOL
  }
  gene <- as.character(gene)
  
  PlotData1 <- exprs(fibrgset)[fData(fibrgset)$gene_name == gene,]
  PlotData1 <-   PlotData1 [!is.na(row.names(PlotData1)),]
  Index=1
  if(!is.null(nrow(PlotData1))) {
    MAX <- apply(PlotData1, 1, max)
    Index <- which(MAX==max(MAX))
    PlotData1 <- PlotData1 [Index,]
    }
  

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
    labs(x='', y = 'Expression', title=paste0("Fibroblast Cell Expression Data ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  
  return(a)
} 

#PEAK Study
PlotPeak <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("PEAKgset")) {
    load("~/projects/PEAK/RNASeqQC-SJ/Data/PEAK.RData")
    PEAKgset <- gset
    rm(gset, gsetMax, cont.wt1)
  }
  PlotData <- exprs(PEAKgset)[fData(PEAKgset)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(PEAKgset))
  
  my_comparisons <- list( c("fibroid" , "lymphoid")
  )
  
  a <-  ggplot(PlotData_t,aes(x=organism_part,y=Value,label = pathotype,fill = pathotype)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0,
                         dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("PEAK Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  maxval <- max(PlotData_t$Value)
  a <- a + #stat_compare_means(aes(group = pathotype))+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
   return(a)

}

##Publci Synovail
PlotSynovial <- function(gene="CDH11"){
  require(ggplot2)
  require(ggpubr)
  if(!exists("Publicgset")) {
    load("~/projects/SynoviumPublic/RNASeqQC-SJ/Data/PublicSynovial.RData")
    Publicgset <- gset
    rm(gset, gsetMax, cont.wt1)
  }
  PlotData <- exprs(Publicgset)[fData(Publicgset)$gene_name==gene,]
 # if(dim(PlotData)[1]==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  PlotData_t <-data.frame(Value=PlotData, pData(Publicgset))
  PlotData_t$Group <- factor(PlotData_t$Group, levels=c( "Normal_Base", "OA", "Arthralgia_Base",  "UnDiffArth", "RA_Early_Base", "RA_Early_6Mont",  "RA_Est"))
 
  my_comparisons <- list( c("Normal_Base", "OA"), 
                          c("Normal_Base", "Arthralgia_Base"),
                          c("Normal_Base", "UnDiffArth"),
                          c("Normal_Base", "RA_Early_Base"),
                          c("Normal_Base", "RA_Early_6Mont"),
                          c("Normal_Base", "RA_Est")
                          )
  maxval <- max(PlotData_t$Value)
   a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = disease,fill = disease)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.15, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("Public Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
 
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+5)     # Add global p-value
  
   return(a)
  
}

#Recombine
PlotRecombine <- function(gene="CDH11"){
    require(ggplot2)
    if(!exists("Recombinegset")) {
      load("~/projects/COMBINE_Synovial/Synovial/RNASeqQC-SJ/Data/Synovial.RData")
      Recombinegset <- gset
      rm(gset, gsetMax, cont.wt1)
    }
    PlotData <- exprs(Recombinegset)[fData(Recombinegset)$hgnc_symbol==gene,]
    PlotData <- PlotData[!is.na(row.names(PlotData)),]
    PlotData_t <-data.frame(Value=PlotData, pData(Recombinegset))
    PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
    PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
                                                          "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
    my_comparisons <- list( c("MTXBiopsy01", "MTXBiopsy02"), 
                            c("iTNFBiopsy01", "iTNFBiopsy02"),
                            c ("BioBiopsy01" , "BioBiopsy02"))
    maxval <- max(PlotData_t$Value)
    
    a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Treatment.name,fill = Treatment.name)) +
      geom_boxplot(outlier.shape = NA)+
      geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                               dodge.width = 0.75, seed = NA))+
      #geom_point(position=position_dodge(width=.75))+
      labs(x='', y = 'Expression', title=paste0("Synovial Expression Data ",gene))+
      theme_bw()+
      #theme(legend.position="none")+
      theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
    a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
      stat_compare_means(label.y = maxval+1)     # Add global p-value
    
    return(a)
    
}

#GSE55235
PlotGSE55235 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE55235gset")) {
    GSE55235gset <- readRDS("~/Cytoreason/GSE55235.rds")
  }
  PlotData <- data.frame(Value=exprs(GSE55235gset )[fData(GSE55235gset)$SYMBOL==gene,])
    if(dim(PlotData)[1]==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  PlotData_t <-data.frame(Value=PlotData, pData(GSE55235gset))
  PlotData_t$Group = PlotData_t$condition
                                                    
  my_comparisons <- list( c("healthy", "osteoarthritis"), 
                          c("healthy", "rheumatoid arthritis"),
                          c ("osteoarthritis" , "rheumatoid arthritis"))
  maxval <- max(PlotData_t$Value)
  
   a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
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

#GSE55235
PlotGSE55457 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE55457gset")) {
    GSE55457gset <- readRDS("~/Cytoreason/GSE55457.rds")
  }
  PlotData <- exprs(GSE55457gset)[fData(GSE55457gset)$SYMBOL==gene,]
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(GSE55457gset ))
  PlotData_t$Group = PlotData_t$condition
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
  #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
  my_comparisons <- list( c("Healthy/Control", "osteoarthritis"), 
                          c("Healthy/Control", "rheumatoid arthritis"),
                          c ("osteoarthritis" , "rheumatoid arthritis"))
  maxval <- max(PlotData_t$Value)
  
  a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE55457 Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
  return(a)
  
}

#GSE1919
PlotGSE1919 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE1919gset")) {
    GSE1919gset <- readRDS("~/Cytoreason/GSE1919.rds")
  }
  PlotData <- exprs(GSE1919gset)[fData(GSE1919gset)$SYMBOL==gene,]
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(GSE1919gset))
  PlotData_t$Group = PlotData_t$condition
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
  #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
  
  my_comparisons <- list( c("Healthy/Control", "osteoarthritis"), 
                          c("Healthy/Control", "rheumatoid arthritis"),
                          c ("osteoarthritis" , "rheumatoid arthritis"))
  maxval <- max(PlotData_t$Value)
  
   a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE1919gset Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())

  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
  return(a)
  
}

#GSE77298
PlotGSE7729 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE7729gset")) {
    GSE7729gset <- readRDS("~/Cytoreason/GSE7729.rds")
  }
  PlotData <- exprs(GSE7729gset)[fData(GSE7729gset)$SYMBOL==gene,]
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(GSE7729gset))
  PlotData_t$Group = PlotData_t$condition
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
  #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
  
  my_comparisons <- list( c("Healthy/Control", "rheumatoid arthritis") )
  maxval <- max(PlotData_t$Value)
  a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE7729 Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  # Visualize: Specify the comparisons you want
    a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
return(a)
  
}

#SDY1299
PlotSDY1299 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("SDY1299gset")) {
    load("~/projects/SynovialImmport/RNASeqQC-SJ/SynovialImmport.RData")
    SDY1299gset <- gset
    rm(gset,gsetMax)  
  }
  df <- data.frame()
  a <-  ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100)
  
  if (gene %in% fData(SDY1299gset)$gene_name){
      PlotData <- exprs(SDY1299gset)[fData(SDY1299gset)$gene_name==gene,]
        
      #PlotData <- PlotData[!is.na(row.names(PlotData)),]
      PlotData_t <-data.frame(Value=PlotData, pData(SDY1299gset))
      #PlotData_t$Group = PlotData_t$condition
      #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
      #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
      #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
      
      my_comparisons <- list( c("OA", "RA") )
      maxval <- max(PlotData_t$Value)
      a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
        geom_boxplot(outlier.shape = NA)+
        geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                                 dodge.width = 0.75, seed = NA))+
        #geom_point(position=position_dodge(width=.75))+
        labs(x='', y = 'Expression', title=paste0("SDY1299 Synovial Expression Data ",gene))+
        theme_bw()+
        #theme(legend.position="none")+
        theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
      # Visualize: Specify the comparisons you want
      a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
        stat_compare_means(label.y = maxval+1)     # Add global p-value
  }#end if  
      
  return(a)
  
}

#GSE24742
PlotGSE24742 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE24742gset")) {
    GSE24742gset <- readRDS("Data/GSE24742.rds")
  }
  # PlotData <- exprs(GSE24742gset)[fData(GSE24742gset)$SYMBOL==gene,]
  # 
  # PlotData <- PlotData[!is.na(row.names(PlotData)),]
  # if(nrow(PlotData)==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  
  
  PlotData <- data.frame(Value=exprs(GSE24742gset )[fData(GSE24742gset)$SYMBOL==gene,])
  PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData <- PlotData[!is.na(PlotData[,1]),]
  PlotData <- PlotData[is.numeric(PlotData[,1]),]
  if(dim(PlotData)[1]==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=t(PlotData), pData(GSE24742gset))
  names(PlotData_t)[1] <- "Value" 
  PlotData_t$Group = PlotData_t$Treatment
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
  #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
  
  my_comparisons <- list( c("baseline", "12 weeks of RTX therapy"))
  maxval <- max(PlotData_t$Value)
  a <- ggplot(PlotData_t,aes(x=Response,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE24742 Synovial Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  # Visualize: Specify the comparisons you want
  a <- a + stat_compare_means(aes(group=Response), label = "p.format")+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
  return(a)
  
}

#GSE45867
PlotGSE45867 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE45867gset")) {
    GSE45867gset <- readRDS("Data/GSE45867.rds")
  }
  # PlotData <- exprs(GSE45867gset)[fData(GSE45867gset)$SYMBOL==gene,]
  # PlotData <- PlotData[!is.na(row.names(PlotData)),]
  #   if(nrow(PlotData)==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
  
  PlotData <- data.frame(Value=exprs( GSE45867gset )[fData( GSE45867gset)$SYMBOL==gene,])
    PlotData <- PlotData[!is.na(row.names(PlotData)),]
    PlotData <- PlotData[!is.na(PlotData[,1]),]
    #PlotData <- PlotData[is.numeric(PlotData[,1]),]
    if(dim(PlotData)[1]==0){PlotData=data.frame(Value=rep(0,ncol(PlotData)))}
    
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=t(PlotData), pData(GSE45867gset))
  names(PlotData_t)[1] <- "Value"
  PlotData_t$Group = PlotData_t$time
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("MTXBiopsy01",  "MTXBiopsy02", "iTNFBiopsy01", "iTNFBiopsy02", 
  #                                                      "BioBiopsy01", "BioBiopsy02", "BioBiopsy03" ))
  
  my_comparisons <- list( c("W0", "W12"))
  maxval <- max(PlotData_t$Value)
  a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE45867 Synovial + TCZ ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  # Visualize: Specify the comparisons you want
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+1)     # Add global p-value
  
  return(a)
  
}


#GSE36700
PlotGSE36700 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE36700gset")) {
    GSE36700gset <- readRDS("Data/GSE36700.rds")
  }
  gene <- as.character(gene)
  fData(GSE36700gset)$SYMBOL <- as.character(fData(GSE36700gset)$SYMBOL)
  PlotData <- exprs(GSE36700gset)[fData(GSE36700gset)$SYMBOL==gene,]
    if(class(PlotData)!="numeric"){
      PlotData <- PlotData[!is.na(row.names(PlotData)),]
      maxTable <- apply(PlotData, 1, max)
      maxTableIndex <- which(maxTable==max(maxTable))
      PlotData <- PlotData[maxTableIndex,]
    }

  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(GSE36700gset))
  PlotData_t$title <- as.character(PlotData_t$title)
  PlotData_t$Group = PlotData_t$group
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  PlotData_t$Group <- factor(PlotData_t$Group, levels=c("Rheumatoid.arthritis", "Seronegative.arthritis", "Osteoarthritis", 
                                                        "Microcrystalline.arthritis", "Systemic.lupus.erythematosus" ))
  
  my_comparisons <- list( c("Rheumatoid.arthritis", "Seronegative.arthritis"),
                          c("Rheumatoid.arthritis", "Osteoarthritis"),
                          c("Rheumatoid.arthritis", "Microcrystalline.arthritis"),
                          c("Rheumatoid.arthritis", "Systemic.lupus.erythematosus"))
  maxval <- max(PlotData_t$Value)
  a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE36700 Synovial ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  # Visualize: Specify the comparisons you want
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+2)     # Add global p-value
  a
  
  return(a)
  
}

#GSE116899
PlotGSE116899 <- function(gene="CDH11"){
  require(ggplot2)
  if(!exists("GSE116899gset")) {
    load("Data/GSE116899.rdata")
    GSE116899gset <- gset 
  }
  gene <- as.character(gene)
  fData(GSE116899gset)$SYMBOL <- fData(GSE116899gset)$Gene.name
  fData(GSE116899gset)$SYMBOL <- as.character(fData(GSE116899gset)$SYMBOL)
  PlotData <- exprs(GSE116899gset)[fData(GSE116899gset)$SYMBOL==gene,]
  PlotData <- PlotData[(!is.na(PlotData[,1])),]
  # if(class(PlotData)!="numeric"){
  #   PlotData <- PlotData[!is.na(row.names(PlotData)),]
  #   maxTable <- apply(PlotData, 1, max)
  #   maxTableIndex <- which(maxTable==max(maxTable))
  #   PlotData <- PlotData[maxTableIndex,]
  # }
  
  #PlotData <- PlotData[!is.na(row.names(PlotData)),]
  PlotData_t <-data.frame(Value=PlotData, pData(GSE116899gset))
  #PlotData_t$title <- as.character(PlotData_t$title)
  PlotData_t$Group = paste(PlotData_t$Tissue, PlotData_t$Disease, sep = "_")
  #PlotData_t <- PlotData_t[!PlotData_t$Group %in% c("Biopsy02", "OtherBiopsy01", "OtherBiopsy02", "NANA"),]
  #PlotData_t$Group <- factor(PlotData_t$Group, levels=c("Rheumatoid.arthritis", "Seronegative.arthritis", "Osteoarthritis", 
  #                                                      "Microcrystalline.arthritis", "Systemic.lupus.erythematosus" ))
  
  my_comparisons <- list( c("Blood_Healthy", "Blood_RheumatoidArthritis"),
                          c("Blood_Healthy", "SynovialFluid_RheumatoidArthritis"))
                         
  maxval <- max(PlotData_t$Value)
  a <-  ggplot(PlotData_t,aes(x=Group,y=Value,label = Group,fill = Group)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point( position=position_jitterdodge(jitter.width = 0.25, jitter.height = 0,
                                             dodge.width = 0.75, seed = NA))+
    #geom_point(aes(shape=Donor), position=position_dodge(width=.75))+
    labs(x='', y = 'Expression', title=paste0("GSE116899 Neutrophils ",gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  # Visualize: Specify the comparisons you want
  a <- a + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = maxval+2)     # Add global p-value
  a
  
  return(a)
  
}
