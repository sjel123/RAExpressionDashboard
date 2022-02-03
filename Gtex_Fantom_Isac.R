

########################   
########GTEX Data
#######################


PlotGTEXFunction <- function(gene="HAPLN3") {  
  require(ggplot2)
  #### Load Data
  #### Data is a flat file downloaded from gtexportal.org
  #### Verison GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct
  ### Added some annotation from human protein atlas (Tissue specificity columns)
  #
  if(!exists("Gtex")) {
    load("./Data/MembraneDataAll.RData")
  }
  
  #Select gene of interest
  index <- grep(paste0("^", gene, "$"), Gtex$Description)
  #Select data for plotting
  #If more than 1 row has the same gene name I select the first one
  # Should be updated to select the one with the highest expression
  plotData <-data.frame(t(data.frame(Gtex[index[1],3:55]))) # Removed data on tissue specificity (Columns 56-58)
  plotData$Tissue <- row.names(plotData)
  colnames(plotData) <- c("Value", "Tissue")
  
  a <- ggplot(plotData[,], aes(x=Tissue, y=as.numeric(Value), fill = Tissue))+
    geom_col()+
    theme(legend.position="none")+ 
    ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
  
  #Clean up labels
  plotData$Tissue <- gsub("Cells...", "", plotData$Tissue)
  plotData$Tissue <- gsub("...", "", plotData$Tissue, fixed = T)
  plotData$Tissue2 <- strtrim(plotData$Tissue, 20)
  
  #Plot data sorted from highest expression to lowest    
  b <- ggplot(plotData[-52,],aes(x=reorder(Tissue2,-as.numeric(Value)),y=as.numeric(Value),label = Tissue,fill = Tissue))+
    geom_bar(stat = 'identity')+
    labs(x='', y = 'Expression') + ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 85, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  return(list(a,b))      
} 







########################   
########Fantom5 data
#######################
#Load Data
Fantom <- read.csv("~/Shiny/Fantom5/GST_FANTOM/TPM_133.cell_median_TPM.csv")
Fantom_meta <- read.csv("~/Shiny/Fantom5/GST_FANTOM/TPM.Cells.only.samples.csv")
###Clean up Cell Labels
    TissueType     <- Fantom_meta$Sample_Type
    TissueType_sub <- gsub("\\.\\.\\.|breast\\.\\.|perirenal|omental","",TissueType)
    TissueType_sub <- gsub("subcutaneous|Bronchial|Alveolar\\.|cerebellum\\.\\.|CD14\\.2b\\.","",TissueType_sub)
    TissueType_sub <- gsub("Amniotic\\.|cerebral\\.cortex\\.\\.|CD14\\.2bCD16\\.2b\\.","",TissueType_sub)
    TissueType_sub <- gsub("CD14\\.2bCD16\\.\\.|CD14\\.CD16\\.2b\\.|CD19\\.2b\\.","",TissueType_sub)
    TissueType_sub <- gsub("\\.|Cells|Cell","",TissueType_sub)
    TissueType_sub <- gsub("CD42b|Cardiac|diff|Ciliary|Corneal|monocyteimmaturederived|Aortic|Lymphatic|Artery","",TissueType_sub)
    TissueType_sub <- gsub("CD82b|Microvascular|Thoracic|Umbilicalvein|Vein|Esophageal|Adventitial|ChoroidPlexus|Conjunctival","",TissueType_sub)
    TissueType_sub <- gsub("Mammary|PeriodontalLigament|PulmonaryVillousMesenchymal|skinnormal|Gingival","",TissueType_sub)
    TissueType_sub <- gsub("Dermal|Lung|Pulmonary|VillousMesenchymal|HairFollicle|HepaticSinusoidal|Hepatic28lipocyte29|IrisPigment","",TissueType_sub)
    TissueType_sub <- gsub("monocytederived|dark|light|derived|","",TissueType_sub)
    TissueType_sub[197:219] <- "MesenchymalStemCell"
    TissueType_sub <- gsub("Olfactory|Placental|Prostate|visceral|28polarized29|Renal","",TissueType_sub)
    TissueType_sub <- gsub("Cortical|Glomerular|ProximalTubular|SmallAirway|Bladder|","",TissueType_sub)
    TissueType_sub[319:377] <- "SmoothMuscle"
    TissueType_sub <- gsub("Tracheal","",TissueType_sub)
    TissueType_sub <- gsub("Intestinal|28lipocyte29|Vertebral|Pancreatic|RetinalPigment|breast|Lens","",TissueType_sub)
    TissueType_sub <- gsub("stromal","Stromal",TissueType_sub)
    TissueType_sub <- gsub("MesenchymalStem","MesenchymalStem Cell",TissueType_sub)
    TissueType_sub <- gsub("Smoothmuscle","SmoothMuscle",TissueType_sub) 
    TissueType_sub <- gsub("epithelial","Epithelial",TissueType_sub) 

      Fantom_meta$CellType <- TissueType_sub
      Fantom_meta$Sample_Type <- gsub("\\.\\.\\.|\\.\\.","\\.", Fantom_meta$Sample_Type)
      Fantom_meta$Sample_Type <- gsub("Smooth.Muscle.Cells.Umbilical.artery.","Smooth.Muscle.Cells.Umbilical.Artery.", Fantom_meta$Sample_Type)

      unique(Fantom_meta$Sample_Type) %in%  colnames(Fantom)
      Fantom_meta <- Fantom_meta[(!duplicated(Fantom_meta$Sample_Type)),]
      
      ####Calculated mean expression values per cell type
        MeanCluster_Fantom <- apply(Fantom[,-1],1,function(x) tapply(x, Fantom_meta$CellType, mean))
          MeanCluster_Fantom <- data.frame(MeanCluster_Fantom)
          MeanCluster_Fantom_t <- as.data.frame(t(MeanCluster_Fantom))
          MeanCluster_Fantom_t$GeneSymbol <- Fantom$X
            head(MeanCluster_Fantom_t)
            
      ####Select max value for each gene
            # loading library
            require("dplyr")
            MeanCluster_Fantom_t$Symbol <- gsub("p[0-9]+@","", MeanCluster_Fantom_t$GeneSymbol)
            MeanCluster_Fantom_t$Max <- apply(MeanCluster_Fantom_t[,1:51], 1, max)
            MeanCluster_Fantom_t_Mean <- MeanCluster_Fantom_t %>% group_by(Symbol) %>% slice(which.max(Max))

###Plot Data

PlotFantonFunction2 <- function(gene="HAPLN3") {         
  index <- grep(gene, MeanCluster_Fantom_t_Mean$Symbol) 
  plotData <- data.frame(t(MeanCluster_Fantom_t_Mean[index,]))
    plotData[] <- lapply(plotData, as.character)
    colnames(plotData) <- (plotData[53,])
      plotData$Cell <- row.names(plotData)
      colnames(plotData) <- make.names(colnames(plotData))
     # plotData$Cell <- (row.names(plotData))

       colnames(plotData)[1] <- "Value"
  
  a <- ggplot(plotData[1:51,], aes(x=Cell, y=as.numeric(Value), fill = Cell))+geom_col()+
    theme(legend.position="none")+ ggtitle(paste0("Fantom5 Expression Data ", gene))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
  
  b <- ggplot(plotData[1:51,],aes(x=reorder(Cell,-as.numeric(Value)),y=as.numeric(Value),label = Cell,fill = Cell))+
    geom_bar(stat = 'identity')+
    labs(x='', y = 'Expression', title=paste0("Fantom5 Expression Data ", gene))+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  return(list(a,b))      
}



PlotFantonFunction(gene="THY1") 
PlotFantonFunction2(gene="THY1") 
