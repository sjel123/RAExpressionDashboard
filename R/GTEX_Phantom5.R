library(Biobase)
library(limma)
library(ggplot2)
library(cowplot)


  
####Realtive Expression in othe tissues using gtex data
load("~/projects/CDH11_membrance/CDH11/Data/combine-expression.rda")#Loaded as dataset
library(biomaRt)
tx <- rownames(dataset[[1]]$expressionData)
  ensembl = useEnsembl(biomart="ensembl")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    res <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "description"),
      filters = "ensembl_gene_id",
      values = tx,
      mart = ensembl)
    head(res)
      #########
      #Tissue enriched data
      #########
      Tissue_Spec <- dataset[["Protein-Atlas"]]$tissueSpecificGenes
      res2 <- merge(Tissue_Spec, res, by.x="Gene", by.y="ensembl_gene_id")
        head(res2)
        res2$description <- gsub("\\ \\[Source\\:HGNC\\ Symbol\\;Acc\\:HGNC\\:\\d+\\]", "", res2$description)
          Gtex <- merge(dataset$`Protein-Atlas`$expressionData, res2, by.x=0, by.y="Gene")
            Specific <- Gtex[,c(37:40)]
            
            #### GETEX DATA 2 (55 Tissues)
            GTEX <- read.table("~/projects/CDH11_membrance/CDH11/Data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", skip=2, header=TRUE, sep="\t")
          #  GTEX_Merge <-  merge(GTEX, Gtex[,31:34], by.x="Description", by.y="external_gene_name", all.x=TRUE)
            
            #Update Gtex data
            Gtex <- merge(GTEX, Specific, by.x="Description", by.y="external_gene_name", all.x=TRUE)
    
          PlotGTEXFunction <- function(gene="HAPLN3") {         
            index <- grep(paste0("^", gene, "$"), Gtex$Description) 
              plotData <-data.frame(t(data.frame(Gtex[index[1],3:55])))
              plotData$Tissue <- row.names(plotData)
                colnames(plotData) <- c("Value", "Tissue")
            
            a <- ggplot(plotData[,], aes(x=Tissue, y=as.numeric(Value), fill = Tissue))+geom_col()+
              theme(legend.position="none")+ ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
              theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
            
            b <- ggplot(plotData[-52,],aes(x=reorder(Tissue,-as.numeric(Value)),y=as.numeric(Value),label = Tissue,fill = Tissue))+
              geom_bar(stat = 'identity')+
              labs(x='', y = 'Expression') + ggtitle(label=paste0("GTEX Expression Data ", gene),subtitle = paste0(Gtex[index,56], Gtex[index,57]))+
              theme_bw()+
              theme(legend.position="none")+
              theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
            return(list(a,b))      
          }      
          
  

########################   
########Fantom5 data
#######################
                Fantom <- read.csv("~/Shiny/Fantom5/GST_FANTOM/TPM_133.cell_median_TPM.csv")
        Fantom_meta <- read.csv("~/Shiny/Fantom5/GST_FANTOM/TPM.Cells.only.samples.csv")
          TissueType <- Fantom_meta$Sample_Type
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
            table(TissueType_sub )
            
          
            Fantom_meta$CellType <- TissueType_sub
            Fantom_meta$Sample_Type <- gsub("\\.\\.\\.|\\.\\.","\\.", Fantom_meta$Sample_Type)
            Fantom_meta$Sample_Type <- gsub("Smooth.Muscle.Cells.Umbilical.artery.","Smooth.Muscle.Cells.Umbilical.Artery.", Fantom_meta$Sample_Type)
            
              unique(Fantom_meta$Sample_Type) %in%
              colnames(Fantom)
              Fantom_meta <- Fantom_meta[(!duplicated(Fantom_meta$Sample_Type)),]

              
              MeanCluster_Fantom <- apply(Fantom[,-1],1,function(x) tapply(x, Fantom_meta$CellType, mean))
              MeanCluster_Fantom <- data.frame(MeanCluster_Fantom)
              MeanCluster_Fantom_t <- as.data.frame(t(MeanCluster_Fantom))
              MeanCluster_Fantom_t$GeneSymbol <- Fantom$X
              head(MeanCluster_Fantom_t)
              
###Plot Data
  
    PlotFantonFunction <- function(gene="HAPLN3") {         
        index <- grep(paste0(gene,"$"), Fantom$X) 
          plotData <- data.frame(t(MeanCluster_Fantom_t[index,]))
          plotData[] <- lapply(plotData, as.character)
          colnames(plotData) <- (plotData[52,])
          plotData$Cell <- row.names(plotData)
            colnames(plotData) <- make.names(colnames(plotData))
            plotData$Cell <- (row.names(plotData))
              index2 <- grep("p1\\.|Cell", colnames(plotData))
              plotData <- plotData[,index2]
              
              colnames(plotData) <- gsub("p1\\.", "", colnames(plotData))
                index3 <- grep(paste0("^", gene, "$"), colnames(plotData) )
                colnames(plotData)[index3] <- "Value"
                     
          a <- ggplot(plotData[-52,], aes(x=Cell, y=as.numeric(Value), fill = Cell))+geom_col()+
            theme(legend.position="none")+ ggtitle(paste0("Fantom5 Expression Data ", gene))+
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))
          
          b <- ggplot(plotData[-52,],aes(x=reorder(Cell,-as.numeric(Value)),y=as.numeric(Value),label = Cell,fill = Cell))+
            geom_bar(stat = 'identity')+
            labs(x='', y = 'Expression', title=paste0("Fantom5 Expression Data ", gene))+
            theme_bw()+
            theme(legend.position="none")+
            theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
       return(list(a,b))      
    }
    
    
    
    PlotFantonFunction(gene="THY1") 
    PlotGTEXFunction(gene="THY1")
    AMPSCPlotFunction(gene="THY1")
   
    

    LowInput_All <- fData(lowinputGSET)
    save (Gtex, Fantom, MeanCluster_Fantom_t, file="Data/MembraneDataAll.RData")
    load("Data/MembraneDataAll.RData")
    