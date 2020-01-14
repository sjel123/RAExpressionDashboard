#HPA Plots

HPABloodPlot <- function(gene="IL6"){
  require(ggplot2)
  if(!exists("HPABloodgset")) {
    load("/app/Shiny/RAExpressionDashboard/Data/HPABloodGset.rdata")
  }
  
  index <- grep( gene, fData(HPABloodgset)$Gene.name) 
  plotData <-(data.frame(exprs(HPABloodgset)[index[1],]))
  plotData$Cell <- pData(HPABloodgset)$Cell
  plotData$State <- pData(HPABloodgset)$Donor
  colnames(plotData) <- c("Value", "Cell", "Donor")
  
  a <- ggplot2::ggplot(plotData[,], aes(x=as.character(Cell), y=as.numeric(Value)))+geom_boxplot()+
    ggtitle(label=paste0("Immune Cell Protein Expression Data ", gene)) +
    geom_point(position=position_dodge(width=.75))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
  return(a)      
}      


SchmiedelBloodPlot <- function(gene="IL6"){
  require(ggplot2)
  if(!exists("Schmiedelgset")) {
    load("/app/Shiny/RAExpressionDashboard/Data/Schmiedelgset.rdata")
  }
  
  index <- grep( gene, fData(Schmiedelgset)$Gene.name) 
  plotData <-(data.frame(exprs(Schmiedelgset)[index[1],]))
  plotData$Cell <-Schmiedelgset$Cell
  # plotData$State <- pData(HPABloodgset)$Donor
  colnames(plotData) <- c("Value", "Cell")
  
  a <- ggplot2::ggplot(plotData[,], aes(x=as.character(Cell), y=as.numeric(Value)))+geom_boxplot()+
    ggtitle(label=paste0("Immune Cell Protein Expression Data ", gene)) +
    geom_point(position=position_dodge(width=.75))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
  return(a)      
}      

MonacoBloodPlot <- function(gene="IL6"){
  require(ggplot2)
  if(!exists("Schmiedelgset")) {
    load("/app/Shiny/RAExpressionDashboard/Data/Monacogset.rdata")
  }
  
  index <- grep( gene, fData(Monacogset)$Gene.name) 
  plotData <-(data.frame(exprs(Monacogset)[index[1],]))
  plotData$Cell <-Monacogset$Cell
  # plotData$State <- pData(HPABloodgset)$Donor
  colnames(plotData) <- c("Value", "Cell")
  
  a <- ggplot2::ggplot(plotData[,], aes(x=as.character(Cell), y=as.numeric(Value)))+geom_boxplot()+
    ggtitle(label=paste0("Immune Cell Protein Expression Data ", gene)) +
    geom_point(position=position_dodge(width=.75))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size = textsize))
  return(a)      
} 

ImmunePlot <- function(gene="AZU1"){
  require(cowplot)
  A=HPABloodPlot(gene)
  B=SchmiedelBloodPlot(gene)
  C=MonacoBloodPlot(gene)
  top_row <- plot_grid(A,B, labels = c( 'A', "B"), label_size = 12)
  plot_grid(top_row, C, ncol=1)
}