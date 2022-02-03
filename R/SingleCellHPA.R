#######
#Single Cell Data from HPA
######

SingleCelllot <- function(Gene ="MET"){
  require(dplyr)
  require(ggplot2)
  require(cowplot)
  if(!exists("SingleDF")){
    SingleDF <- read.table("Data/rna_single_cell_type.tsv", header = T, sep="\t", stringsAsFactors = F )
    SingleMeta <- read.table("Data/SingleCellMeta.txt",  header = T, sep="\t", stringsAsFactors = F)
  }
    SingleDF1 <- filter(SingleDF, Gene.name == UQ(Gene))
      SingleDF1 <- merge(SingleDF1, SingleMeta, by.x="Cell.type", by.y="Cell.type")
      SingleDF1$Cell.type <- factor(SingleDF1$Cell.type , levels=(SingleDF1$Cell.type[order(SingleDF1$Cell.type.group)]))
     p <- ggplot(SingleDF1, aes(x=Cell.type, y=NX, fill=Cell.type.group))+geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(Gene)
      p
     return(p)
      }


SingleCelllot(Gene = "MET")
SingleCelllot(Gene = "CD74")
