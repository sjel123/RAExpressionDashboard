#Scatterplot of adjusted vs unadjusted bulk RA Synovium

# REad Data
Cyto <- read.table("CytoreasonRABulk.txt", header=T, sep="\t", stringsAsFactors = F)
  names(Cyto) <- c("gene_name", "FDRUnadj", "Model", "Value_Unadj", "FDRadj", "Model1", "Value_adj", "gene_name_1")
#Scatterplot
  require(ggplot2)
  library(cowplot)
  Cyto$Color <- ifelse(Cyto$FDRUnadj>1.3 & Cyto$FDRadj>1.3,"green",
                       ifelse(Cyto$FDRUnadj>1.3,2,
                              ifelse(Cyto$FDRadj>1.3,3,4)))
ggplot(Cyto, aes(x=FDRUnadj, y=FDRadj))+
  geom_point(aes(color=Color),show.legend = F)+
  geom_hline(yintercept = 1.3)+
  geom_vline(xintercept = 1.3)+
  annotate(geom="text",x=7.5,y=7.5,label="N=1717")+
  annotate(geom="text",x=0.5,y=7.5,label="N=1776")+
  annotate(geom="text",x=7.5,y=1,label="N=1450")+
  ggtitle("RA Control vs Disease")
  #geom_text(x=7.5,y=7.5,label="N=1717")

Cyto[order(Cyto$FDRadj,decreasing = T),][1:10,]
Cyto[order(Cyto$FDRUnadj,decreasing = T),][1:10,]
