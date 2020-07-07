####HPA Protein
HPAProtein <- read.table("Data/normal_tissue4.tsv", header=T, sep="\t")


head(HPAProtein)

Protein <- "SLC7A9"

PlotProteinHPA <- function(Protein=Protein){
  ProteinIndex <- which(HPAProtein$Gene.name==Protein)
  HPAProtein.f <- HPAProtein[ProteinIndex, ]
ProteinIndex2 <- which(HPAProtein.f$Level!="Not detected")
  HPAProtein.f <- HPAProtein.f[ProteinIndex2, ]
  HPAProtein.f$Level <- factor(HPAProtein.f$Level, levels=c("Low", "Medium", "High"))
  
p <- ggplot(HPAProtein.f, aes(x=interaction(Tissue, Cell.type), y=Level, fill=Cell.type))+geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + ggtitle(paste0(Protein, " Protein Expression in HPA"))
p <- p + scale_fill_manual(
  values = c("red", rep(c("blue", "green", "red"), length(unique(HPAProtein.f$Cell.type))-1))) + theme(legend.position = "none") 
p <- p + xlab("Disease:Cell")


png(filename = paste0(Protein, ".png"),
    width = 750, height = 480)
print(p)
dev.off()

 return(p)
}

PList <-c("SLC7A9", "SLC26A3", "LRRK2", "ACACA", "ACACB", "SMAD3", "TRAF6", "CCRL2","MERTK")
Protein <- "SLC7A9"
Protein <- "SLC26A3"
Protein <- "LRKK2"
Protein <- "ACACA"
Protein <- "SLC3A1"
for (i in 1:9){
  PlotProteinHPA(Protein=PList[i])
}


