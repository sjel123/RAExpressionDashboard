# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Nov 15 14:39:43 EST 2019


Server: www.ncbi.nlm.nih.gov
Query: acc=GSE93777&platform=GPL570&type=txt&groups=All&colors=dfeaf4&selection=0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000&padj=fdr&logtransform=auto&columns=ID&columns=adj.P.Val&columns=P.Value&columns=F&columns=Gene+symbol&columns=Gene+title&num=250&annot=ncbi

# Unable to generate script analyzing differential expression.
#      Invalid input: at least two groups of samples should be selected.

################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE93777", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


gset$Tissue <- paste0(gset$`tissue:ch1`, gset$`cell type:ch1`)
gset$Tissue <- gsub("NA", "", gset$Tissue)
gset$Disease <- paste0(gset$`disease:ch1`, gset$`disease state:ch1`)
gset$Disease <- gsub("NA", "", gset$Disease)

library(hgu133plus2.db)

## Bimap interface:
x <- hgu133plus2GENENAME
# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the GENENAME for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
Description <- data.frame(unlist(xx))


## Bimap interface:
x <- hgu133plus2SYMBOL
# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the GENENAME for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
Symbol <- data.frame(unlist(xx))

Annotation <- data.frame(Symbol, Description)

  names(Annotation) <- c("Symbol", "Description")
  Annotation$Probe <- row.names(Annotation)
  table(table(row.names(Annotation)))

  
   library(dplyr)
  ProbeId <- data.frame(Probe=row.names(exprs(gset)))
 Annotation2 <-left_join(ProbeId, Annotation)
    head(Annotation2)
  table(row.names(exprs(gset))== ProbeId$Probe)
    fData(gset) <- Annotation2
    
 Prob  <-  row.names(rbind(Annotation[which(Annotation$Symbol=="MAP3K8"),],
        Annotation[which(Annotation$Symbol=="IRAK4"),]))
 index <- which (row.names(exprs(gset))%in% Prob) 
 PlotData <- data.frame(t(exprs(gset[index,])), pData(gset))

library(ggplot2)
library(cowplot)
 ggplot(PlotData, aes(y=X205027_s_at, color=Disease,x=Tissue, fill=Disease))+geom_boxplot(outlier.shape = NA)+
   geom_point(position=position_dodge(width = .7))+
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("TPL2 is RA") #MAP3K8
 
 ggplot(PlotData, aes(y=X219618_at, color=Disease, x=Tissue, fill=Disease))+geom_boxplot(outlier.shape = NA)+
   geom_point(position=position_dodge(width = .7))+
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("IRAK4 is RA") #Irak4
 
 
 ggplot(PlotData, aes(y=X235421_at, color=Disease, x=interaction(Tissue, Disease)))+geom_point()+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

 
 
 # group names for all samples in a series
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "000000000000000000000000000000000000000000000000")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  set group names

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("All")

# set parameters and draw the plot
palette(c("#dfeaf4", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE93777", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
