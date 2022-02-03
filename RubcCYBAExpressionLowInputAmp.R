Gene1="LRRC15"
Gene1="MEGF10"
Gene1="RUBCN"
Gene1="IRAK2"
library(cowplot)
library(Biobase)
p <- PlotLowInputAMP(gene = Gene1)
p <- p + geom_hline(yintercept =4.5, color="red", linetype = "dashed")
p <-  p + geom_hline(yintercept =1, color="red", linetype = "dashed")
p


p <- PlotFIBData(gene = Gene1)
p <- p + geom_hline(yintercept =0.5, color="red", linetype = "dashed")
p


p <- PlotSynovial(gene = Gene1)

p


AMPPahse1 <- function(gene="RUBCN"){
  require(cowplot)
  require(Biobase)
  require(ggpubr)
  if(!exists("lowinputGSET")) {
    load("~/projects/CDH11_membrance/CDH11/Data/ProteinClass/Data/LowInput.RData")
  }
  
  PlotData <- exprs(lowinputGSET)[fData(lowinputGSET)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(lowinputGSET))
  a <-  ggplot(PlotData_t,aes(x=Disease,y=Value,label = DiseaseTissue,fill = Disease)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_point(pch=21, position = position_jitterdodge(jitter.width = 0.2))+
    labs(x='', y = 'Expression', title=paste0("AMP Low Input Expression Data ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a <- a + facet_wrap(~Cell.type, nrow = 1, strip.position="bottom")+
    theme(strip.text.x = element_text(angle=0))+stat_compare_means(aes(group=Disease), method= "t.test")+
     theme(axis.text.x=element_blank())+  theme(
      #aspect.ratio = 41,
      strip.text = element_text(hjust = .5, angle=0, size=16),
      strip.background = element_blank(),
      strip.placement = "outside")
  return(a)
}

AMPPahse1(gene = "RUBCN")
AMPPahse1(gene = "CYBA")
gg <- ggplot(plotdata2, aes(x=Disease, y=value,
                            color = Disease, fill=Disease, palette = "jco"))+ geom_boxplot() + facet_wrap(~Tissue, nrow = 1, strip.position="bottom")+
  theme(strip.text.x = element_text(angle=90))
gg <- gg + stat_compare_means(comparisons = list(c("Healthy", "RA")), label.y = maxval+1)#comparisons = my_comparisons, label.y= maxval) # Add pairwise comparisons p-value
gg <- gg  +  theme (text=element_text(size=16), 
                    axis.text.x = element_text(angle = 0, hjust =1, size=16))+geom_point(position=position_dodge(width=0.7))
gg <- gg + ggtitle(fData(ESET)[generow,2])
gg <- gg + theme(axis.text.x=element_blank())+  theme(
  #aspect.ratio = 41,
  strip.text = element_text(hjust = 1),
  strip.background = element_blank(),
  strip.placement = "outside"
)




RubPlot <- function(gene="RUBCN"){
  require(cowplot)
  require(Biobase)
  require(ggpubr)
  textsize=16
  
  PlotData <- exprs(gset)[fData(gset)$DFnames==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(gset))
  PlotData_t <-  PlotData_t [PlotData_t$Group%in% c("iTNF_1",  "iTNF_2","MTX_1" , "MTX_2"),]
  a <-  ggplot(PlotData_t,aes(x=Time,y=Value,fill = Time)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_point(pch=21, position = position_jitterdodge(jitter.width = 0.2))+
    labs(x='', y = 'Expression', title=paste0("KI ReCombine WholeBlood ",gene))+
    theme_bw()+
    #theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a
  a <- a + facet_wrap(~Treatment, nrow = 1, strip.position="bottom")+
    theme(strip.text.x = element_text(angle=0))+stat_compare_means(aes(group=Time), method= "t.test")+
    theme(axis.text.x=element_blank())+  theme(
      #aspect.ratio = 41,
      strip.text = element_text(hjust = .5, angle=0, size=16),
      strip.background = element_blank(),
      strip.placement = "outside")
  a <- a+ scale_fill_manual(name="Time",
                            labels=c("baseline","3 months"),
                            values=c("#F8766D", "#00BFC4"))
  return(a)
}

PlotGSE24742(gene="RUBCN")
RubPlot(gene="IRAK2")
AMPPahse1(gene="IRAK2")
PlotLowInputAMPII2()


PlotLowInputAMPII2 <- function(gene="CDH11"){
  if(!exists("lowinputAMPPhaseII")) {
    load("~/Projects/AMP_Phase2_RA/Data/AMPPhase2.RData")
    lowinputAMPPhaseII <- gset
    rm(gset)  
  }
  
  PlotData <- exprs(lowinputAMPPhaseII)[fData(lowinputAMPPhaseII)$gene_name==gene,]
  PlotData_t <-data.frame(Value=PlotData, pData(lowinputAMPPhaseII))
  a <-  ggplot(PlotData_t,aes(x=Diagnosis,y=Value,label = Group, fill=Diagnosis)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_point(pch=21, position = position_jitterdodge(jitter.width = 0.2))+
    labs(x='', y = 'Expression')+
    theme_bw()+
    theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  a 
  a <- a + facet_wrap(~cell_type, nrow = 1, strip.position="bottom")+
    theme(strip.text.x = element_text(angle=0))+stat_compare_means(aes(group=Diagnosis), method= "t.test")+
    theme(axis.text.x=element_blank())+  theme(
      #aspect.ratio = 41,
      strip.text = element_text(hjust = .5, angle=0, size=16),
      strip.background = element_blank(),
      strip.placement = "outside")
  a <- a + ggtitle(paste0("AMP2.0 Low input Expression Data ", gene))
  return(a)
}
  
PlotLowInputAMPII2a <- function(gene="CDH11") {
  if(!exists("lowinputAMPPhaseII")) {
    load("~/Projects/AMP_Phase2_RA/Data/AMPPhase2.RData")
    lowinputAMPPhaseII <- gset
    rm(gset)  
  }

PlotData <- exprs(lowinputAMPPhaseII)[fData(lowinputAMPPhaseII)$gene_name==gene,]
PlotData_t <-data.frame(Value=PlotData, pData(lowinputAMPPhaseII))
  PlotData_t$Cohort <- gsub(" APPROVED", "", PlotData_t$Cohort)
  #PlotData_t$Cohort <- gsub("", "OA", PlotData_t$Cohort)
  Index <- grep("Group", PlotData_t$Cohort, invert = T)
  PlotData_t$Cohort[Index] <- "OA"
  
  b <-  ggplot(PlotData_t,aes(x=Cohort,y=Value,fill = Cohort)) + 
    #geom_bar(stat = "summary", fun.y = "median", width=0.8)+
    geom_boxplot(outlier.shape = NA)+ geom_point(pch=21, position = position_jitterdodge(jitter.width = 0.2))+
    labs(x='', y = 'Expression')+
    theme_bw()+
    # theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = textsize),panel.grid.major= element_blank(), panel.grid.minor = element_blank())
  b <- b + facet_grid(~cell_type)
  b <-  b  + scale_x_discrete(labels=c("Group 1" = "Drug Naive", "Group 2" = "MTX Fail",
                                    "Group 3" = "TNF Fail"))
  b <-  b + stat_compare_means(aes(group=Cohort))+ ggtitle( paste0("AMP Phase II Low Input Expression Data ",gene))
 

  
  
  return(b)
}
PlotLowInputAMPII2a()
