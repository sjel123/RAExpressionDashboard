PlotFunction <- function(num=1, meta=Meta){
  require(ggplot2)
  require(cowplot)
  require(ggrepel)
  DF <- meta[num,]
  DF1 <- data.frame()
  index1 <- grep("Size$", names(meta))
  index2 <- grep("FDR", names(meta))
  FC= DF[index1]
  pval = DF[index2]
  DF1 <- as.data.frame(cbind(as.numeric(FC), as.numeric(pval)))
  #DF1$Tissue <- DF1$label
  DF1$Tissue <- names(DF[index1])
  DF1$Tissue <- gsub("_effectSize", "", DF1$Tissue)
  DF1$label<-DF1$Tissue
  DF1$Disease <- DF1$Tissue
  DF1$Tissue[grep('blood', DF1$label)]="Blood"
  DF1$Tissue[grep('blood', DF1$label,invert = T)]="Tissue"
  
  row.names(DF1) <- colnames(FC)
  row.names(DF1) <- gsub("_effectSize", "", row.names(DF1))
  colnames(DF1) <- c("FC", "pVal", "Tissue", "Label", "Disease")
  DF1$label <- rownames(DF1)
  temp1 <- strsplit(DF1$label,split = "_",fixed = TRUE)
  DF1$Disease <- sapply(temp1, "[",1)
  
  DF1$Label2 <- gsub("_colon|_blood|_synovial|_skin|_liver", "",DF1$label)
  # All labels should be to the right of 3.
  x_limits <- c( NA,min(DF1$FC, na.rm = TRUE)+1)
  g <- ggplot(DF1, aes(x=FC, y=-log10(pVal), label = Label2, size=4, shape=Tissue))+
    geom_point(aes(color=Disease,size=4))+
    xlim(min(DF1$FC)-1, NA) +
    geom_text_repel(aes(color=Disease),
                    arrow = arrow(length = unit(0.02, "npc"), type = "open", ends = "last"),
                    force = 3,
                    direction = "y",
                    #xlim=x_limits,
                    nudge_x  = min(DF1$FC, na.rm = T)-1,
                    hjust = 0
                    
    )
  g <- g + geom_hline(yintercept = 1.3, linetype="dashed") 
  g <- g + geom_vline(xintercept = 0)
  g <- g + ggtitle(meta$gene[num]) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  g <- g + scale_size(guide=FALSE) + theme(legend.text = element_text(colour="blue", size=16)) +
    guides(shape = guide_legend(override.aes = list(size = 5)))+
    guides(color = guide_legend(override.aes = list(size = 5)))
  g <- g+facet_wrap(~Disease,nrow = 3)
  #+ theme(legend.key.size = unit(3,"line")
  DF1 <- DF1[order(row.names(DF1)),]
  DF1$Signif <- ifelse(DF1$pVal>0.001|is.na(DF1$pVal),"FALSE","TRUE")
  DF1$Bold   <- ifelse(DF1$Signif, "bold", "plain")
  DF2 <- DF1[c(-1,-2,-3),]
  #Eliminate Matt AD analysis from figure
  gg <- ggplot(DF2, aes(x=FC, y=label, size= -log10(pVal), color=Signif))+geom_point()+ggtitle(meta[num,1])
  gg <- gg + geom_vline(xintercept = 0, color="red", alpha = 0.4, linetype="dashed")
  gg <- gg  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gg <- gg  + theme(legend.text = element_text(colour="blue", size=16),
                    panel.border = element_rect(color = "black", fill=NA, size=1, linetype = 1, inherit.blank = FALSE)) 
  gg <- gg + guides(color=FALSE)
  gg <- gg +  scale_colour_manual(values = c("blue", "red"))
  gg <- gg + theme(axis.text.y = element_text(face = DF2$Bold))
  gg <- gg + geom_hline(yintercept = c(3.5, 5.5,8.5, 12.5, 14.5, 15.5 ), linetype="dashed", alpha = 0.2)
  gg 
  
  
  my_list=list("g"=g, "gg"=gg)
  return(my_list)                                                            
}