PlotFunctionMeta <- function(num=1, meta=Meta){
  require(ggplot2)
  require(cowplot)
  require(ggrepel)
    DF <- meta[num,]
    DF1 <- data.frame()
      index1 <- grep("Size$|Size.[xy]", names(meta))
      index2 <- grep("FDR", names(meta))
      index3 <- grep("Pval", names(meta))
        FC= DF[index1]
        pval = DF[index2]
        het = DF[index3]
    DF1 <- as.data.frame(cbind(as.numeric(FC), as.numeric(pval), as.numeric(het)))
    #DF1$Tissue <- DF1$label
      DF1$Tissue <- names(DF[index1])
        DF1$Tissue <- gsub("_effectSize", "", DF1$Tissue)
      DF1$label<-DF1$Tissue
      DF1$Disease <- DF1$Tissue
        DF1$Tissue[grep('blood', DF1$label)]="Blood"
        DF1$Tissue[grep('blood', DF1$label,invert = T)]="Tissue"
    
      row.names(DF1) <- colnames(FC)
      row.names(DF1) <- gsub("_effectSize", "", row.names(DF1))
      colnames(DF1) <- c("FC", "pVal", "Hetp", "Tissue", "Label", "Disease")
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
    g <- g + xlab("Effect Size")
    
  
  #+ theme(legend.key.size = unit(3,"line")
  DF1 <- DF1[order(row.names(DF1)),]
    DF1$Signif <- ifelse(DF1$pVal>0.01|is.na(DF1$pVal),"black", "red")
    DF1$SignHet <- ifelse(DF1$Hetp>0.01|is.na(DF1$Hetp),"FALSE", "TRUE")
    DF1$Bold   <- ifelse(DF1$Signif=="red", "bold", "plain")
    #DF2 <- DF1[c(-1,-2,-3),]
    DF2=DF1
  #Eliminate Matt AD analysis from figure
  gg <- ggplot(DF2, aes(x=FC, y=label, size= -log10(pVal),color=Signif, shape=SignHet))+geom_point(  )+ggtitle(meta[num,1])
  
    gg 
    gg <- gg + geom_vline(xintercept = 0, color="red", alpha = 0.4, linetype="dashed")
    gg <- gg  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    gg <- gg  + theme(legend.text = element_text(colour="blue", size=16),
                      panel.border = element_rect(color = "black", fill=NA, size=1, linetype = 1, inherit.blank = FALSE)) 
    gg <- gg + guides(color=FALSE)
    gg <- gg +  scale_colour_manual(values = c("blue", "red"))
    gg <- gg + theme(axis.text.y = element_text(face = DF2$Bold))
    gg <- gg + geom_hline(yintercept = c(3.5,  5.5,  9.5, 11.5, 15.5, 16.5, 18.5, 19.5, 21.5,22.5, 24.5, 27.5, 30.5, 31.5, 32.5 ), linetype="dashed", alpha = 0.2)
    gg <- gg + xlab("Effect Size")
  gg
  #Use only Blood Samples
  DF3 <- DF2[DF2$Tissue=="Blood",]
  ggg <- ggplot(DF3, aes(x=FC, y=label, size= -log10(pVal),color=Signif, shape=SignHet))+geom_point(  )+ggtitle(meta[num,1])
    
    ggg 
    ggg <- ggg + geom_vline(xintercept = 0, color="red", alpha = 0.4, linetype="dashed")
    ggg <- ggg  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggg <- ggg  + theme(legend.text = element_text(colour="blue", size=16),
                      panel.border = element_rect(color = "black", fill=NA, size=1, linetype = 1, inherit.blank = FALSE)) 
    ggg <- ggg + guides(color=FALSE)
    ggg <- ggg +  scale_colour_manual(values = c("blue", "red"))
    ggg <- ggg + theme(axis.text.y = element_text(face = DF3$Bold)) + xlab("Effect Size")+ylab("")+
      labs(shape = "SignHet (p<0.01)")
   # ggg <- ggg + geom_hline(yintercept = c(3.5,  5.5,  9.5, 11.5, 15.5, 16.5, 18.5, 19.5, 21.5,22.5, 24.5, 27.5, 30.5, 31.5, 32.5 ), linetype="dashed", alpha = 0.2)
    g3 <-ggg 
    
  #Use Only Tissue Samples
  DF3 <- DF2[DF2$Tissue=="Tissue",]
    ggg <- ggplot(DF3, aes(x=FC, y=label, size= -log10(pVal),color=Signif, shape=SignHet))+geom_point(  )+ggtitle(meta[num,1])
    
    ggg 
    ggg <- ggg + geom_vline(xintercept = 0, color="red", alpha = 0.4, linetype="dashed")
    ggg <- ggg  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggg <- ggg  + theme(legend.text = element_text(colour="blue", size=16),
                        panel.border = element_rect(color = "black", fill=NA, size=1, linetype = 1, inherit.blank = FALSE)) 
    ggg <- ggg + guides(color=FALSE)
    ggg <- ggg +  scale_colour_manual(values = c("blue", "red"))
    ggg <- ggg + theme(axis.text.y = element_text(face = DF3$Bold))
      BRESKS= as.numeric(cumsum(table(DF3$Disease))+0.5)
    ggg <- ggg + geom_hline(yintercept = BRESKS, linetype="dashed", alpha = 0.2)+ xlab("Effect Size")+ylab("")+
                labs(shape = "SignHet (p<0.01)") + guides(shape = guide_legend(override.aes = list(size=3)))
    g4 <-ggg 
  my_list=list("g"=g, "gg"=gg, "g3"=g3, "g4"=g4)
  return(my_list)                                                            
}





# Selected lowest p per tissue



PlotFunctionMeta(num = 8394)
