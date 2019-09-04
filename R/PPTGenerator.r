library(officer)
source("utils.R")

Gene="PRKAB1"

PPTCreate <- function (GeneName="PRKAB1", RA=TRUE){
  require(ggplot2)
  require(ggpubr)
  library(officer)
  doc <- read_pptx("PPT/IISlide.pptx")
  doc <- remove_slide(doc, 1)
  top=1.64
  print("slide1")
  
  # Add a slide
  
  
  
  doc$slideLayouts$get_xfrm_data()[(doc$slideLayouts$get_xfrm_data()[,11]=="Title Slide"),]
  
  doc <- add_slide(doc, layout = "Title Slide", master = "IandItheme")
  doc <- ph_with(x = doc, value = paste0 ("CREDIT report Gene ", Gene), location = ph_location_type(type = "ctrTitle") ) 
  doc <- ph_with(x = doc, value = CurrentDate,  location = ph_location(left = 1, top =3, width = 4.25, height = 1))
  #doc <- ph_with(x = doc, value = CurrentDate, location = ph_location_type(type= "body") )
  
  
  for (i in GeneName){
    Gene=i
    print(Gene)
  p <- PlotGTEXFunction(gene = Gene)
  p1 <- PlotFantonFunction(gene = Gene)
  p2 <- plotDataFunction(Gene=Gene)
  p3 <- plotImmuneProteomics(gene=Gene)
  
  ggsave("Gtex.png", p[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Fantom.png", p1[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Dice.png", p2, device = "png", height = 1.75 , width = 4.25, scale=3)
  ggsave("Protein.png", p3, device = "png", height = 1.75 , width =4.25, scale=3)
  
  if(RA=="TRUE"){
    q <- PlotAMP(gene = Gene)
    q1 <- PlotLowInputAMP(gene = Gene)
    q2 <- PlotSynovial(gene=Gene) #Disease Progression
    q3 <- PlotGSE55235(gene = Gene) #Healthy OA RA
    
    ggsave("AMP.png",q, device = "png", height = 2 , width = 5.65, scale=3)
    ggsave("LOWAMP.png", q1, device = "png", height = 2 , width = 5.65, scale=3)
    ggsave("Synovial.png", q2, device = "png", height = 1.75 , width = 4.25, scale=3)
    ggsave("GSE55235.png", q3, device = "png", height = 1.75 , width =4.25, scale=3)
    
    q4 <- PlotPeak(gene = Gene)
    q5 <- PlotSDY1299(gene = Gene)
    
    ggsave("Peak.png",q4, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("SDY1299.png", q5, device = "png", height = 1.75, width = 4.25, scale=3)
    
    
    q6 <- PlotRecombine(gene = Gene)
    q7 <- PlotGSE24742(gene = Gene)
    q8 <- PlotGSE45867(gene = Gene)
    
    
    ggsave("Recombine.png",q6, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("GSE24742.png", q7, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("GSE45867.png",q8, device = "png", height = 1.75, width = 4.25, scale=3)
  }
  
  CurrentDate <- format(Sys.time(), "%a %b %d %Y")
  

  #   layout_summary(doc)

  
    
# add a "Two Content" slide and then content ----
  doc <- add_slide(doc, layout = "Title and Content", master = "IandItheme")
    doc <- ph_with(x = doc, value = "Body Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Gtex.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top, width = 5.65, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("Fantom.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top + 2, width = 5.65, height = 2))
    doc <- ph_with(x = doc, value = external_img("Dice.png"), location = ph_location(left = 5.70, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("Protein.png" ), location = ph_location(left = 5.70, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 

if(RA==TRUE){      
  doc <- add_slide(doc, layout = "Title and Content", master = "IandItheme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("AMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top, width = 5.65, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("LOWAMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top + 2, width = 5.65, height = 2))
    doc <- ph_with(x = doc, value = external_img("Synovial.png"), location = ph_location(left = 5.70, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE55235.png" ), location = ph_location(left = 5.70, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
  
    
  doc <- add_slide(doc, layout = "Title and Content", master = "IandItheme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Peak.png"), location = ph_location(left = 0, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("SDY1299.png" ), location = ph_location(left = 0, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
    
    
  doc <- add_slide(doc, layout = "Title and Content", master = "IandItheme")
    doc <- ph_with(x = doc, value = "Treatment Effects", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Recombine.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE24742.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top + 2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = external_img("GSE45867.png"), location = ph_location(left = 5.70, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
}#End if RA==TRUE  
  }# end i in genename    
  #doc <- ph_with(x = doc, value = iris[1:4, 3:5], location = ph_location_right() )
    #   slide_summary(doc)
    target <- paste0("PPT/RAExpression_", gene, ".pptx")
  print(doc, target = target)
}


PPTCreate(GeneName = "CTCS", RA=FALSE)
PPTCreate(GeneName = list("SIK1", 'SIK2', "SIK3", "PFKFB3", "CTCS", "ITK", "TMEM173", "FAP", 
                          "RASGRP1", "IRAK1", "IRAK4", "VAV1"), RA=FALSE)   
