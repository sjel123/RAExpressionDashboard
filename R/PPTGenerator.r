library(officer)
source("utils.R")

Gene="PRKAB1"

PPTCreate <- function (GeneName="PRKAB1", Index=TRUE, RA=TRUE){
  require(ggplot2)
  require(ggpubr)
  require(Biobase)
  library(officer)
  #doc <- read_pptx("PPT/IISlide.pptx")
  doc <- read_pptx("PPT/IISlide2.pptx")
  doc <- remove_slide(doc, 1)
  doc <- remove_slide(doc, 1)
  top=1.64
  print("slide1")
  
  # Add a slide
  
  CurrentDate <- format(Sys.time(), "%a %b %d %Y")
  
  doc$slideLayouts$get_xfrm_data()[(doc$slideLayouts$get_xfrm_data()[,11]=="Title Slide"),]
  #doc$slideLayouts$get_xfrm_data()[(doc$slideLayouts$get_xfrm_data()[,11]=="Title Slide"),]
  
  doc <- add_slide(doc, layout = "Title Slide 1a", master = "Office Theme")
  doc <- ph_with(x = doc, value = paste0 ("CREDIT report Gene ", GeneName), location = ph_location_type(type = "ctrTitle") ) 
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
    
    ggsave("Recombine.svg",q6, device = "svg", height = 1.75, width = 4.25, scale=3)
    
    
    ggsave("Recombine.png",q6, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("GSE24742.png", q7, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("GSE45867.png",q8, device = "png", height = 1.75, width = 4.25, scale=3)
  }
  
  CurrentDate <- format(Sys.time(), "%a %b %d %Y")
  

  #   layout_summary(doc)

  
    
# add a "Two Content" slide and then content ----
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = paste0("Body Expression ", Gene), location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Gtex.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 1.18, width = 7, height = 2.5 ))
    doc <- ph_with(x = doc, value = external_img("Fantom.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 3.58, width = 7, height = 2.5))
    doc <- ph_with(x = doc, value = external_img("Dice.png"), location = ph_location(left = 7.86, top = 1.04, width = 5, height = 2.5 ))
    doc <- ph_with(x = doc, value = external_img("Protein.png" ), location = ph_location(left = 7.86, top =3.68, width = 5, height = 2.5))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=5.21, top=6.07, width = 4.52, height = .5)) 


if(Index==TRUE){
    # for Metanalysis 1.0
    # if (!file.exists("Meta")) {load("/app/Shiny/IndexTable2/Data/Meta.Rdata")}
    # source("/app/Shiny/RAExpressionDashboard/R/IndexPlotFunction.R")
    # generow1 <- which(Meta$gene==Gene)
    # p <- (PlotFunction(num = generow1, meta = Meta)[[2]])
  
  # for Metanalysis 2.0
      if (!file.exists("Meta")) {load("/app/Shiny/IndexTable2/Data/MetaAnalysis2/Meta2_0.rds")}
      Meta <- Meta1
    source("/app/Shiny/IndexTable2/R/Factor.R")
  generow1 <- which(Meta$gene==Gene)
  p <- (PlotFunction2(num = generow1, meta = Meta)[[2]])
    ggsave("Index.png", p, device = "png", height = 5.56 , width = 6.25, scale=1)
  
  # add a "Two Content" slide and then content ----
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
  doc <- ph_with(x = doc, value = paste0("Disease Expression ", Gene), location = ph_location_type(type = "title") )
  doc <- ph_with(x = doc, value = external_img("Index.png",width = 6.25, height = 5.56 ), location = ph_location(left = 0, top = top-0.5, width = 6.25, height = 5.56 ))
}    
    
    
if(RA==TRUE){      
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("AMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top, width = 5.65, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("LOWAMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top + 2, width = 5.65, height = 2))
    doc <- ph_with(x = doc, value = external_img("Synovial.png"), location = ph_location(left = 5.70, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE55235.png" ), location = ph_location(left = 5.70, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
  
    
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Peak.png"), location = ph_location(left = 0, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("SDY1299.png" ), location = ph_location(left = 0, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
    
    
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Treatment Effects", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Recombine.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE24742.png",width = 5.65, height = 2 ), location = ph_location(left = 0, top = top + 2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = external_img("GSE45867.png"), location = ph_location(left = 5.70, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=6.07, width = 4.52, height = .5)) 
}#End if RA==TRUE  
  }# end i in genename    
  #doc <- ph_with(x = doc, value = iris[1:4, 3:5], location = ph_location_right() )
    #   slide_summary(doc)
    target <- paste0("PPT/RAExpression_", Gene, ".pptx")
  print(doc, target = target)
}


PPTCreate(GeneName = list("SIK1", 'SIK2', "SIK3", "PFKFB3", "CTCS", "ITK", "TMEM173", "FAP", 
                          "RASGRP1", "IRAK1", "IRAK4", "VAV1"), RA=FALSE)   

PPTCreate(GeneName = list("SLC3A1","IRF5", "CDK4", "SLC3A1", "STAT6", "BHLHE40", "BATF"), RA=FALSE)  

PPTCreate(GeneName = list("VAV1","VAV2","VAV3", "RASGRP1","RASGRP2","RASGRP3","RASGRP4"), RA=FALSE) 
PPTCreate(GeneName = list("IL17A"), RA=FALSE, Index = TRUE) 

PPTCreate(GeneName = list("ART3","BMPR1B","ERVMER34-1", "GJA5","KCNH1","MEGF10","NCAM2","NLGN4X","PALM2","PRTG"), RA=FALSE) 

PPTCreate(GeneName = list("IL23A", "IL18", "OSM", "ILR13RA2","IL4R", "IL13",   "IL18BP", "IL18R1", "IL18RP",  "IL23R", 
                          "OSMR"), RA=FALSE) 

PPTCreate(list("IL1RAP","IL1A","IL1B","IL1F10","IL1R1","IL33","IL36A","IL36B","IL36G","IL36RN"), Index = TRUE, RA = FALSE)
PPTCreate( list("EMGF10","ERVMER34-1"), Index = TRUE, RA = FALSE)
PPTCreate( list("IL17A","CSF2"), Index = TRUE, RA = FALSE)
PPTCreate( list("MEGF10"), Index = TRUE, RA = FALSE)
PPTCreate( list("MAP3K8", "IRAK4"), Index = TRUE, RA = TRUE)
PPTCreate( list("MAP3K8", "IRAK4"), Index = TRUE, RA = TRUE)

