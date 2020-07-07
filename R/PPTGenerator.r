library(officer)
source("utils.R")

Gene="PRKAB1"

PPTCreate <- function (GeneName="PRKAB1", Index=TRUE, RA=TRUE){
  require(ggplot2)
  require(ggpubr)
  require(Biobase)
  require(gplots)
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
  p4 <- HPABloodPlot(gene=Gene)
  
  ggsave("Gtex.png", p[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Fantom.png", p1[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Dice.png", p2, device = "png", height = 1.75 , width = 4.25, scale=3)
  ggsave("Protein.png", p3, device = "png", height = 1.75 , width =4.25, scale=3)
  ggsave("HPA.png", p4, device = "png", height = 1.75 , width =4.25, scale=3)
  
  if(RA=="TRUE"){
    q  <- try(PlotAMP(gene = Gene), silent=T)
     if(class(q) == "try-error"){
      q <- qplot(x=1,y=1, label="Gene Not Detected", geom="text")
      }

    q1 <- try(PlotLowInputAMP(gene = Gene), silent=T)
      if(class(q1) == "try-error"){
        q1 <- qplot(x=1,y=1, label="Gene Not Detected", geom="text")
      }
    q2 <- PlotSynovial(gene=Gene) #Disease Progression
    q3 <- PlotGSE55235(gene = Gene) #Healthy OA RA
    
    ggsave("AMP.png",q, device = "png", height = 2 , width = 5.65, scale=3)
    ggsave("LOWAMP.png", q1, device = "png", height = 2 , width = 5.65, scale=3)
    ggsave("Synovial.png", q2, device = "png", height = 1.75 , width = 4.25, scale=3)
    ggsave("GSE55235.png", q3, device = "png", height = 1.75 , width =4.25, scale=3)
    
    q4 <- try(PlotPeak(gene = Gene), silent = T)
    if(class(q4) == "try-error"){
      q4 <- qplot(x=1,y=1, label="Gene Not Detected", geom="text")
    }
    q5 <- PlotSDY1299(gene = Gene)
    
    ggsave("Peak.png",q4, device = "png", height = 1.75, width = 4.25, scale=3)
    ggsave("SDY1299.png", q5, device = "png", height = 1.75, width = 4.25, scale=3)
    
    
    q6 <- try(PlotRecombine(gene = Gene), silent = T)
    if(class(q6) == "try-error"){
      q6 <- qplot(x=1,y=1, label="Gene Not Detected", geom="text")
    }
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
  # first define a fpar ----
  # fp <- fpar(
  #   ftext(c("Un titre", "Deux titre"), fp_text(bold = TRUE, font.size = 20))
  # )
  # first define a fpar ----
  fp <- unordered_list(
    level_list = c(1, 1),
    str_list = c("Level1", "Level2"),
    style = fp_text(color = "red", font.size = 20) )
  small_text <- fp_text(color = "red", font.size = 18)
  smaller_text <- fp_text(color = "red", font.size = 12)
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
  doc <- ph_with(x = doc, value = paste0("Body Expression ", Gene), location = ph_location_type(type = "title") )
  
  
  
  
  doc <- ph_with(x = doc, value = external_img("Gtex.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 1.18, width = 7, height = 2.5 ))
  doc <- ph_with(x = doc, value = external_img("Fantom.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 3.58, width = 7, height = 2.5))
  doc <- ph_with(x = doc, value = external_img("Dice.png"), location = ph_location(left = 7.86, top = .73, width = 4, height = 2. ))
  doc <- ph_with(x = doc, value = external_img("Protein.png" ), location = ph_location(left = 7.86, top = 2.69, width = 4, height = 2.))
  doc <- ph_with(x = doc, value = external_img("HPA.png" ), location = ph_location(left = 7.86, top =4.73, width = 4, height = 2.))
  doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 
  
  styles <- fp_text(color = "black", font.size = 16, font.family = "Arial", bold=FALSE)
  ul <- unordered_list(
    level_list = c(1),
    str_list = c("RNA"),
    style = styles)
  
  # fp1 <- fpar(
  #   ftext("RNA", fp_text(bold = FALSE,color = "black", font.size = 16)))
  
  
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 6.66, top = 1.13, width = 0.73, height = 0.4 ))
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 6.66, top = 3.33, width = 0.73, height = 0.4 ))
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 11.88, top = 0.83, width = 0.73, height = 0.4 ))
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 11.88, top = 4.85, width = 0.73, height = 0.4 ))
  
  ul$str <- "Protein"
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 11.88, top = 2.95, width = 0.98, height = .59 ))
  
  ul <- unordered_list(
    level_list = c(1),
    str_list = c("Level1"),
    style = fp_text(color = "black", font.size = 12, font.family = "Arial", bold=FALSE))
  ul$str="https://dice-database.org/"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 1.24, width = 1.5, height = .5))
  ul$str="Nature Immunology vol 18, pages583–593(2017)"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 3.47, width =1.5, height = .83))
  ul$str="https://www.proteinatlas.org/"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 5.44, width =1.5, height = .45))
  doc <- ph_with(x = doc, value = external_img("GTEXLogo.png", width = 1.43, height = .19 ), location = ph_location(left = 5.15, top = 1.24, width = 1.43, height = .19))
  doc <- ph_with(x = doc, value = external_img("FantomLogo.png", width = .39, height =0.48), location = ph_location(left = 5.87, top =3.26, width = .39, height =0.48))
  
if(Index==TRUE){
    # for Metanalysis 1.0
    # if (!file.exists("Meta")) {load("/app/Shiny/IndexTable2/Data/Meta.Rdata")}
    # source("/app/Shiny/RAExpressionDashboard/R/IndexPlotFunction.R")
    # generow1 <- which(Meta$gene==Gene)
    # p <- (PlotFunction(num = generow1, meta = Meta)[[2]])
  
  # for Metanalysis 2.0
      if (!file.exists("Meta")) {load("/app/Shiny/IndexTable2/Data/MetaAnalysis2/Meta2_0.rds")}
      Meta <- Meta1
    #source("/app/Shiny/IndexTable2/R/Factor.R")
      source("MetaDataPlots.R")
  generow1 <- which(Meta$gene==Gene)
  p <- (PlotFunctionMeta(num = generow1, meta = Meta)[[3]])
  p1 <- (PlotFunctionMeta(num = generow1, meta = Meta)[[4]])
    ggsave("MetaTissue.png", p1, device = "png", height = 5.56 , width = 6.25, scale=1)
    ggsave("MetaBlood.png", p, device = "png", height = 5.56 , width = 6.25, scale=1)
  # add a "Two Content" slide and then content ----
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
  doc <- ph_with(x = doc, value = paste0("Disease Expression ", Gene), location = ph_location_type(type = "title") )
  doc <- ph_with(x = doc, value = external_img("MetaTissue.png",width = 6.25, height = 5.56 ), location = ph_location(left = 0, top = 1.32, width = 6.25, height = 5.56 ))
  doc <- ph_with(x = doc, value = external_img("MetaBlood.png",width = 6.25, height = 5.56 ), location = ph_location(left = 6.61, top = 1.32, width = 6.25, height = 5.56 ))
  doc <- ph_with(x = doc, value = c("Tissue"), location = ph_location(left=2.65, top=.91, width = 4.52, height = .5)) 
  doc <- ph_with(x = doc, value = c("Blood"), location = ph_location(left=9.17, top=0.91, width = 4.52, height = .5)) 
  doc <- ph_with(x = doc, value =  fp, location = ph_location(left=4.45, top=6.0, width = 4.17, height = 1.45)) 
}  
    
if(RA==TRUE){      
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("AMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top, width = 5.65, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("LOWAMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top + 2, width = 5.65, height = 2))
    doc <- ph_with(x = doc, value = external_img("Synovial.png"), location = ph_location(left = 7.08, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE55235.png" ), location = ph_location(left = 7.08, top =top+2, width = 4.25, height = 2))
    #doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
    doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 
    
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Peak.png"), location = ph_location(left = 0.5, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("SDY1299.png" ), location = ph_location(left = 0.5, top =top+2, width = 4.25, height = 2))
    #doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
    doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 
    
  doc <- add_slide(doc, layout = "1_One Column Text", master = "Office Theme")
    doc <- ph_with(x = doc, value = "Treatment Effects", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Recombine.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE24742.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top + 2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = external_img("GSE45867.png"), location = ph_location(left = 7.08, top = top, width = 4.25, height = 2 ))
   # doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=6.07, width = 4.52, height = .5)) 
    doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 
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
PPTCreate( list("CD109"), Index = TRUE, RA = TRUE)
PPTCreate( list("SMAD3", "ACACA", "CCRL2", "SLC3A1","SLC7A9"), Index = TRUE, RA = TRUE)
PPTCreate( list("SLC7A9", "SLC26A7", "LRRK2", "ACACA", "ACACB", "SMAD3", "TRAF6", "CCRL2","MERTK"), Index = TRUE, RA = TRUE)
PPTCreate( list("SLC3A1"), Index = TRUE, RA = FALSE)
PPTCreate( list("TNFSF10", "TNFRSF10A","TNFRSF10B","TNFRSF10C","TNFRSF10D", "TNFRSF6B"), Index = TRUE, RA = FALSE)
PPTCreate( list("KIAA0226"), Index = TRUE, RA = TRUE)
PPTCreate( list("TNFRSF1B"), Index = TRUE, RA = TRUE)
PPTCreate( list("PADI2", "PADI4"), Index = TRUE, RA = TRUE)
warnings()
