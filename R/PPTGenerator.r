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
  doc <- read_pptx("PPT/IISlide3.pptx")
  doc <- remove_slide(doc, 1)
  doc <- remove_slide(doc, 1)
  top=1.64
  print("slide1")
  
  # Add a slide
  
  CurrentDate <- format(Sys.time(), "%a %b %d %Y")
  
  doc$slideLayouts$get_xfrm_data()[(doc$slideLayouts$get_xfrm_data()[,11]=="Title Slide"),]
  #doc$slideLayouts$get_xfrm_data()[(doc$slideLayouts$get_xfrm_data()[,11]=="Title Slide"),]
  
  doc <- add_slide(doc, layout = "Title Slide", master = "II2021")
  doc <- ph_with(x = doc, value = paste0 ("CREDIT report Gene ", GeneName), location = ph_location_type(type = "ctrTitle") ) 
  doc <- ph_with(x = doc, value = CurrentDate,  location = ph_location(left = 1, top =3, width = 4.25, height = 1))
  #doc <- ph_with(x = doc, value = CurrentDate, location = ph_location_type(type= "body") )
  
  #doc <- move_slide(doc, index = 1, to = 2)
  
  for (i in GeneName){
    Gene=i
    print(Gene)
  p <- PlotGTEXFunction(gene = Gene)
  p1 <- PlotFantonFunction(gene = Gene)
  p2 <- plotDataFunction(Gene=Gene)
  p3 <- plotImmuneProteomics(gene=Gene)
  p4 <- HPABloodPlot(gene=Gene)
  p5 <- SingleCelllot(gene=Gene)
  
  ggsave("Gtex.png", p[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Fantom.png", p1[[2]], device = "png", height = 2 , width = 5.65, scale=3)
  ggsave("Dice.png", p2, device = "png", height = 1.75 , width = 4.25, scale=3)
  ggsave("Protein.png", p3, device = "png", height = 1.75 , width =4.25, scale=3)
  ggsave("HPA.png", p4, device = "png", height = 1.75 , width =4.25, scale=3)
  ggsave("Single.png", p5, device = "png", height = 1.75 , width =4.25, scale=3)
  
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
    style = fp_text(color = "red", font.size = 20))
  small_text <- fp_text(color = "red", font.size = 18)
  smaller_text <- fp_text(color = "red", font.size = 12)
  doc <- add_slide(doc, layout = "1_One Column Text", master = "II2021")
  doc <- ph_with(x = doc, value = paste0("Body Expression ", Gene), location = ph_location_type(type = "title") )
  
  
  doc <- ph_with(x = doc, value = external_img("Gtex.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 0.92, width = 7, height = 2.5 ))
  doc <- ph_with(x = doc, value = external_img("Fantom.png",width = 5.65, height = 2 ), location = ph_location(left = 0.47, top = 3.09, width = 7, height = 2.5))
  doc <- ph_with(x = doc, value = external_img("Dice.png"), location = ph_location(left = 7.86, top = .31, width = 4, height = 2. ))
  doc <- ph_with(x = doc, value = external_img("Protein.png" ), location = ph_location(left = 7.86, top = 2.38, width = 4, height = 2.))
  doc <- ph_with(x = doc, value = external_img("HPA.png" ), location = ph_location(left = 7.86, top =4.25, width = 4, height = 2.))
  doc <- ph_with(x = doc, value = external_img("Single.png" ), location = ph_location(left = 0.47, top =4.73, width = 7, height = 2.5))
  doc <- ph_with(x = doc, value =  fp, location = ph_location(left=7.63, top=6.19, width = 6.83, height = .5)) 
  
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
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 11.88, top = 4.37, width = 0.73, height = 0.4 ))
  
  ul$str <- "Protein"
  doc <- ph_with(x = doc, value = ul, location = ph_location(left = 11.88, top = 2.26, width = 0.98, height = .59 ))
  
  ul <- unordered_list(
    level_list = c(2),
    str_list = c("Level1"),
    style = fp_text(color = "black", font.size = 12, font.family = "Arial", bold=FALSE))
  ul$str="https://dice-database.org/"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 1.24, width = 1.5, height = .5))
  ul$str="Nature Immunology vol 18, pages583–593(2017)"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 2.78, width =1.5, height = .83))
  ul$str="https://www.proteinatlas.org/"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 11.88, top = 4.77, width =1.5, height = .45))
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
  doc <- add_slide(doc, layout = "1_One Column Text", master = "II2021")
  doc <- ph_with(x = doc, value = paste0("Disease Expression ", Gene), location = ph_location_type(type = "title") )
  doc <- ph_with(x = doc, value = external_img("MetaTissue.png",width = 6.25, height = 5.56 ), location = ph_location(left = 0, top = 1.32, width = 6.25, height = 5.56 ))
  doc <- ph_with(x = doc, value = external_img("MetaBlood.png",width = 6.25, height = 5.56 ), location = ph_location(left = 6.61, top = 1.32, width = 6.25, height = 5.56 ))
  doc <- ph_with(x = doc, value = c("Tissue"), location = ph_location(left=2.65, top=.91, width = 4.52, height = .5)) 
  doc <- ph_with(x = doc, value = c("Blood"), location = ph_location(left=9.17, top=0.91, width = 4.52, height = .5)) 
  doc <- ph_with(x = doc, value =  fp, location = ph_location(left=4.45, top=6.0, width = 4.17, height = 1.45)) 
  ul$str="Meta analysis is a summary of multiple independent studies"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 6.08, top = 6.91, width =4.11, height = .71))
  ul$str="Matt Maciejewski"
  doc <- ph_with(x = doc, value = ul,
                 location = ph_location(left = 6.08, top = 6.91, width =2.2, height = 0.04))
}  
    
if(RA==TRUE){      
  doc <- add_slide(doc, layout = "1_One Column Text", master = "II2021")
  #Add 4 plots
    doc <- ph_with(x = doc, value = paste0("Disease Expression ", Gene), location = ph_location_type(type = "title"))
      doc <- ph_with(x = doc, value = external_img("AMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top, width = 5.65, height = 2 ))
      doc <- ph_with(x = doc, value = external_img("LOWAMP.png",width = 5.65, height = 2 ), location = ph_location(left = 0.5, top = top + 2, width = 5.65, height = 2))
      doc <- ph_with(x = doc, value = external_img("Synovial.png"), location = ph_location(left = 7.64, top = top, width = 4.25, height = 2 ))
      doc <- ph_with(x = doc, value = external_img("Peak.png"), location = ph_location(left = 7.64, top = top+2, width = 4.25, height = 2 ))
      #Add AMP Logo
      doc <- ph_with(x = doc, value = external_img("AMPLogo.png", width = 0.92, height =0.78), location = ph_location(left = 0, top =5.59, width = 0.92, height =0.78))
      #Add text boxes
      # first define a fpar ----
      fp1 <- unordered_list(
        level_list = c(2),
        str_list = c("Level2"),
        style = fp_text(color = "red", font.size = 18))
      fp2 <- unordered_list(
        level_list = c(2,2),
        str_list = c("Level2", "Level2"),
        style = fp_text(color = "red", font.size = 18))
      fp2$str <- c("Based on AMP Data", "...between OA and RA or between healthy and RA")
      doc <- ph_with(x = doc, value = fp2, 
                     location = ph_location(left=0.81, top=5.56, width = 6.28, height = 1.32)) 
      fp1$str <- c("Public Synovial data shows ....")
      doc <- ph_with(x = doc, value = fp1, 
                     location = ph_location(left=7.27, top=3.07, width = 4.76, height = 0.71)) 
      fp1$str <- c("Expression in PBMC vs Synovial .....")
      doc <- ph_with(x = doc, value = fp1, 
                     location = ph_location(left=7.27, top=5.91, width = 5.36, height = 0.71)) 
      #Add References
      ul$str="PLoS One 2017 Sep 1;12(9): e0183928."
      doc <- ph_with(x = doc, value = ul,
                     location = ph_location(left = 11.88, top = 1.73, width =1.15, height = 0.86))
      ul$str="Cell Rep. 2019 Aug 27; 28(9): 2455–2470.e5"
      doc <- ph_with(x = doc, value = ul,
                     location = ph_location(left = 11.88, top = 4.16, width =1.15, height = 0.86))
    #doc <- ph_with(x = doc, value =  fp, location = ph_location(left=6.91, top=4.97, width = 6.83, height = .5)) 
    
  doc <- add_slide(doc, layout = "1_One Column Text", master = "II2021")
    doc <- ph_with(x = doc, value = "Disease Expression", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("GSE55235.png" ), location = ph_location(left = 7.08, top =top+2, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = external_img("Peak.png"), location = ph_location(left = 0.5, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("SDY1299.png" ), location = ph_location(left = 0.5, top =top+2, width = 4.25, height = 2))
    #doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=5.83, width = 4.52, height = .5)) 
    doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 

    #ADD SLIDE FOR TREATMENT EFFECTS      
  doc <- add_slide(doc, layout = "1_One Column Text", master = "II2021")
    doc <- ph_with(x = doc, value = "Treatment Effects", location = ph_location_type(type = "title") )
    doc <- ph_with(x = doc, value = external_img("Recombine.png",width = 5.65, height = 2 ), 
                   location = ph_location(left = 0.5, top = top, width = 4.25, height = 2 ))
    doc <- ph_with(x = doc, value = external_img("GSE24742.png",width = 5.65, height = 2 ), 
                   location = ph_location(left = 4.75, top = top, width = 4.25, height = 2))
    doc <- ph_with(x = doc, value = external_img("GSE45867.png"), 
                   location = ph_location(left = 8.82, top = top, width = 4.25, height = 2 ))
    #Add References
    ul$str="KI collaboration-Synovial Tissue"
    doc <- ph_with(x = doc, value = ul,
                   location = ph_location(left = 0.76, top = 3.66, width =2.98, height = 0.34))
    ul$str="Arthritis Rheum 2011 May;63(5):1246-54"
    doc <- ph_with(x = doc, value = ul,
                   location = ph_location(left = 4.75, top = 3.66, width =3.47, height = 0.33))
    ul$str = "Arthritis Rheumatol 2014 Jan;66(1):15-23"
    doc <- ph_with(x = doc, value = ul,
                   location = ph_location(left = 9.2, top = 3.66, width =4.36, height = 0.33))
    
   # doc <- ph_with(x = doc, value = c("Un titre", "Deux titre"), location = ph_location(left=3.63, top=6.07, width = 4.52, height = .5)) 
    doc <- ph_with(x = doc, value =  fp, location = ph_location(left=0.64, top=6.0, width = 6.83, height = .5)) 

    
    
    }#End if RA==TRUE  
  }# end i in genename    
  #doc <- ph_with(x = doc, value = iris[1:4, 3:5], location = ph_location_right() )
    #   slide_summary(doc)
  doc <- move_slide(doc, index = 1, to = 2)
    target <- paste0("PPT/RAExpression_", Gene, ".pptx")
  print(doc, target = target)
}


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
PPTCreate( list("ERN1", "ERN2"), Index = TRUE, RA = FALSE)
PPTCreate( list("CH25H"), Index = TRUE, RA = FALSE)
PPTCreate( list("IL13","IL17A", "IL12B"), Index = TRUE, RA = FALSE)
PPTCreate( list("CXCL16","CXCR6", "CX3CL1", "CX3CR1"), Index = TRUE, RA = TRUE)
PPTCreate( list("PIM1","PIM3"), Index = TRUE, RA = TRUE)
PPTCreate( list("PKM2"), Index = TRUE, RA = TRUE)
PPTCreate(GeneName = list("SIK1", 'SIK2', "SIK3", "PFKFB3", "CTCS", "ITK", "TMEM173", "FAP")) 
PPTCreate( list("DPEP1"), Index = TRUE, RA = TRUE)
PPTCreate( list("NOX1","CYBB","NOX3","NOX4"), Index = TRUE, RA = FALSE)
PPTCreate( list("JAK1", "JAK3", "IL6R"), Index = TRUE, RA = TRUE)
PPTCreate( list("FURIN"), Index = TRUE, RA = TRUE)
PPTCreate( list("LAG3", "HMGB1", "CD209", "BTLA", "TIGIT", "CD226", "PVRIG", "PVRL2", "LAIR2", 
                "LAIR1","PCBP1", "ILDR2"), Index = TRUE, RA = TRUE)
PPTCreate( list("TNFRSF18", "TNFSF18", "TNFRSF1B"), Index = TRUE, RA = TRUE)
PPTCreate( list("GSDMD"), Index = TRUE, RA = TRUE)
PPTCreate( list("HERC5"), Index = TRUE, RA = TRUE)
PPTCreate( list("GSDMA", "GSDMB", "GSDMC", "GSDMD", "GSDME"), Index = TRUE, RA = FALSE)
PPTCreate( list("IRAK2"), Index = TRUE, RA = TRUE)
PPTCreate( list("CMKLR1"), Index = TRUE, RA = TRUE)
PPTCreate( list("HGF", "MET"), Index = TRUE, RA = TRUE)
PPTCreate( list("CD74", "SPPL2A"), Index = TRUE, RA = TRUE)
PPTCreate( list("MOSPD2"), Index = TRUE, RA = TRUE)
PPTCreate( list("ANXA5"), Index = TRUE, RA = TRUE)
PPTCreate( list("SLC40A1","HAMP", "FTH1", "FTL"), Index = TRUE, RA = TRUE)
PPTCreate( list("CD109"), Index = TRUE, RA = TRUE)
PPTCreate( list("IL4", "IL13", "PDCD1", "PDCD2"), Index = TRUE, RA = TRUE)
PPTCreate( list("EPAS1", "SLC11A2", "CYBRD1", "SLC40A1", "GAPDH"), Index = TRUE, RA = FALSE)
PPTCreate( list("CD300A","CD300C", "CD300E", "CD300LB", "CD300LD", "CD300LF", "CD300LG"), Index = FALSE, RA = FALSE)
PPTCreate( list("CD200E"), Index = TRUE, RA = TRUE)

PPTCreate( list("CCL18", "PITPNM3", "CCR8"), Index = TRUE, RA = TRUE)
PPTCreate( list("ANXA1", "FPR2", "FPR1"), Index = TRUE, RA = TRUE)
PPTCreate( list("MLKL"), Index = TRUE, RA = TRUE)
PPTCreate( list("PTPN11"), Index = TRUE, RA = TRUE)
PPTCreate( list("BTK", "BMX", "ITK", "TXK", "TEC"), Index = TRUE, RA = TRUE)
  options(error=NULL)
  options(error=recover)
  PPTCreate( list("CXCR2"), Index = TRUE, RA = TRUE)
  PPTCreate( list("CD63", "GPNMB"), Index = TRUE, RA = TRUE)
  PPTCreate( list("FLT4", "VEGFC", "NRP2", "PROX1", "LYVE1", "CCBE1", "ADAMTS3", "CTSD"), Index = TRUE, RA = TRUE)
  
  
  PPTCreate( list("P4HA2"), Index = TRUE, RA = TRUE)
  PPTCreate( list("MRGPRX2"), Index = TRUE, RA = TRUE)
  PPTCreate( list("IL32"), Index = TRUE, RA = TRUE)
  PPTCreate( list("SH2D1B", "SLAMF7", "PTPN6", "SH2D1A"), Index = TRUE, RA = TRUE)
  PPTCreate( list("PRKAB1"), Index = TRUE, RA = FALSE)
  PPTCreate( list("TCF21"), Index = TRUE, RA = TRUE)
  
  PPTCreate( list("PDGFRA", "PDGFRB", "PDGF", "PDGFB", "PDGFC", "PDGFD"), Index = TRUE, RA = FALSE)
  PPTCreate( list("PDGFA"), Index = TRUE, RA = FALSE)
  
  PPTCreate( list("CD84"), Index = TRUE, RA = TRUE)
  
  PPTCreate( list("IGF1", "IGF1R", "IGFBP3","TMEM219"), Index = TRUE, RA = TRUE)
             PPTCreate( list("IL27RA", "IL6ST"), Index = TRUE, RA = FALSE)
  
             PPTCreate( list("LILRB4"), Index = TRUE, RA = TRUE)
             
             
             PPTCreate( list("USP18"), Index = TRUE, RA = TRUE) 
             
             PPTCreate( list("LTB4R2"), Index = TRUE, RA = TRUE)  
             