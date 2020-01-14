#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


## app.R ##
library(shiny)
library(shinyTree)
library(shinydashboard)
library(ggpubr)
library(htmlwidgets)

source("utils.R")
source("HPAplots.R")
# res <- Gtex[,1:2]
# res <- res[(!duplicated(res$Description)),]
# colnames(res) <- c("Gene_Name", "Gene_id")
load("/app/Shiny/RAExpressionDashboard/Data/res1.rData")

js <- HTML("
$(document).on('shiny:connected', function(event) {
  $('#dataset1').on('hover_node.jstree ', function (e, data) {
    var nodetext = data.node.text;
     Shiny.setInputValue('hoverednode', nodetext + ' is hovered on');
  });
});
")


ui <- dashboardPage(
        title="RA Dashboard",
          header=dashboardHeader(
            tags$li(class = "dropdown",
                    (tags$img(height = "200px", alt="SNAP Logo", src="CSI_3.jpg",style = "padding-top:0px; padding-bottom:50px;")
                    )),
              titleWidth = '100%',

                   title = span(
                     tags$img(src="CSI_3.jpg", hieght = "200px", style="repeat-y"),
                     column(12, class="title-box",
                            tags$h1(class="primary-title", style='margin-top:10px;', 'Rheumatoid Arthritis Synovial'),
                            tags$h1(class="primary-subtitle", style='margin-top:10px;', 'Tissue Expression')
                     )#end Column row
                   ), #End Tit;e
                  dropdownMenuOutput("helpMenu")
                  ),#End dashboard header
                body = dashboardBody(
                  tags$style(type="text/css", "
                /*    Move everything below the header */
                             .content-wrapper {
                             margin-top: 50px;
                             }
                             .content {
                             padding-top: 60px;
                             }
                             /*    Format the title/subtitle text */
                             .title-box {
                             position: absolute;
                             text-align: center;
                             top: 50%;
                             left: 50%;
                             transform:translate(-50%, -50%);
                             }
                             @media (max-width: 590px) {
                             .title-box {
                             position: absolute;
                             text-align: center;
                             top: 10%;
                             left: 10%;
                             transform:translate(-5%, -5%);
                             }
                             }
                             @media (max-width: 767px) {
                             .primary-title {
                             font-size: 1.1em;
                             }
                             .primary-subtitle {
                             font-size: 1em;
                             }
                             }
                             /*    Make the image taller */
                             .main-header .logo {
                             height: 125px;
                             }
                             /*    Override the default media-specific settings */
                             @media (max-width: 50px) {
                             .main-header {
                             padding: 0 0;
                             position: relative;
                             }
                             .main-header .logo,
                             .main-header .navbar {
                             width: 100%;
                             float: none;
                             }
                             .main-header .navbar {
                             margin: 0;
                             }
                             .main-header .navbar-custom-menu {
                             float: right;
                             }
                             }
                             /*    Move the sidebar down */
                             .main-sidebar {
                             position: absolute;
                             }
                             .left-side, .main-sidebar {
                             padding-top: 175px;
                             }"
                ),#end tag style
                tabItems(
                  # First tab content
                  tabItem(tabName = "dashboard",
                          fluidRow(
                            box(DT::dataTableOutput('x1'),  hr(),fluid=FALSE, width = 12)
                          ),#End Fluid Row
                          # fluidRow(column("",
                          #                 tableOutput("dataset"), width = 12,
                          #                 align = "center")),
                          # Button
                          downloadButton("downloadData", "Download Rainbow Table"),
                          fluidRow(
                            uiOutput("boxes")
                          )#End Fluid Row
                  ), #End TabItem
                  
                  # Second tab content
                  tabItem(tabName = "GraphsOptions",
                          h2("Select Dataset to show")
                  )#End TabItem
                )#End TabItems
                ),#end dashboardBody
  
  
  dashboardSidebar(sidebarMenu(
    menuItem("Dashboard1", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Datasets", tabName = "GraphsOptions", icon = icon("th-list"), startExpanded = TRUE),
  
    
    h2(tags$div("Select Datasets", tags$br(), "to Visualize")),
    tags$head(tags$script(js)),
    shinyTree("dataset1", checkbox = TRUE, wholerow = FALSE),
    br(),
    verbatimTextOutput("hoverinfo"),
    
 
  checkboxGroupInput("Columns", "Select Columns to Filter in Table:",
                     names(res.rainbow)[4:24])
 # # Button
 # downloadButton("downloadData", "Download")                   
   ))#End dashboardSudeBar
)

server <- function(input, output) { 


  res <- res.rainbow
  
 DF1 <- reactive({
   res1 <- res
    index <- which(names(res) %in% input$Columns)
          index2=0
      for(i in index){
          index2 <- which(res1[,i]==1)
          res1 <- res1[index2,]
      }
          res1
 })      
   
 
 output$dataset1 <- renderTree({ 
   sss= list( 
     'Cell/Tissue Expression'   =  structure(list("Fantom5" = "f5",
                                                    "GTEX" = "gtex"),stselected=FALSE),  
     'Immune Cell Expression' = structure(list("ImmuneNavigator" = "IN",
                                               "DICE" ="dice",
                                               "ImmuneProt"="ImmuneProt",
                                               "HPA"="HPA"),
                                                stselected=FALSE),
     'Fibroblast Expression'     =  structure(list('AMP'='AMP', 
                                                   "AMP2"= "AMP2",
                                                   "AMPSC" ="AMPSC",
                                                   "BrennerFib" = "BrennerFib",
                                                   "BrennerFibMicro" = "BrennerFibMicro",
                                                   "BrennerSingleCell" = "BrennerSingleCell",
                                                   "FibroblastSignature"= "FibroblastSignature"
                                                   ), stselected=FALSE),
                                                
     'Disease Expression'     =  structure(list('AMP'='AMP', 
                                                "AMP2"= "AMP2",
                                                "AMPSC" ="AMPSC",
                                                "PEAC" = "peac",
                                                "DiseaseProgression" = "pub",
                                                 "GSE55235"= "GSE55235",
                                                 "GSE55457" = "GSE55457",
                                                 "GSE1919" = "GSE1919",
                                                 "GSE7729" = "GSE7729",
                                                 "SDY1299" = "SDY1299",
                                                 "GSE36700" = "GSE36700",
                                                 "GSE116899" = "GSE116899"), stselected=FALSE),
     'Treatment Effects'      = structure(list("Recombine" ="recombine",
                                                "GSE24742" = "GSE24742",
                                                "GSE45867" = "GSE45867"), stselected=FALSE)
     ) 
 })

 output$hoverinfo <- renderText({
   req(input$hoverednode)
   input$hoverednode
 })
 
 output$dataset <- renderTable(
   names(as.data.frame(get_selected(input$dataset1, format = "slices"))))
 # names(as.data.frame(get_tree_value(input$dataset1, format = "slices"))))

 dataset <- reactive({
   names(as.data.frame(get_selected(input$dataset1, format = "slices")))
 })
 
 
 #editing the tbl dataset with the filtered rows only
 thedata_filtered <- reactive({
   DF1()[c(input[["x1_rows_all"]]), ]
 })
 
 # Downloadable csv of selected dataset ----
 # output$downloadData <- downloadHandler(
 #   filename = function() {
 #     paste("RainbowTable_filt", ".csv", sep = "")
 #   },
 #   content = function(file) {
 #     write.csv(thedata_filtered(), file, row.names = FALSE)
 #   }
 # )
  
  #output$results <- input$data.set
  
  output$boxes <- renderUI({
    tagList(
      if("Fantom5" %in% dataset()){
        box(title = "Fantom5 Data", status = "primary", solidHeader = TRUE,
            plotOutput("plot1"), width=12)
      },
      if("GTEX" %in% dataset()){
        box(title = "GTEX Data",status = "primary", solidHeader = TRUE,
            plotOutput("plot2"), width=12)
      },
      if("ImmuneNavigator" %in% dataset()){
        box(title = "ImmuneNavigator",status = "primary", solidHeader = TRUE,
            plotOutput("plot3"), width=12)
      },
      if("DICE" %in% dataset()){
        box(title = "DICE2",status = "primary", solidHeader = TRUE,
            plotOutput("plot4"), width=12)
      },
      if("ImmuneProt" %in% dataset()){
        box(title = "Immune Cell Proteomics",status = "primary", solidHeader = TRUE,
            plotOutput("plot4a"), width=12)
      },
      if("HPA" %in% dataset()){
        box(title = "HPA Immune Cell Expression",status = "primary", solidHeader = TRUE,
            plotOutput("plotHPA3"),width=12)
      },
      if("HPA" %in% dataset()){
        box(title = "HPA SchmiedelBloodPlot Immune Cell Expression",status = "primary", solidHeader = TRUE,
            plotOutput("plotHPA1"),width=12)
      },
      if("HPA" %in% dataset()){
        box(title = "Monaco HPA Immune Cell Expression",status = "primary", solidHeader = TRUE,
            plotOutput("plotHPA2"),width=12)
      },

      if("AMP" %in% dataset()){
        box(title = "AMPLowInput",status = "primary", solidHeader = TRUE,
            plotOutput("plot5"), width=12)
      },
      if("AMP2" %in% dataset()){
        box(title = "AMPLeukocyute",status = "primary", solidHeader = TRUE,
            plotOutput("plot51"), width=12)
      },
      if("AMPSC" %in% dataset()){
        box(title = "AMPSingle Cell Mean",status = "primary", solidHeader = TRUE,
            plotOutput("plotAMPSC"), width=12)
      },
  
      if("BrennerFib" %in% dataset()){
        box(title = "BrennerFib",status = "primary", solidHeader = TRUE,
            plotOutput("plot52"), width=12)
      },
      if("BrennerFibMicro" %in% dataset()){
        box(title = "BrennerFibMicro",status = "primary", solidHeader = TRUE,
            plotOutput("plot53"), width=12)
      },
      
       if("BrennerSingleCell" %in% dataset()){
          box(title = "BrennerFibMicro",status = "primary", solidHeader = TRUE,
              plotOutput("plotSSRA"), width=12)        
      },
      if("FibroblastSignature" %in% dataset()){
        box(title = " FibroblastSignature",status = "primary", solidHeader = TRUE,
            plotOutput("plotFibSig"), width=12)        
      },
     
      if("PEAC" %in% dataset()){
        box(title = "PEAC Cohort",status = "primary", solidHeader = TRUE,
            plotOutput("plot6"), width=12)
      },
      if("DiseaseProgression" %in% dataset()){
        box(title = "Public Synovial",status = "primary", solidHeader = TRUE,
            plotOutput("plot7"), width=12)
      },
      if("Recombine" %in% dataset()){
        box(title = "Recombine Synovial",status = "primary", solidHeader = TRUE,
            plotOutput("plot8"), width=12)
      },
      if("GSE55235" %in% dataset()){
        box(title = "GSE55235",status = "primary", solidHeader = TRUE,
            plotOutput("plot9"), width=12)
      },
      if("GSE55457" %in% dataset()){
        box(title = "GSE55457",status = "primary", solidHeader = TRUE,
            plotOutput("plot10"), width=12)
      },
      if("GSE1919" %in% dataset()){
        box(title = "GSE1919",status = "primary", solidHeader = TRUE,
            plotOutput("plot11"), width=12)
      },
      if("GSE7729" %in% dataset()){
        box(title = "GSE7729",status = "primary", solidHeader = TRUE,
            plotOutput("plot12"), width=12)
      },
      if("SDY1299" %in% dataset()){
        box(title = "SDY1299",status = "primary", solidHeader = TRUE,
            plotOutput("plot13"), width=12)
      },
      if("GSE24742" %in% dataset()){
        box(title = "GSE24742",status = "primary", solidHeader = TRUE,
            plotOutput("plot14"), width=12)
      },
      if("GSE45867" %in% dataset()){
        box(title = "GSE45867",status = "primary", solidHeader = TRUE,
            plotOutput("plot15"), width=12)
      },
      if("GSE36700" %in% dataset()){
        box(title = "GSE36700",status = "primary", solidHeader = TRUE,
            plotOutput("plot16"), width=12)
      },
      if("GSE116899" %in% dataset()){
        box(title = "GSE116899",status = "primary", solidHeader = TRUE,
            plotOutput("plot17"), width=12)
      }
      


    )#End List 
  })#End Output$box
  

  
  output$x1 = DT::renderDataTable(DF1(), server = T, escape=T, selection = 'single',
                                  extensions = "FixedColumns",
                                  filter = 'top',
                                  options=list(
                                        lengthMenu = list(c(5, 10, 15, -1), c('5', '10', '15', 'All')),
                                        autoWidth = TRUE,
                                        columnDefs = list(list(width = '20px', targets = c(3:7))),
                                        scrollX = TRUE,
                                        fixedColumns = list(leftColumns = 2)
    )
    )
  
  SW <- reactive({ 
    s = NULL
    #s = res$Gene_Name[input$x1_rows_selected]
    s = DF1()
    s = s$Gene_Name[input$x1_rows_selected]
    print(sprintf("input$x1_rows_selected %s", input$x1_rows_selected))
    if (is.null(s)) s="THY1"  
    if (length(s)==0) s="THY1"  
    print(sprintf("row selected %s", s))
    s
  })
  
  output$plot1 <- renderPlot({
    #data <- histdata[seq_len(input$slider)]
    PlotFantonFunction(gene=SW())
  })
  output$plot2 <- renderPlot({
    #data <- histdata[seq_len(input$slider)]
    PlotGTEXFunction(gene=SW())
  })
  output$plot3 <- renderPlot({
    ImmNavPlot(SW())
  })
  #DICE Data
  output$plot4 <- renderPlot({
    plotDataFunction(SW())
  })
  #Immune Proteomics Data
  output$plot4a <- renderPlot({
    plotImmuneProteomics(SW())
  })
  output$plot5 <- renderPlot({
    PlotLowInputAMP(SW())
  })
  output$plot51 <- renderPlot({
    PlotAMP(gene=SW())
  })
  output$plotAMPSC <- renderPlot({
    PlotAMPSingleCell(gene=SW())
  })

  output$plot52 <- renderPlot({
    PlotBrennerFibro (gene=SW())
  })
  
  output$plot53 <- renderPlot({
    PlotBrennerFibroMicro (gene=SW())
  })
  output$plotSSRA <- renderPlot({
    PlotSSRA (gene=SW())
  })
  output$plotFibSig <- renderPlot({
    PlotFIBData (gene=SW())
  })
  output$plot6 <- renderPlot({
    PlotPeak(gene=SW())
  }) 
  output$plot7 <- renderPlot({
    PlotSynovial(gene=SW())
  }) 
  output$plot8 <- renderPlot({
    PlotRecombine(gene=SW())
  }) 
  output$plot9 <- renderPlot({
    PlotGSE55235(gene=SW())
  }) 
  output$plot10 <- renderPlot({
    PlotGSE55457(gene=SW())
  }) 
  output$plot11 <- renderPlot({
    PlotGSE1919(gene=SW())
  }) 
  output$plot12 <- renderPlot({
    PlotGSE7729(gene=SW())
  })
  output$plot13 <- renderPlot({
    PlotSDY1299(gene=SW())
  })
  output$plot14 <- renderPlot({
    PlotGSE24742(gene=SW())
  }) 
  output$plot15 <- renderPlot({
    PlotGSE45867(gene=SW())
  }) 
  output$plot16 <- renderPlot({
    PlotGSE36700(gene=SW())
  }) 
  output$plot17 <- renderPlot({
    PlotGSE116899(gene=SW())
  }) 
  output$plotHPA1 <- renderPlot({
    SchmiedelBloodPlot(gene=SW())
  }) 
  output$plotHPA2 <- renderPlot({
    MonacoBloodPlot(gene=SW())
  }) 
  output$plotHPA3 <- renderPlot({
    HPABloodPlot(gene=SW())
  }) 


  # # Downloadable csv of selected dataset ----
  # output$downloadData <- downloadHandlerdownloadHandler(
  #   filename <- function() {
  #     paste("output", "zip", sep=".")
  #   },
  #   
  #   content <- function(file) {
  #     file.copy("PPT/ph_with_demo.pptx", file)
  #   },
  #   contentType = "application/zip"
  # )
  
  }
shinyApp(ui, server)
