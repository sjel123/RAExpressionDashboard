library(Biobase)

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