library(cowplot)
Genes <- c("CDK2", "CDK4", "CDK6")
textsize=10

plotfunctionname <- c("PlotFantonFunction","PlotGSE1919", "PlotGSE55235", 
                      "PlotGSE55457", "PlotGSE7729", "PlotGTEXFunction",
                      "PlotLowInputAMP", "PlotPeak", "PlotRecombine",
                      "PlotSDY1299", "PlotSynovial")      
     

SinglePlot <- function(x="PlotPeak"){
    Plot <- list()
      for (i in Genes){
        Plot[[i]] <- print(eval(call(x, i)))
        Plot[[i]] <- Plot[[i]] + theme(legend.position  = "none") +
          theme(plot.title = element_text(hjust = 0.5,size = 9)) +
          ggtitle("")
      }
    p <- plot_grid(plotlist=Plot,
              labels = Genes,#c('Fig A','Fig B','Fig C'),
              label_x = 0.1,
              ncol = 3)
    print(p)
    return(p)
}


SinglePlot("PlotPeak")
SinglePlot("PlotGSE55235")
SinglePlot("plotDataFunction")
SinglePlot("PlotLowInputAMP")
SinglePlot("PlotSynovial")
SinglePlot(plotfunctionname[4])
