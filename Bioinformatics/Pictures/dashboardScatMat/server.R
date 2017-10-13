library(shinydashboard)
library(plotly)
library(ggplot2)
library(shiny)
library(htmlwidgets)
library(utils)
library(tidyr)
library(stats)
library(hexbin)
library(stringr)
library(dplyr)
library(data.table)
library(GGally)

xbins = 15

server <- function(input, output, session) {

  load("p.rda")
  load("pS.rda")
  load("ggPS.rda")
  load("data.rda")
  load("renderedPlot.rda")
  
  output$scatMatPlot <- renderPlotly({renderedPlot})
  selID <- reactive(input$selID)
  
  pcpDat <- reactive(data[which(data$ID %in% selID()), c(1:(p$nrow+1))])
  output$selectedValues <- renderPrint({pcpDat()$ID})
  colNms <- colnames(data[, c(2:(p$nrow+1))])
  
  boxDat <- data[, c(1:(p$nrow+1))] %>% gather(key, val, -c(ID))
  colnames(boxDat) <- c("ID", "Sample", "Count")
  BP <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot()
  ggBP <- ggplotly(BP, width=700)
  
  output$boxPlot <- renderPlotly({
    ggBP %>% onRender("
    function(el, x, data) {
    
    var Traces = [];
    
    var dLength = data.pcpDat.length
    var vLength = data.nVar
    var cNames = data.colNms
    
    for (a=0; a<dLength; a++){
    xArr = [];
    yArr = [];
    for (b=0; b<vLength; b++){
    xArr.push(b+1)
    yArr.push(data.pcpDat[a][cNames[b]]);
    }
    
    var traceHiLine = {
    x: xArr,
    y: yArr,
    mode: 'lines',
    line: {
    color: 'orange',
    width: 1.5
    },
    opacity: 0.9,
    }
    Traces.push(traceHiLine);
    }
    Plotly.addTraces(el.id, Traces);
    
    }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})
}