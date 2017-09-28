library(shinydashboard)
library(shiny)
library(plotly)
library(hexbin)
library(htmlwidgets)
library(dplyr)
library(tidyr)
library(data.table)
library(DESeq2)

dat <- readRDS("RLogDat.Rds")
datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]

myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
metrics <- readRDS("topGenes_limma.Rds")
for (i in 1:(length(myPairs)-1)){
  for (j in (i+1):length(myPairs)){
    setDT(metrics[[paste0(myPairs[i],"vs",myPairs[j])]], keep.rownames = TRUE)[]
    colnames(metrics[[paste0(myPairs[i],"vs",myPairs[j])]])[1] <- "ID"
    metrics[[paste0(myPairs[i],"vs",myPairs[j])]]$ID <- as.factor(metrics[[paste0(myPairs[i],"vs",myPairs[j])]]$ID)
    metrics[[paste0(myPairs[i],"vs",myPairs[j])]] <- as.data.frame(metrics[[paste0(myPairs[i],"vs",myPairs[j])]])
  }
}
myMetrics <- colnames(metrics[[1]])[-which(colnames(metrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

shinyServer(function(input, output, session){
  observeEvent(c(input$goButton1, input$goButton2), values$x <- values$x + 1)
  observeEvent(c(input$selPair1, input$selPair2), values$x <- 0)
  observeEvent(c(input$selMetric1, input$selMetric2), values$x <- 0)
  observeEvent(c(input$selOrder1, input$selOrder2), values$x <- 0)
  observeEvent(input$binSize, values$x <- 0)
  
  observeEvent(input$selPair1, values$selPair <- input$selPair1)
  observeEvent(input$selPair2, values$selPair <- input$selPair2)
  observeEvent(input$selMetric1, values$selMetric <- input$selMetric1)
  observeEvent(input$selMetric2, values$selMetric <- input$selMetric2)  
  observeEvent(input$selOrder1, values$selOrder <- input$selOrder1)
  observeEvent(input$selOrder2, values$selOrder <- input$selOrder2)
  
  observe({x <- values$selPair
  updateSelectizeInput(session, "selPair1", "Treatment pairs:", choices = myPairs, options = list(maxItems = 2), selected = x)
  updateSelectizeInput(session, "selPair2", "Treatment pairs:", choices = myPairs, options = list(maxItems = 2), selected = x)})
  
  observe({x <- values$selMetric
  updateSelectizeInput(session, "selMetric1", "Metric:", choices = myMetrics, selected = x)
  updateSelectizeInput(session, "selMetric2", "Metric:", choices = myMetrics, selected = x)})
  
  observe({x <- values$selOrder
  updateSelectizeInput(session, "selOrder1", "Metric order:", choices = c("Increasing", "Decreasing"), selected = x)
  updateSelectizeInput(session, "selOrder2", "Metric order:", choices = c("Increasing", "Decreasing"), selected = x)}) 
  
  # Create data subset based on two letters user chooses
  datSel <- eventReactive(values$selPair, {
    validate(need(length(values$selPair) == 2, "Select a pair of treatments."))
    sampleIndex <- reactive(which(sapply(colnames(dat), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1], values$selPair[2])))
    dat[,c(1, sampleIndex())]
  }, ignoreNULL = FALSE)
  
  metricDF <- eventReactive(c(values$selPair, values$selMetric, values$selOrder), {
    metricDF <- metrics[[paste0(values$selPair[1], "vs", values$selPair[2])]]
    if (!is.null(metricDF[[values$selMetric]])){
      metricDF <- metricDF[order(metricDF[[values$selMetric]]),]
      if (values$selOrder == "Decreasing"){
        metricDF <- metricDF[order(-metricDF[[values$selMetric]]),]
      }
    }
    metricDF
  })
  
  currMetric <- eventReactive(values$x, {
    validate(need(values$x > 0, "Plot a case."))
    metricDF()[values$x, ]})
  currID <- eventReactive(currMetric(), {as.character(currMetric()$ID)})
  currGene <- eventReactive(currID(), {unname(unlist(datSel()[which(datSel()$ID == currID()), -1]))})
  output$info1 <- renderPrint({ currMetric() })
  output$info2 <- renderPrint({ currMetric() })
  
  output$hexPlot <- renderPlotly({
    
    sampleIndex1 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1]))
    sampleIndex2 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[2]))
    
    minVal = min(datSel()[,-1])
    maxVal = max(datSel()[,-1])
    maxRange = c(minVal, maxVal)
    xbins= input$binSize
    buffer = (maxRange[2]-maxRange[1])/(xbins/2)
    x <- c()
    y <- c()
    for (i in 1:length(sampleIndex1)){
      for (j in 1:length(sampleIndex2)){
        x <- c(x, unlist(datSel()[,(sampleIndex1[i])]))
        y <- c(y, unlist(datSel()[,(sampleIndex2[j])]))
      }
    }
    
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    
    my_breaks <- c(2, 4, 6, 8, 20, 1000)    
    p <- reactive(ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = values$selPair[1], y = values$selPair[2]) + coord_fixed(xlim = c(-0.5, (maxRange[2]+buffer)), ylim = c(-0.5, (maxRange[2]+buffer))) + theme(aspect.ratio=1) + scale_fill_gradient(name = "count", trans = "log", breaks = my_breaks, labels = my_breaks, guide="legend"))
    
    plotlyHex <- reactive(ggplotly(p(), height = 400, width = 400) %>% config(displayModeBar = F))
    
    # Use onRender() function to draw x and y values of selected row as orange point
    plotlyHex() %>% onRender("
     function(el, x, data) {
     noPoint = x.data.length;
     Shiny.addCustomMessageHandler('points', function(drawPoints) {
     if (x.data.length > noPoint){
     Plotly.deleteTraces(el.id, x.data.length-1);
     }
     var Traces = [];
     var trace = {
     x: drawPoints.geneX,
     y: drawPoints.geneY,
     mode: 'markers',
     marker: {
     color: 'orange',
     size: drawPoints.pointSize
     },
     hoverinfo: 'none'
     };
     Traces.push(trace);
     Plotly.addTraces(el.id, Traces);
     });}")
  })
  
  observe({
    # Get x and y values of selected row
    sampleIndex1 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1]))
    sampleIndex2 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[2]))
    
    geneX <- c()
    geneY <- c()
    for (i in 1:length(sampleIndex1)){
      for (j in 1:length(sampleIndex2)){
        geneX <- c(geneX, currGene()[sampleIndex1[i]-1])
        geneY <- c(geneY, currGene()[sampleIndex2[j]-1])
      }
    }
    
    pointSize <- input$pointSize
    
    # Send x and y values of selected row into onRender() function
    session$sendCustomMessage(type = "points", message=list(geneX=geneX, geneY=geneY, pointSize = pointSize))
  })
  
  output$boxPlot <- renderPlotly({
    nVar = reactive(ncol(datSel()))
    colNms <- reactive(colnames(datSel()[, c(2:nVar())]))
    
    boxDat <- eventReactive(datSel(), {
      boxDat <- datSel()[, c(1:nVar())] %>% gather(key, val, -c(ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      boxDat
    })
    
    BP <- reactive(ggplot(boxDat(), aes(x = Sample, y = Count)) + geom_boxplot())
    ggBP <- reactive(ggplotly(BP(), width=600, height = 400) %>% config(displayModeBar = F, staticPlot = T))
    
    observe({
      session$sendCustomMessage(type = "lines", currGene())
    })
    
    ggBP() %>% onRender("
      function(el, x, data) {
      
      console.log(['x', x])
      console.log(['x.data.length', x.data.length])
      noPoint = x.data.length;
      
      function range(start, stop, step){
      var a=[start], b=start;
      while(b<stop){b+=step;a.push(b)}
      return a;
      };
      
      Shiny.addCustomMessageHandler('lines',
      function(drawLines) {
      
      i = x.data.length
      if (i > 1){
      while (i > 1){
      Plotly.deleteTraces(el.id, (i-1));
      i--;
      }
      }
      
      var dLength = drawLines.length
      
      var Traces = [];
      var traceLine = {
      x: range(1, dLength, 1),
      y: drawLines,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1
      },
      opacity: 0.9,
      hoverinfo: 'none'
      }
      Traces.push(traceLine);
      Plotly.addTraces(el.id, Traces);
      })
      }")})
})
