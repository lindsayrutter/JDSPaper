
datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
metrics[[1]] <- metrics[[1]][which(metrics[[1]]$PValue<0.01),]
metrics[[1]] <- metrics[[1]][which(metrics[[1]]$ID %in% dat$ID),]
dat = dat[which(dat$ID %in% metrics[[1]]$ID),]
myMetrics <- colnames(metrics[[1]])[-which(colnames(metrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

sidebar <- dashboardSidebar(
  width = 180,
  hr(),
  sidebarMenu(id="tabs",
              menuItem("Binned scatterplot", tabName="hexPlot", selected=TRUE), # icon=icon("line-chart"),
              menuItem("Parallel coordinates", tabName = "boxPlot") # icon = icon("file-text-o")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "hexPlot",
            fluidRow(
              column(width = 4, 
                     box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
                         selectizeInput("selPair1", "Treatment pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
                         selectInput("selMetric1", "Metric:", choices = myMetrics),
                         selectInput("selOrder1", "Metric order:", choices = c("Increasing", "Decreasing")),
                         numericInput("binSize", "Hexagon size:", value = 10, min = 1),
                         numericInput("pointSize", "Point size:", value = 6, min = 1),
                         actionButton("goButton1", "Plot case!"))),
              column(width = 8,
                     box(width = NULL, plotlyOutput("hexPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),
            fluidRow(
              column(width = 8, offset = 4,
                     box(width = NULL, verbatimTextOutput("info1"), collapsible = TRUE, title = "Case metrics", status = "primary", solidHeader = TRUE)))),
    
    tabItem(tabName = "boxPlot",
            fluidRow(
              column(width = 4, 
                     box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
                         selectizeInput("selPair2", "Treatment pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
                         selectInput("selMetric2", "Metric:", choices = myMetrics),
                         selectInput("selOrder2", "Metric order:", choices = c("Increasing", "Decreasing")),
                         actionButton("goButton2", "Plot case!"))),
              column(width = 8,
                     box(width = NULL, plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))),
            fluidRow(
              column(width = 8, offset = 4,
                     box(width = NULL, verbatimTextOutput("info2"), collapsible = TRUE, title = "Case metrics", status = "primary", solidHeader = TRUE))))
  )
)

dashboardPage(
  dashboardHeader(title = "Overlaying cases", titleWidth = 180),
  sidebar,
  body
)

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
    
    # By default, groups into six equal-sized bins
    hexdf$countColor <- cut2(hexdf$counts, g=6)
    hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor), function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE), ",")[[1]][1], 2))))
    
    hexdf$countColor2 <- factor(hexdf$countColor2, levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))
    
    for (i in 1:(length(levels(hexdf$countColor2))-1)){
      levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],"-",levels(hexdf$countColor2)[i+1])
    }
    levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <- paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")
    
    #my_breaks = c(2, 4, 6, 8, 20, 1000)
    my_breaks = levels(hexdf$countColor2)
    clrs <- brewer.pal(length(my_breaks)+3, "Blues")
    clrs <- clrs[3:length(clrs)]
    
    print("my_breaks")
    print(str(my_breaks))
    print("clrs")
    print(str(clrs))
    print("hexdf$countColor2")
    print(str(hexdf$countColor2))
    print("hexdf")
    print(str(hexdf))
    print("HEAD")
    print(head(hexdf))
    
    p <- reactive(ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Cases count") + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = paste0("Read count ", "(", values$selPair[1], ")"), y = paste0("Read count ", "(", values$selPair[2], ")")) + coord_fixed(ratio=1))
    
    # coord_fixed(xlim = c(-0.5, (maxRange[2]+buffer)), ylim = c(-0.5, (maxRange[2]+buffer)))
    # theme(aspect.ratio=1)
    
    gP <- eventReactive(p(), {
      gP <- ggplotly(p(), height = 400) #  height = 400
      for (i in 1:(length(gP$x$data)-1)){
        info <- gP$x$data[i][[1]]$text
        info2 <- strsplit(info,"[<br/>]")
        myIndex <- which(startsWith(info2[[1]], "counts:"))
        gP$x$data[i][[1]]$text <- info2[[1]][myIndex]
      }
      gP$x$data[length(gP$x$data)][[1]]$text <- NULL
      gP
    })
    
    plotlyHex <- reactive(gP() %>% config(displayModeBar = F))
    
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
       hoverinfo: 'none',
       showlegend: false
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
    
    BP <- reactive(ggplot(boxDat(), aes(x = Sample, y = Count)) + geom_boxplot() + labs(y = "Read count"))
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

shinyApp(ui = ui, server = server)