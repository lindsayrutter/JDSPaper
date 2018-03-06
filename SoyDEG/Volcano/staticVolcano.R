library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
#library(GGally)
#library(edgeR)
library(plotly)
library(htmlwidgets)
library(shinyBS)
library(shiny)
library(DESeq2)

ui <- shinyUI(pageWithSidebar(
  headerPanel("Click the button"),
  sidebarPanel(
    uiOutput("selInput"), ########## choose all pairs
    uiOutput("slider"),
    sliderInput("threshP", "P-value:", min = 0, max = 1, value=0.1, step=0.05),
    uiOutput("uiExample") #, ########## choose action GoButton
    #uiOutput("testPair")
    #actionButton("goButton", "Go!")
  ),
  mainPanel(
    plotlyOutput("plot1"),
    #verbatimTextOutput("click"),
    #verbatimTextOutput("selectedValues"),
    plotlyOutput("boxPlot")
  )
))

dds <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[1]]
rld <- readRDS("/Users/lindz/BeeVirusDiet/beeDataDDSRLD.rds")[[2]]

myLevels <- unique(sapply(colnames(assay(rld)), function(x) unlist(strsplit(x,"[.]"))[1]))
myPairs <- list()

# Runs exact test on all pairs of groups and saves in list
k=1
for (i in 1:(length(myLevels)-1)){
  for (j in (i+1):(length(myLevels))){
    myPairs[[k]] <- paste0(myLevels[i], " and ", myLevels[j])
    k=k+1
  }
}

dat <- readRDS("/Users/lindz/BeeVirusDiet/beeVolcanoData.rds")

nCol = ncol(dat)
datFCP = dat[,(nCol-2*length(myLevels)+1):nCol]

# x-axis FC, y-axis pval
xMax = max(datFCP[,seq(1,ncol(datFCP),by=2)], na.rm=TRUE)
xMin = min(datFCP[,seq(1,ncol(datFCP),by=2)], na.rm=TRUE)
yMax = max(datFCP[,seq(2,ncol(datFCP),by=2)], na.rm=TRUE)
yMin = min(datFCP[,seq(2,ncol(datFCP),by=2)], na.rm=TRUE)
fcMax = ceiling(max(exp(xMax), 1/exp(xMin)))

server <- shinyServer(function(input, output) {
  #make dynamic slider
  output$slider <- renderUI({
    sliderInput("threshFC", "Fold change:", min=0, max=fcMax, value=ceiling((fcMax)/3), step=0.5)
  })
  
  output$selInput <- renderUI({
    selectInput("selPair", "Pairs:", myPairs, selected = unlist(myPairs[1]))
  })
  
  pairNum <- reactive(as.numeric(which(myPairs==input$selPair)))
  col1 <- reactive(colnames(dat)[nCol-2*length(myPairs)+2*pairNum()])
  col2 <- reactive(colnames(dat)[nCol-2*length(myPairs)+2*pairNum()-1])
  
  # datInput only validated once the go button is clicked
  datInput <- eventReactive(input$goButton, {
    dat[ which(dat[isolate(col1())] > -1* log10(input$threshP) & exp(abs(dat[isolate(col2())])) > input$threshFC), ]
  })
  
  output$uiExample <- renderUI({
    tags$span(
      tipify(actionButton("goButton", "Go!"), "Choose low p-value and high fold change", "This button is pointless!")
    )
  })
  
  output$plot1 <- renderPlotly({
    # will wait to render until datInput is validated
    plot_dat <- datInput()
    p <- qplot(plot_dat[[isolate(col2())]], plot_dat[[isolate(col1())]], xlim = c(xMin, xMax), ylim=c(yMin, yMax)) + xlab("log2(Fold change)") + ylab("-log10(p-value)")
    ggplotly(p)
  })
  
  d <- reactive(event_data("plotly_selected"))
  output$click <- renderPrint({
    if (is.null(d())){
      "Click on a state to view event data"
    }
    else{
      datInput()[d()$pointNumber+1,] #Working now
    }
  })
  
  pcpDat <- reactive(datInput()[d()$pointNumber+1,1:(ncol(dat)-2*length(myPairs))])
  #pcpDat <- eventReactive(input$goButton, {data.frame()})
  colNms <- colnames(dat[, 2:(ncol(dat)-2*length(myPairs))])
  nVar <- length(2:(ncol(dat)-2*length(myPairs)))
  
  boxDat <- dat[, 1:(ncol(dat)-2*length(myPairs))] %>% gather(key, val, -c(ID))
  
  #output$selectedValues <- renderPrint({str(boxDat)})
  colnames(boxDat)[2:3] <- c("Sample","Counts")
  BP <- ggplot(boxDat, aes(x = Sample, y = Counts)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90, hjust=1))
  ggBP <- ggplotly(BP)
  
  output$boxPlot <- renderPlotly({
    ggBP %>% onRender("
      function(el, x, data) {
      
      console.log(data.pcpDat)
      
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
      
      var tracePCPLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'blue',
      width: 1
      },
      opacity: 0.9,
      }
      Traces.push(tracePCPLine);
      }
      Plotly.addTraces(el.id, Traces);
      
      }", data = list(pcpDat = pcpDat(), nVar = nVar, colNms = colNms))})

  })

shinyApp(ui, server)
