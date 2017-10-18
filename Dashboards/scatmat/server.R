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
library(RColorBrewer)
library(Hmisc)

load("soybean_cn.rda")
data <- soybean_cn
data <- data[,c(1:7)]
#data("soybean_cn_metrics")
#metrics <- soybean_cn_metrics
datCol <- colnames(data)[-which(colnames(data) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
rm(soybean_cn)
#metrics[[1]] <- metrics[[1]][which(metrics[[1]]$PValue<0.01),]
#metrics[[1]] <- metrics[[1]][which(metrics[[1]]$ID %in% data$ID),]
#data = data[which(data$ID %in% metrics[[1]]$ID),]
#myMetrics <- colnames(metrics[[1]])[-which(colnames(metrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

server <- function(input, output, session) {

output$scatMatPlot <- renderPlotly({
withProgress(message = 'Making plot', detail = 'This could take about 30 seconds even after this message disappears', value = 0, {

  ################################ Prepare scatterplot matrix
  ###########################################################
  
  maxVal = max(abs(data[,-1]))
  maxRange = c(-1*maxVal, maxVal)
  
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    h <- hexbin(x=x, y=y, xbins=15, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-0.5, maxRange[2]+0.5), ylim = c(-0.5, maxRange[2]+0.5)) + theme(legend.position="none")
    p
  }
  
  p <- ggpairs(data[,-1], lower = list(continuous = my_fn))
  pS <- p
  
incProgress(1/4)
  
  ggPS <- ggplotly(pS, width=700, height=600)

incProgress(1/4)  
    
  myLength <- length(ggPS[["x"]][["data"]])
  for (i in 1:myLength){
    item =ggPS[["x"]][["data"]][[i]]$text[1]
    if (!is.null(item)){
      if (!startsWith(item, "co")){
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
      }}
    hexHover = ggPS[["x"]][["data"]][[i]]$text
    if (!is.null(hexHover) && grepl("hexID", hexHover)){
      ggPS[["x"]][["data"]][[i]]$text <- strsplit(hexHover, "<")[[1]][1]
      ggPS[["x"]][["data"]][[i]]$t2 <- hexHover
      ggPS[["x"]][["data"]][[i]]$hoverinfo <- "text"
    }
  }
  
  for(i in 2:(p$nrow)) {
    for(j in 1:(p$nrow-1)) {
      data[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
    }
  }
  
incProgress(1/4)
ggPS2 <- ggPS %>% onRender("
    function(el, x, data) {
    
    function range(start, stop, step){
    var a=[start], b=start;
    while(b<stop){b+=step;a.push(b)}
    return a;
    };
    
    len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
    AxisNames = [];
    for (i = 1; i < (len+1); i++) {
    AxisNames.push(Object.keys(data[0])[i]);
    }
    noPoint = x.data.length;
    
    el.on('plotly_click', function(e) {
    
    if (x.data.length > noPoint){
    Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
    }
    
    xVar = (e.points[0].xaxis._id).replace(/[^0-9]/g,'')
    if (xVar.length == 0) xVar = 1
    yVar = (e.points[0].yaxis._id).replace(/[^0-9]/g,'')
    if (yVar.length == 0) yVar = 1
    myX = len + 1 - (yVar - len * (xVar - 1))
    myY = xVar
    cN = e.points[0].curveNumber
    split1 = (x.data[cN].text).split(' ')
    hexID = (x.data[cN].t2).split(':')[2]
    counts = split1[1].split('<')[0]
    var selRows = [];
    data.forEach(function(row){
    if(row[myX+'-'+myY]==hexID) selRows.push(row);
    });
    selID = []
    for (a=0; a<selRows.length; a++){
    selID.push(selRows[a]['ID'])
    }
    // Save selected row IDs for PCP
    Shiny.onInputChange('selID', selID);
    
    var Traces = [];
    var i=0;
    var k=1;
    while ((i*len+k)<=Math.pow((len-1),2)) {
    var xArr = [];
    for (a=0; a<selRows.length; a++){
    xArr.push(selRows[a][AxisNames[i]])
    }
    while ((i+k)<len){
    var yArr = [];
    for (a=0; a<selRows.length; a++){
    yArr.push(selRows[a][AxisNames[(len-k)]])
    }
    //console.log(['xArr', xArr])
    //console.log(['yArr', yArr])
    var trace = {
    x: xArr,
    y: yArr,
    mode: 'markers',
    marker: {
    color: 'orange',
    size: 6
    },
    xaxis: 'x' + (i+1),
    yaxis: 'y' + (i*len+k),
    hoverinfo: 'none'
    };
    Traces.push(trace);
    k++;
    }
    i++;
    k=1;
    }
    Plotly.addTraces(el.id, Traces);
    })}
    ", data = data)
incProgress(1/4)
ggPS2
})
    
})

  selID <- reactive(input$selID)
  
  pcpDat <- reactive(data[which(data$ID %in% selID()), c(1:7)])
  output$selectedValues <- renderPrint({pcpDat()$ID})
  colNms <- colnames(data[, c(2:7)])
  print(c('colNms', colNms))
  
  boxDat <- data[, c(1:7)] %>% gather(Sample, Count, -c(ID))
  print(c('boxDat', str(boxDat)))
  #colnames(boxDat) <- c("ID", "Sample", "Count")
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
      
      }", data = list(pcpDat = pcpDat(), nVar = 7, colNms = colNms))}) # should be p$nrow not 7

}