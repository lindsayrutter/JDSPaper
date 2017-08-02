library(shinydashboard)
library(shiny)
library(plotly)
library(hexbin)
library(htmlwidgets)
library(dplyr)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(Hmisc)
library(GGally)

load("bindataL120.Rda")
dat <- bindata
dat[,2:ncol(dat)] = log2(dat[,2:ncol(dat)]+1)
rm(bindata)
datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
values <- reactiveValues(x=0, selPair=NULL, goButton=NULL)

shinyServer(function(input, output, session){
  observeEvent(input$selPair1, values$selPair <- input$selPair1, ignoreNULL = FALSE)
  observeEvent(input$selPair2, values$selPair <- input$selPair2, ignoreNULL = FALSE)
  
  observeEvent(input$goButton1, values$goButton <- input$goButton1)
  observeEvent(input$goButton2, values$goButton <- input$goButton2)
  
  observe({x <- values$selPair
  updateSelectizeInput(session, "selPair1", "Treatment(s):", choices = myPairs, options = list(maxItems = 2), selected = x)
  updateSelectizeInput(session, "selPair2", "Treatment(s):", choices = myPairs, options = list(maxItems = 2), selected = x)})
  
  # Create data subset based on two letters user chooses
  datSel2 <- eventReactive(values$selPair, {
    validate(need(length(values$selPair) == 1 || length(values$selPair) == 2, "Select at least one treatment."))
    sampleIndex <- reactive(which(sapply(colnames(dat), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1], values$selPair[2])))
    dat[,c(1, sampleIndex())]
  }, ignoreNULL = FALSE) #, ignoreNULL = FALSE)
  
  datSel <- eventReactive(values$goButton, {
    #req(input$selPair1)
    #req(input$selPair2)
    datSel2()
  })
  
  binSize <- eventReactive(values$goButton, {
    input$binSize
  })
  
  #observeEvent(values$goButton, { 
  
  # output$text1 <- renderText({
  #   if (input$num < 5) {
  #     hide("text1")
  #   } else {
  #     show("text1")
  #   }
  #   paste0("The number is ", input$num)
  # })
  
  output$hexPlot <- renderPlotly({
    
    
    withProgress(message = 'Creating plot', value = 0, {
      # Number of times we'll go through the "creating plot" loop
      n <- 2
      
      maxVal = max(datSel()[,-1])
      minVal = min(datSel()[,-1])
      maxRange = c(minVal, maxVal)
      xbins=binSize()
      buffer = maxRange[2]/(xbins)
      
      my_fn <- function(data, mapping, ...){
        xChar = as.character(mapping$x)
        yChar = as.character(mapping$y)
        x = data[,c(xChar)]
        y = data[,c(yChar)]
        h <- hexbin(x=x, y=y, xbins=binSize(), shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
        hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
        attr(hexdf, "cID") <- h@cID
        
        hexdf$countColor <- cut2(hexdf$counts, g=6)
        hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor), function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE), ",")[[1]][1], 2))))
        hexdf$countColor2 <- factor(hexdf$countColor2, levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))
        for (i in 1:(length(levels(hexdf$countColor2))-1)){
          levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],"-",levels(hexdf$countColor2)[i+1])
        }
        levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <- paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")
        my_breaks = levels(hexdf$countColor2)
        clrs <- brewer.pal(length(my_breaks)+3, "Blues")
        clrs <- clrs[3:length(clrs)]
        
        p <- ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Cases count") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_fixed(xlim = c(-1*buffer, maxRange[2]+buffer), ylim = c(-1*buffer, maxRange[2]+buffer))
        p
      }
      
      p <- reactive(ggpairs(datSel()[,2:ncol(datSel())], lower = list(continuous = my_fn)))
      
      # Update progress bar
      incProgress(1/n)    
      
      gP <- eventReactive(p(), {
        gP <- ggplotly(p(), height = 600, width = 600) #  height = 400
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
      
      # Update progress bar
      incProgress(1/n)    
      plotlyHex()
      
    })
    
  })
  
  output$boxPlot <- renderPlotly({
    
    withProgress(message = 'Creating plot', value = 0, {
      # Number of times we'll go through the "creating plot" loop
      n <- 2
      
      nVar = reactive(ncol(datSel()))
      colNms <- reactive(colnames(datSel()[, c(2:nVar())]))
      
      boxDat <- eventReactive(datSel(), {
        boxDat <- datSel()[, c(1:nVar())] %>% gather(key, val, -c(ID))
        colnames(boxDat) <- c("ID", "Sample", "Count")
        boxDat
      })
      
      BP <- reactive(ggplot(boxDat(), aes(x = Sample, y = Count)) + geom_boxplot() + labs(y = "Read count"))
      # Update progress bar
      incProgress(1/n) 
      ggBP <- reactive(ggplotly(BP(), width=600, height = 400) %>% config(displayModeBar = F, staticPlot = T))
      # Update progress bar
      incProgress(1/n) 
      ggBP()
    })
  })
})