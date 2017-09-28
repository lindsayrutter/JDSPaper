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