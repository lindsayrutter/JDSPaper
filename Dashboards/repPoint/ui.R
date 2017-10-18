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

load("soybean_cn.rda")
dat <- soybean_cn
dat <- dat[,1:7]
datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
load("soybean_cn_metrics.rda")
metrics <- soybean_cn_metrics
dat = dat[which(dat$ID %in% metrics[[1]]$ID),]
myMetrics <- colnames(metrics[[1]])[-which(colnames(metrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

sidebar <- dashboardSidebar(
  width = 180,
  hr(),
  sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="hexPlot", selected=TRUE) #, #hexPlot
    #shinydashboard::menuItem("About", tabName = "about")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "hexPlot",
      fluidRow(
        column(width = 4, 
         box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
           selectizeInput("selPair", "Treatment pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
           selectInput("selMetric", "Metric:", choices = myMetrics),
           selectInput("selOrder", "Metric order:", choices = c("Increasing", "Decreasing")),
           numericInput("binSize", "Hexagon size:", value = 10, min = 1),
           numericInput("pointSize", "Point size:", value = 6, min = 1),
           actionButton("goButton", "Plot gene!"))),
        column(width = 8,
           box(width = NULL, plotlyOutput("hexPlot"), collapsible = FALSE, background = "black", title = "Replicate point plot", status = "primary", solidHeader = TRUE))),
    
      fluidRow(
        column(width = 12,
         box(width = NULL, plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))),
      
      fluidRow(
        column(width = 12,
               box(width = NULL, verbatimTextOutput("info1"), collapsible = TRUE, title = "Gene metrics", status = "primary", solidHeader = TRUE))))))

dashboardPage(
  dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)