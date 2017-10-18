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
    shinydashboard::menuItem("Application", tabName="hexPlot", selected=TRUE),
    shinydashboard::menuItem("About", tabName = "about")
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
               box(width = NULL, verbatimTextOutput("info1"), verbatimTextOutput("info2"), collapsible = TRUE, title = "Gene metrics", status = "primary", solidHeader = TRUE)))),
  
  shinydashboard::tabItem(tabName = "about",
    shiny::fluidRow("This application allows you to examine the relationship between all variables in your dataset with an interactive scatterplot matrix. Plotting points can obscure the number of genes in a given area due to overplotting. As a result, we use hexagon bins in the scatterplot matrix. If you hover over a given hexagon bin of interest, you can determine the number of genes in its area (Figure 1).", style='padding:10px;'),
    br(),
    shiny::fluidRow("You can also click on a given hexagon bin of interest to overlay the genes it contains across all scatterplots as orange points (Figure 2). Doing so will automatically overlay these same genes as orange parallel coordinate lines across a side-by-side boxplot of your data immediately below (Figure 3). Moreover, beneath that, you will see an output of the IDs of theses selection genes (Figure 4).", style='padding:10px;'),
    br(),
    shiny::fluidRow("The four figures below were created in simulated data drawn from the normal distribution for didactic purposes. We hope to improve upon this application by allowing users more customizing options, such as selecting hexagon bin size, changing color mappings, and providing a clear color legend.", style='padding:10px;'),
      br(),
      br(),
      div(p('Figure 1'), style="text-align: center;"),
      div(img(src='Figure1.png', style="width: 75%; height: 75%"), style="text-align: center;"),
      br(),
      br(),
      div(p('Figure 2'), style="text-align: center;"),
      div(img(src='Figure2.png', style="width: 75%; height: 75%"), style="text-align: center;"),
      br(),
      br(),
      div(p('Figure 3'), style="text-align: center;"),
      div(img(src='Figure3.png', style="width: 75%; height: 75%"), style="text-align: center;"),
      br(),
      br(),
      div(p('Figure 4'), style="text-align: center;"),
      div(img(src='Figure4.png', style="width: 75%; height: 75%"), style="text-align: center;")
  )))
  

dashboardPage(
  dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)