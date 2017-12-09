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

sidebar <- shinydashboard::dashboardSidebar(
  width = 180,
  shiny::hr(),
  shinydashboard::sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="scatMatPlot"), #hexPlot
    shinydashboard::menuItem("About", tabName = "about", selected=TRUE) #boxPlot
  )
)

body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    
    shinydashboard::tabItem(tabName = "scatMatPlot",
      shiny::fluidRow(
        shiny::column(width = 12, shinydashboard::box(width = 660, height = 660, plotly::plotlyOutput("scatMatPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),
      
      shiny::fluidRow(
        shiny::column(width = 12, shinydashboard::box(width = NULL, plotly::plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Boxplot", status = "primary", solidHeader = TRUE))),
      
      shiny::fluidRow(
        shiny::column(width = 12, shinydashboard::box(width = NULL, shiny::verbatimTextOutput("selectedValues"), collapsible = TRUE, title = "Selected Gene IDs", status = "primary", solidHeader = TRUE)))),

shinydashboard::tabItem(tabName = "about",
      shiny::fluidRow("This application allows you to examine the relationship between all variables in your dataset with an interactive scatterplot matrix. Plotting points can obscure the number of genes in a given area due to overplotting. As a result, we use hexagon bins in the scatterplot matrix. If you hover over a given hexagon bin of interest, you can determine the number of genes in its area (Figure 1).", style='padding:10px;'),
      br(),
      shiny::fluidRow("You can also click on a given hexagon bin of interest to overlay the genes it contains across all scatterplots as orange points (Figure 2). Doing so will automatically overlay these same genes as orange parallel coordinate lines across a side-by-side boxplot of your data immediately below (Figure 3). Moreover, beneath that, you will see an output of the IDs of theses selection genes (Figure 4).", style='padding:10px;'),
      br(),
      shiny::fluidRow("The four figures below were created in simulated data drawn from the normal distribution for didactic purposes. We hope to improve upon this application by allowing users more customizing options, such as selecting hexagon bin size, changing color mappings, and providing a clear color legend.", style='padding:10px;'),
      br(),
      shiny::fluidRow("Please use the 'Application' tab to interact with published RNA-seq data of soybean developmental stages (Brown and Kudson, 2015).", strong("When you hover over a hexagon bin, if it appears to contain thousands of genes, please avoid clicking on it."), "Plotting such a large number of genes could crash the program. In any case, hexagon bins that contain that many genes are unlikely to be of interest anyway because they usually occur close to the", em("x=y"), "line. Instead, we are typically interested in hexagon bins that deviate from the" ,em("x=y"), "line, and these bins tend to only contain a handful of genes in the first place.", style='padding:10px;'),
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
    ))
)

ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)