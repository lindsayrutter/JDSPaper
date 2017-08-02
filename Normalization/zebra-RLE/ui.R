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

load("zebraData.Rda")
dat <- zebraData
dat[,2:ncol(dat)] = log2(dat[,2:ncol(dat)]+1)
rm(zebraData)

# Add pseudocount, and do RLE (mean of means is not 0. It is -0.1452344; mean(rowMeans(dat[,2:ncol(dat)]))). Scatterplot looks odd. PCP is centered at 0, but lots of outliers.
#dat[,2:ncol(dat)] <- dat[,2:ncol(dat)] + 0.1
#dat[,2:ncol(dat)] <- log(dat[,2:ncol(dat)]/rowMedians(as.matrix(dat[,2:ncol(dat)])))

# Same problem
#dat[,2:ncol(dat)] <- dat[,2:ncol(dat)] + 0.000001
#dat[,2:ncol(dat)] <- log((dat[,2:ncol(dat)])/rowMedians(as.matrix(dat[,2:ncol(dat)])))

# Remove any rows that have any columns with value of 0
dat <- dat[which(apply(dat[,2:ncol(dat)]!=0,1,all)),]
dat[,2:ncol(dat)] <- log(dat[,2:ncol(dat)]/rowMedians(as.matrix(dat[,2:ncol(dat)])))

#dat[,2:ncol(dat)] <- log((dat[,2:ncol(dat)]+0.1)/rowMedians(as.matrix(dat[,2:ncol(dat)])))

datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
values <- reactiveValues(x=0, selPair=NULL)

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
        column(width = 3, 
           box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
             selectizeInput("selPair1", "Treatment(s):", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
             numericInput("binSize", "Hexagon size:", value = 10, min = 1),
             actionButton("goButton1", "Plot case!")),
             textOutput("text1")),
        column(width = 9,
           box(width = NULL, height = 675, plotlyOutput("hexPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE)))),
    
    tabItem(tabName = "boxPlot",
      fluidRow(
        column(width = 3, 
           box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
             selectizeInput("selPair2", "Treatment(s):", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
             actionButton("goButton2", "Plot case!"))),
        column(width = 9,
           box(width = NULL, plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))))
  )
)

dashboardPage(
  dashboardHeader(title = "Plotting", titleWidth = 180),
  sidebar,
  body
)