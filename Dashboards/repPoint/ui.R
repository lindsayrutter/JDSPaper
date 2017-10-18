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
    shiny::fluidRow("This application allows users to superimpose a differentially expressed gene of interest onto a replicate point plot. In the replicate point plot, each gene can be plotted once for each combination of sample pairs between treatment groups. For example, the data we use below contains two treatments (S1 and S2), each with three replicates. Hence, there are six ways to pair a replicate from one treatment group with a replicate from the other treatement group (S1.1 and S2.1, S1.1 and S2.2, S1.1 and S2.3, S1.2 and S2.1, S1.2 and S2.2, and S1.2 and S2.3).", style='padding:10px;'),
    br(),
    shiny::fluidRow("Each gene for this dataset could be plotted as six points to construct the replicate point plot. However, with 73,320 genes in this dataset, we would have 439,920 points. In interactive versions of the plot, this would reduce the speed of the funcionality as well as cause overplotting problems. As a result, we use hexagon bins to construct the background of the replicate point plot as is shown in the right side of Figure 1.", style='padding:10px;'),
    br(),
    shiny::fluidRow("This application comes with several input fields as shown on the left side of Figure 1. The user must choose exactly two treatment groups in the 'Treatment Pairs' tab. They must choose an order (increasing or decreasing) in which to scroll through genes by a metric of choice. We see in Figure 1 that the user chose to superimpose the genes by increasing order of FDR values between S1 and S2.", style='padding:10px;'),
    br(),
    shiny::fluidRow("Upon making these decisions, the user can then select the 'Plot gene!' button to superimpose each gene one by one onto the replicate point plot. In Figure 1, we see this as six orange points that show high values for S2 and low values for S1. This gene is also superimposed as an orange parallel coordinate line on top of a box plot as shown in Figure 2. Moreover, the gene ID and its metric values are output as shown in Figure 3. We can determine that the gene we are viewing ranks third by our metric and order parameters. This  means the user has hit the 'Plot gene!' button three times now for this set of parameters and that this gene has the third lowest FDR value between S1 and S2 for this dataset.", style='padding:10px;'),
    
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
    div(img(src='Figure3.png', style="width: 75%; height: 75%"), style="text-align: center;")
  )))
  

dashboardPage(
  dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)