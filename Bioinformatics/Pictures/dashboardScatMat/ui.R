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

set.seed(3)
data = data.frame(ID = paste0("ID", 1:1999), A = abs(rnorm(1999)), B = abs(rnorm(1999)), C=abs(rnorm(1999)))
data$ID = as.character(data$ID)
xbins = 15

sidebar <- shinydashboard::dashboardSidebar(
  width = 180,
  shiny::hr(),
  shinydashboard::sidebarMenu(id="tabs",
  shinydashboard::menuItem("Application", tabName="scatMatPlot", selected=TRUE), #hexPlot
  shinydashboard::menuItem("About", tabName = "about") #boxPlot
)
)

body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    
shinydashboard::tabItem(tabName = "scatMatPlot",
shiny::fluidRow(
shiny::column(width = 12, shinydashboard::box(width = 660, height = 660, plotly::plotlyOutput("scatMatPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),

shiny::fluidRow(
shiny::column(width = 12,
shinydashboard::box(width = NULL, plotly::plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Boxplot", status = "primary", solidHeader = TRUE))),

shiny::fluidRow(
shiny::column(width = 12,
shinydashboard::box(width = NULL, shiny::verbatimTextOutput("selectedValues"), collapsible = TRUE, title = "Selected Gene IDs", status = "primary", solidHeader = TRUE)))),

shinydashboard::tabItem(tabName = "about",
shiny::fluidRow("This application allows you to examine the relationship between all variables in your dataset with an interactive scatterplot matrix. Plotting points can obscure the number of genes in a given area due to overplotting. As a result, we use hexagon bins in the scatterplot matrix. If you hover over a given hexagon bin of interest, you can determine the number of genes in its area (Figure 1).", style='padding:10px;'),
br(),
shiny::fluidRow("You can also click on a given hexagon bin of interest to overlay the genes it contains across all scatterplots as orange points (Figure 2). Doing so will automatically overlay these same genes as orange parallel coordinate lines across a side-by-side boxplot of your data immediately below (Figure 3). Moreover, beneath that, you will see an output of the IDs of theses selection genes (Figure 4).", style='padding:10px;'),
br(),
shiny::fluidRow("The four figures below were created in simulated data drawn from the normal distribution for didactic purposes. We hope to improve upon this application by allowing users more customizing options, such as selecting hexagon bin size, changing color mappings, and providing a clear color legend.", style='padding:10px;'),
br(),
br(),
div(p('Figure 1'), style="text-align: center;"),
div(img(src='Figure1.png', style="width: 60%; height: 60%"), style="text-align: center;"),
br(),
br(),
div(p('Figure 2'), style="text-align: center;"),
div(img(src='Figure2.png', style="width: 60%; height: 60%"), style="text-align: center;"),
br(),
br(),
div(p('Figure 3'), style="text-align: center;"),
div(img(src='Figure3.png', style="width: 60%; height: 60%"), style="text-align: center;"),
br(),
br(),
div(p('Figure 4'), style="text-align: center;"),
div(img(src='Figure4.png', style="width: 60%; height: 60%"), style="text-align: center;")
))
)

ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Overlaying cases", titleWidth = 180),
  sidebar,
  body
)