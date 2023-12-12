# A shiny application for RNA-seq data analysis : ui code
# Author: Komlan Dieu-Donné TOTO
#email: komlan-dieu-donne.toto@univ-rouen.fr
# Date: 24/11/2023

# Infos: To run this application, you need to install the following packages:
# shiny, shinydashboard, shinythemes, ggplot2, plotly, DT, data.table

library(shinydashboard)
library(shinythemes)
library(ggplot2)
library(plotly)
library(DT)
library(data.table)


#   This function is used to get the organism choices sorted alphabetically.
getOrganismChoices <- function() {
  organism_choices <- c("Homo sapiens", "Escherichia coli", "Mus musculus", "Rattus norvegicus", "Drosophila melanogaster", "Caenorhabditis elegans", "Danio rerio")   
  return(sort(organism_choices))
}


# The UI definition for the Shiny dashboard page.
# This code defines the layout and components of the user interface.
ui <- dashboardPage(
  dashboardHeader(title = 'My Dashboard'),
  dashboardSidebar(
    sidebarMenu( # definition of the menu of the interface
      menuItem(
        "Home",
        tabName = "home",
        icon = icon("home")
      ),
      fileInput("file", label = h6("Select a CSV file :")),
      selectInput(
        "select",
        label = h6("Select an organism name :"),
        choices = getOrganismChoices(),  # call for the getOrganismChoices() function
        selected = 1 # only one choice can be made 
      ),
      menuItem(
        "Whole Data Inspection",
        tabName = "Whole Data Inspection",
        icon = icon("database") 
      ),
      menuItem(
        "Go Term Enrichment",
        tabName = "Go Term Enrichment",
        icon = icon("sitemap")
      ),
      menuItem(
        "Pathway Enrichment",
        tabName = "Pathway Enrichment",
        icon = icon("chart-pie")
      ),
      fluidRow(
        box( # definition of the box
          title = NULL,
          height = 50, 
          width = 11,
          background = "black"
        )
      ),
      fluidPage(
        menuItem(
          "About",
          tabName = "About",
          icon = icon("th-large")
        )
      )
    )
  ),
  dashboardBody( # definition of the body of the interface
    h2("Whole Data Inspection", class = "text-center"),
    fluidRow(
      box(title = "Volcano Plot", plotlyOutput("volcanoPlot"), height = 455, status = "primary", solidHeader = TRUE), # box for the volcano plot
      box( # box for the options
        title = "Options", height = 455, status = "primary", solidHeader = TRUE,
        sliderInput("Pvalue_cutoff_from_input", "P-value cutoff from input:", min = 0, max = 1, value = 0.31, step = 0.01),
        br(), br(),  # 2 line breaks
        sliderInput("log2_FoldChange_cutoff_from_input", "log2 FoldChange cutoff from input:", min = 0, max = 5, value = 0.5, step = 0.1),
        br(),br(),br(),br(),br(),br(), # 6 line breaks
        downloadButton("downloadVolcanoPlot", "Download Volcano Plot")
      )
    ),
    fluidRow( # fluid row for the table
      box(
        title = "Table",
        DT::dataTableOutput("table"), # output of the table
        height = 500,
        width = 12,
        status = "primary",
        solidHeader = TRUE
      )
    )
  )
)

