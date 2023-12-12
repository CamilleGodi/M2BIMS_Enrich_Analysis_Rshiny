# Rshiny application produced by Victor BAILLEUL (Université de Rouen Normandie)
# Last updated: 24/11/2023

# This application aimed to facilitate functional enrichment analysis from
# differential expression results

# Contact me at: victor.bailleul@univ-rouen.fr

library(shiny)
library(shinyjs)
library(shinydashboardPlus)
library(shinydashboard)
library(DT)      
library(plotly)   


################################################################################
################################################################################
################################################################################


dashboardPage( skin = "green",
  
  ### Header ###
  dashboardHeader(title = "VictoRshiny"),
  
  ### Sidebar ###
  dashboardSidebar( 
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("home")),
      
      fileInput(
        "inputFile", 
        "Select a CSV file:", 
        multiple = FALSE, 
        accept = c(".csv"),
        width = '100%'
      ),
      
      selectInput(
        "selectOrganism",
        "Select organism name:",
        c("Arabidopsis thaliana", "Caenorhabditis elegans","Danio rerio", "Drosophila melanogaster","Homo sapiens", "Mus musculus", "Saccharomyces cervisiae", "Xenopus laevis"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      
      sidebarMenu(
        menuItem("Whole Data Inspection", tabName = "WholeData", icon = icon("table")),
        menuItem("GO Term Enrichment", tabName = "GOTermEnrichment", icon = icon("sitemap"),
                 menuSubItem("ORA", tabName = "GOTermORA"),
                 menuSubItem("GSEA", tabName = "GOTermGSEA")),
        menuItem("Pathway Enrichment", tabName = "PathwayEnrichment", icon = icon("project-diagram"),
                 menuSubItem("ORA", tabName = "pathORA"),
                 menuSubItem("GSEA", tabName = "pathGSEA")),
        menuItem("About", tabName = "About", icon = icon("info-circle"))
      )
    )
  ),
  
  ##############################################################################
  ##############################################################################
  
  ### Body ###
  dashboardBody(
    tags$head(
      tags$style(HTML("
        #shiny-notification-panel {
          font-size: 50px;
          text-align: center;
          top: 50%;
          bottom: unset;
          left: 0;
          right: 0;
          margin-left: auto;
          margin-right: auto;
          width: 100%;
          max-width: 700px;
        }
      "))
    ),
    
    tabItems(
      
      ##########################################################################
      
      ### Home ###
      tabItem(
        tabName = "Home",
        h2("Home", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### Whole Data Inspection ###
      tabItem(
        tabName = "WholeData",
        h2("Whole data inspection", style = "text-align: center"),
        # Ligne 1
        fluidRow(
          # Plot
          box(
            title = "Volcano Plot",
            status = "warning", # Seulement a but esthétique
            solidHeader = TRUE, 
            collapsible = FALSE,
            plotlyOutput("volcanoPlot", height = "400px"),
            width = 8
          ),
          # Sliders
          box(
            title = "Options",
            status = "warning", # Seulement a but esthétique
            solidHeader = TRUE, 
            collapsible = FALSE,
            sliderInput("pValueCutoff", "P-Value cutoff", min = 0, max = 1, value = 0.05),
            sliderInput("logFC", "log2 FoldChange cutoff", min = 0, max = 10, value = 1),
            downloadButton("downloadFilteredTable", label = "Download filtered data table", icon = shiny::icon("download")),
            width = 4
          )
        ),
        # Ligne 2
        fluidRow(
          # Aperçu tableau
          box(
            title = "Filtered table preview",
            status = "warning", # Seulement a but esthétique
            DTOutput("dataPreview"),
            width = 12,
            solidHeader = TRUE 
            
          )
        )
      ),
      
      ##########################################################################
      
      ### GO Term Enrichment : ORA ###
      tabItem(
        tabName = "GOTermORA",
        h2("GO Term Enrichment : ORA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### GO Term Enrichment : GSEA ###
      tabItem(
        tabName = "GOTermGSEA",
        h2("GO Term Enrichment : GSEA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### Pathway Enrichment : ORA ###
      tabItem(
        tabName = "pathORA",
        h2("Pathway Enrichment : ORA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### Pathway Enrichment : GSEA ###
      tabItem(
        tabName = "pathGSEA",
        h2("Pathway Enrichment : GSEA", style = "text-align: center")
      )
      
      ##########################################################################
    )
  )
)

################################################################################
################################################################################
################################################################################
