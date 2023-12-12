################################################################################
####################### M2.1 BIMS - PROJET R SHINY #############################
##################### Camille GODI - Promo 2022-2025 ###########################
################################################################################

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyalert)

################################################################################

dashboardPage( skin = "purple",

  ### HEADER ###
  dashboardHeader(title = "DE Viz'Tools"),
  
  ### SIDEBAR ###
  dashboardSidebar(
    
    useShinyjs(),
    
    sidebarMenu(
      # Main menu / home
      menuItem("Home", tabName = "tab_data_inspection", icon = shiny::icon("home"))),
    
      # Organism selection TBA: list from ENSEMBL (?)
      selectizeInput("organism", "1/ Select the organism :", 
                     options = list(placeholder = "", onInitialize = I('function() {this.setValue(""); }') ),
                     choices=c( "Caenorhabditis_elegans", "Homo_sapiens", "Mus_musculus", "Pan_troglodytes", "Rattus_norvegicus", "Saccharomyces_cerevisiae" ) ),
    
      # Input File (csv) selection 
      disabled(fileInput(inputId = "csv_upload", 
                label = "2/ Select a csv file :")),

      # Sub-menus
      sidebarMenu(
        menuItem("Whole data Inspection", tabName = "tab_data_inspection", icon = shiny::icon("database") ),
        menuItem("GO Term Enrichment", tabName = "tab_GO_term", icon = shiny::icon("network-wired") ),
        menuItem("Pathway Enrichment", tabName = "tab_pathway", icon = shiny::icon("pie-chart") ),
        menuItem("About", tabName = "tab_about", icon = shiny::icon("question") )
      )
      
    ),
  
  ### BODY ###
  dashboardBody(
    tabItems(
      ### TAB : DATA INSPECTION / HOME
      tabItem( tabName = "tab_data_inspection", 
        # Dashboard body title
        h2( "Whole Data Inspection", align="center" ),
        
        fluidRow(
          # Volcano Plot box
          box(title = "Volcano Plot", width = 8,
              solidHeader = TRUE, status = "primary",
              plotOutput("volcano_plot", height = 500, 
                         dblclick = "volcano_plot_dblclick",
                         brush = brushOpts( id = "volcano_plot_brush", resetOnNew = TRUE ))
              ),
          
          # Options and download button
          box(title = "Options and info", width = 4,
              solidHeader = TRUE, status = "primary",
              
              # Option sliders for p-value and log2 cutoffs
              h4("Plot thresholds :"),
              sliderInput("p_value_slider", "p-value cutoff from input", 
                          min = 1e-10, value = 1e-3, max = 5e-2, step = 1e-3),
              sliderInput("log2_slider", "log2 FoldChange cutoff from input", 
                          min = 0, max = 5, value = 1, step = 0.1),
              
              # Download button
              h4("Downloads :"),
              disabled( downloadButton("download_volcano", label = "Download volcano plot", icon = shiny::icon("download")) ),
              disabled( downloadButton("download_table", label = "Download data table", icon = shiny::icon("download")) ),
              
              # Control infos
              h4("Controls :"),
              h5("Drag mouse on the plot to select a zone, double-click to zoom on it, double-click again to reset.")
          )
          
        ),
        
        # Display csv table once one is loaded
        fluidPage(
          box(
            title = "Table", solidHeader = TRUE, status = "primary",
            width = 12, collapsible = TRUE,
            dataTableOutput("input_csv_table")
            )
          )
        ),
      
      ### TAB : GO TERM ENRICHMENT
      tabItem(
        tabName = "tab_GO_term",
        h2( "GO Term Enrichment", align="center" )    
      ),
      
      ### TAB : PATHWAY ENRICHMENT
      tabItem(
        tabName = "tab_pathway",
        h2( "Pathway Enrichment", align="center" )    
      ),
      
      ### TAB : ABOUT
      tabItem(
        tabName = "tab_about",
        h2( "About", align="center" )    
      )
    )
  )
  
)
