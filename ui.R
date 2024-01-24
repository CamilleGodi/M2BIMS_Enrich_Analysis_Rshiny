# Rshiny application produced by Victor BAILLEUL (Université de Rouen Normandie)
# Last updated: 24/11/2023

# This application aimed to facilitate functional enrichment analysis from
# differential expression results

# Contact me at: victor.bailleul@univ-rouen.fr

source("global.R")

################################################################################
################################################################################
################################################################################

dashboardPage( skin = "purple",
  
  ### Header ###
  dashboardHeader(title = "VictoRshiny"),
  
  ### Sidebar ###
  dashboardSidebar( 
    sidebarMenu(
      menuItem("Home", tabName = "home_tab", icon = icon("home")),
      
      fileInput(
        "input_file", 
        "Select a CSV file:", 
        multiple = FALSE, 
        accept = c(".csv"),
        width = '100%'
      ),
      
      selectInput(
        "select_organism",
        "Select organism name:",
        c("Arabidopsis thaliana", "Caenorhabditis elegans","Danio rerio", "Drosophila melanogaster","Homo sapiens", "Mus musculus", "Saccharomyces cervisiae", "Xenopus laevis"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
      
      sidebarMenu(
        menuItem("Whole Data Inspection", tabName = "whole_data_tab", icon = icon("table")),
        menuItem("GO Term Enrichment", tabName = "GO_term_enrichment_tab", icon = icon("sitemap"),
                 menuSubItem("ORA", tabName = "GO_term_ORA_subtab"),
                 menuSubItem("GSEA", tabName = "GO_term_GSEA_subtab")),
        menuItem("Pathway Enrichment", tabName = "pathway_enrichment_tab", icon = icon("project-diagram"),
                 menuSubItem("ORA", tabName = "pathway_ORA_subtab"),
                 menuSubItem("GSEA", tabName = "pathway_GSEA_subtab")),
        menuItem("About", tabName = "about_tab", icon = icon("info-circle"))
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
        tabName = "home_tab",
        h2("Home", style = "text-align: center")
      ),
      
      
      ##########################################################################
      
      ### Whole Data Inspection ###
      tabItem(
        tabName = "whole_data_tab",
        h2("Whole data inspection", style = "text-align: center"),
        # Line 1
        fluidRow(
          # Plot
          box(
            title = "Volcano Plot",
            status = "warning", # Cosmetic purpose only
            solidHeader = TRUE, 
            collapsible = FALSE,
            plotlyOutput("volcano_plot", height = "400px"),
            width = 8
          ),
          # Sliders
          box(
            title = "Options",
            status = "warning", # Cosmetic purpose only
            solidHeader = TRUE, 
            collapsible = FALSE,
            sliderInput("padj_cutoff", "P-Value cutoff", min = 0, max = 1, value = 0.05),
            sliderInput("fc_cutoff", "log2 FoldChange cutoff", min = 0, max = 10, value = 1),
            downloadButton("download_filtered_table", label = "Download filtered data table", icon = shiny::icon("download")),
            width = 4
          )
        ),
        # Line 2
        fluidRow(
          # Table preview
          box(
            title = "Filtered table preview",
            status = "warning", # Cosmetic purpose only 
            DTOutput("data_preview_table"),
            width = 12,
            solidHeader = TRUE 
            
          )
        )
      ),
      
      ##########################################################################
      
      ### GO Term Enrichment : ORA ###
      tabItem(
        tabName = "GO_term_ORA_subtab",
        h2("GO Term Enrichment : ORA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### GO Term Enrichment : GSEA ###
      tabItem(
        tabName = "GO_term_GSEA_subtab",
        h2("GO Term Enrichment : GSEA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### Pathway Enrichment : ORA ###
      tabItem(
        tabName = "pathway_ORA_subtab",
        h2("Pathway Enrichment : ORA", style = "text-align: center")
      ),
      
      ##########################################################################
      
      ### Pathway Enrichment : GSEA ###
      tabItem(
        tabName = "pathway_GSEA_subtab",
        h2("Pathway Enrichment : GSEA", style = "text-align: center")
      )
      
      ##########################################################################
    )
  )
)

################################################################################
################################################################################
################################################################################
