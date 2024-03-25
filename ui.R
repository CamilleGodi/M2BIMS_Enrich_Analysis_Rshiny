# Rshiny application produced by Victor BAILLEUL (Université de Rouen Normandie)
# Last updated: 24/11/2023

# This application aimed to facilitate functional enrichment analysis from
# differential expression results

# Contact me at: victor.bailleul@univ-rouen.fr

source("global.R")
################################################################################
################################################################################
################################################################################
dashboardPage(
  skin = "purple",
  
  ### Header ###
  dashboardHeader(title = "Rshiny"),
  
  ### Sidebar ###
  dashboardSidebar( 
    sidebarMenu(
      menuItem("Home", tabName = "home_tab", icon = icon("home")),
      
      fileInput(
        "input_file", 
        "Select a CSV file:", 
        multiple = FALSE, 
        accept = c(".csv"),
        width = '100%'),
      
      selectInput(
        "select_organism",
        "Select organism name:",
        c("Arabidopsis thaliana", "Caenorhabditis elegans","Danio rerio", "Drosophila melanogaster","Homo sapiens", "Mus musculus", "Saccharomyces cervisiae", "Xenopus laevis"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL),
      
      sidebarMenu(
        menuItem("Whole Data Inspection", tabName = "whole_data_tab", icon = icon("table")),
        menuItem("GO Term Enrichment",    tabName = "GO_term_enrichment_tab", icon = icon("sitemap"),
                 menuSubItem("ORA",  tabName = "GO_term_ORA_subtab"),
                 menuSubItem("GSEA", tabName = "GO_term_GSEA_subtab")),
        menuItem("Pathway Enrichment", tabName = "pathway_enrichment_tab", icon = icon("project-diagram"),
                 menuSubItem("ORA",  tabName = "pathway_ORA_subtab"),
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
    "))),
    
    tabItems(
      
      ### Home ###
      tabItem(
        tabName = "home_tab",
        h2("Home", style = "text-align: center"),
        box(
          title = "Tutorial",
          status = "warning", # Cosmetic purpose only : orange box
          solidHeader = TRUE, 
          width = 12,
          h4("1/ Select a CSV (or CSV2) file. It must at least have the following columns : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'."),
          h4("2/ Select the organism from which the data originates."),
          h4("3/ Explore your data through the 'Whole data inspection' tab, and/or perform desired analysis.")
        )
      ),
      
      tabItem(
        tabName = "about_tab",
        h2("About the project", style = "text-align: center"),
        box(
          width = 12,
          strong("Goal :"),
          p("The goal of this application is to facilitates functional enrichment analysis from differential expression results, and allow for quick and interactive visualization."),
          br(),
          strong("Authors :"),
          p("Victor BAILLEUL ( victor.bailleul@univ-rouen.fr )"),
          p("Camille GODI ( camille.godi@univ-rouen.fr )"),
          p("Benjamin MARSAC ( benjamin.marsac@univ-rouen.fr )"),
          p("Komlan Dieu-Donné TOTO ( komlan-dieu-donne.toto@univ-rouen.fr )"),
          br(),
          p("This app is the result of a group work in second year of Bioinformatics Master's Degree, 'BIMS', year 2023-2024s, Université de Rouen Normandie ( URN ).")
        )
      ),
      
      ##########################################################################
      
      ### Whole Data Inspection ###
      tabItem(
        tabName = "whole_data_tab",
        h2("Whole data inspection", style = "text-align: center"),
        fluidRow(
          box(
            title = "Volcano Plot",
            status = "warning", # Cosmetic purpose only
            solidHeader = TRUE, 
            collapsible = FALSE,
            plotlyOutput("volcano_plot", height = "400px"),
            width = 6
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
            width = 6)),
        fluidRow(
          box(
            title = "Filtered table preview",
            status = "warning", # Cosmetic purpose only 
            DTOutput("data_preview_table"),
            width = 12,
            solidHeader = TRUE))),
      
      ### GO Term Enrichment : ORA ###
      tabItem(
        tabName = "GO_term_ORA_subtab",
        h2("GO Term Enrichment : ORA", style = "text-align: center"),
        
        fluidRow(
          # GENE ONTOLOGY SETTINGS
          column(
            width = 8,
            box(
              title = "GENE ONTOLOGY SETTINGS",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  selectInput("goAnnotationORA", "Select a GO Annotation:",
                              choices = c("Biological Process", "Molecular Function", "Cellular Component"),
                              selected = "Biological Process"))))),
          
          # GO Level Selection
          column(
            width = 8,
            box(
              title = "GO Level Selection",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("levelSelectionORA", "Select a GO Level:", 
                               choices = c("One-level GO ORA" = "OneLevelORA", "All-level GO ORA" = "AllLevelORA")))),
              # Conditional Panel for OneLevelORA
              conditionalPanel(
                condition = "input.levelSelectionORA == 'OneLevelORA'",
                fluidRow(
                  column(
                    width = 12,
                    box(
                      status = "primary",
                      solidHeader = TRUE,
                      width = 12,
                      fluidRow(
                        column(
                          width = 12,
                          sliderInput("levelSliderORA", "Select a GO level:", min = 1, max = 7, value = 1))))))))),
          
          # DEG Selection
          column(
            width = 8,
            box(
              title = "DEG profile Selection",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DEGSelection", "Select the correct DEG processing:", 
                               choices = c("Over expressed DEG only" = "OverDEG", "Under expressed DEG only" = "UnderDEG", "Both" = "BothDEG"))))))),
        
        # Parameters
        fluidRow(
          column(
            width = 8,
            box(
              title = "Parameters",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                    sliderInput("PValueORA", "Select a P-Value:", min = 0, max = 1, value = 0.05)),
                column(
                  width = 12,
                    sliderInput("QValueORA", "Select a Q-Value:", min = 0, max = 1, value = 0.05))))))),
      
      ### GO Term Enrichment : GSEA ###
      tabItem(
        tabName = "GO_term_GSEA_subtab",
        h2("GO Term Enrichment : GSEA", style = "text-align: center"),
        
        fluidRow(
          # GENE ONTOLOGY SETTINGS
          column(
            width = 8,
            box(
              title = "GENE ONTOLOGY SETTINGS",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  selectInput("goAnnotationGSEA", "Select a GO Annotation:",
                              choices = c("Biological Process", "Molecular Function", "Cellular Component"),
                              selected = "Biological Process"))))),
          
          # GO Level Selection
          column(
            width = 8,
            box(
              title = "GO Level Selection",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("levelSelectionGSEA", "Select a GO Level:", 
                               choices = c("One-level GO GSEA" = "OneLevelGSEA", "All-level GO GSEA" = "AllLevelGSEA")))),
              # Conditional Panel for OneLevelGSEA
              conditionalPanel(
                condition = "input.levelSelectionGSEA == 'OneLevelGSEA'",
                fluidRow(
                  column(
                    width = 12,
                    box(
                      status = "primary",
                      solidHeader = TRUE,
                      width = 12,
                      fluidRow(
                        column(
                          width = 12,
                          sliderInput("levelSliderGSEA", "Select a GO level:", min = 1, max = 7, value = 1))))))))),
          
          # DEG Selection
          column(
            width = 8,
            box(
              title = "DEG profile Selection",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DEGSelection", "Select the correct DEG processing:", 
                               choices = c("Over expressed DEG only" = "OverDEG", "Under expressed DEG only" = "UnderDEG", "Both" = "BothDEG")))))),
          
        # Parameters
          column(
            width = 8,
            box(
              title = "Parameters",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  sliderInput("PValueCutoffGSEA", "Select a P-Value:", min = 0, max = 1, value = 0.05)),
                column(
                  width = 12,
                  sliderInput("QValueGSEA", "Select a Q-Value:", min = 0, max = 1, value = 0.05)),
                column(
                  width = 12,
                  checkboxInput("metricAbsoluteValGSEA", "Absolute value", value = TRUE))))))),
      
      ### Pathway Enrichment : ORA ###
      tabItem(
        tabName = "pathway_ORA_subtab",
        h2("Pathway Enrichment : ORA", style = "text-align: center"),
        
        fluidRow(
          # PATHWAY SETTINGS
          column(
            width = 8,
            box(
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DBSelection", "Select a Database", 
                               choices = c("KEGG" = "kegg", "REACTOME" = "reactome")))))),
          
          # DEG Selection
          column(
            width = 8,
            box(
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DEGSelectionORA", "Select the correct DEG processing:", 
                               choices = c("Over expressed DEG only" = "OverDEGORA", "Under expressed DEG only" = "UnderDEGORA", "Both" = "BothDEG")))))),
          
          # Adjusted P-Value Cutoff
          column(
            width = 8,
            box(
              title = "Adjusted P-Value Cutoff",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  sliderInput("adjustedPValueCutoffORA", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05))))))),
      
      ### Pathway Enrichment : GSEA ###
      tabItem(
        tabName = "pathway_GSEA_subtab",
        h2("Pathway Enrichment : GSEA", style = "text-align: center"),
        
        fluidRow(
          column(
            width = 8,
            box(
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DBSelection", "Select a Database", 
                               choices = c("KEGG" = "kegg", "REACTOME" = "reactome")))))),
          
          # Adjusted P-Value Cutoff
          column(
            width = 8,
            box(
              title = "Adjusted P-Value Cutoff",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  sliderInput("adjustedPValueCutoffGSEA", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05))))))))))



################################################################################
################################################################################
################################################################################