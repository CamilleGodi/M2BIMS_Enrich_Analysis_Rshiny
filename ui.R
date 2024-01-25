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
        width = '100%'),
      
      selectInput(
        "selectOrganism",
        "Select organism name:",
        c("Arabidopsis thaliana", "Caenorhabditis elegans","Danio rerio", "Drosophila melanogaster","Homo sapiens", "Mus musculus", "Saccharomyces cervisiae", "Xenopus laevis"),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL),
      
      sidebarMenu(
        menuItem("Whole Data Inspection", tabName = "WholeData", icon = icon("table")),
        menuItem("GO Term Enrichment", tabName = "GOTermEnrichment", icon = icon("sitemap"),
                 menuSubItem("ORA", tabName = "GOTermORA"),
                 menuSubItem("GSEA", tabName = "GOTermGSEA")),
        menuItem("Pathway Enrichment", tabName = "PathwayEnrichment", icon = icon("project-diagram"),
                 menuSubItem("ORA", tabName = "pathORA"),
                 menuSubItem("GSEA", tabName = "pathGSEA")),
        menuItem("About", tabName = "About", icon = icon("info-circle"))))),
  
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
        tabName = "Home",
        h2("Home", style = "text-align: center")),
      
      ### Whole Data Inspection ###
      tabItem(
        tabName = "WholeData",
        h2("Whole data inspection", style = "text-align: center"),
        fluidRow(
          box(
            title = "Volcano Plot",
            status = "warning", # Cosmetic purpose only
            solidHeader = TRUE, 
            collapsible = FALSE,
            plotlyOutput("volcanoPlot", height = "400px"),
            width = 6),
          box(
            title = "Options",
            status = "warning", # Cosmetic purpose only
            solidHeader = TRUE, 
            collapsible = FALSE,
            sliderInput("pValueCutoff", "P-Value cutoff", min = 0, max = 1, value = 0.05),
            sliderInput("logFC", "log2 FoldChange cutoff", min = 0, max = 10, value = 1),
            downloadButton("downloadFilteredTable", label = "Download filtered data table", icon = shiny::icon("download")),
            width = 6)),
        fluidRow(
          box(
            title = "Filtered table preview",
            status = "warning", # Cosmetic purpose only 
            DTOutput("dataPreview"),
            width = 12,
            solidHeader = TRUE))),
      
      ### GO Term Enrichment : ORA ###
      tabItem(
        tabName = "GOTermORA",
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
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  radioButtons("DEGSelectionORA", "Select the correct DEG processing:", 
                               choices = c("Over expressed DEG only" = "OverDEGORA", "Under expressed DEG only" = "UnderDEGORA", "Both" = "BothDEG"))))))),
        
        # Adjusted P-Value Cutoff
        fluidRow(
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
      ### GO Term Enrichment : GSEA ###
      tabItem(
        tabName = "GOTermGSEA",
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
                  sliderInput("adjustedPValueCutoffGSEA", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05))))))),
      
      ### Pathway Enrichment : ORA ###
      tabItem(
        tabName = "pathORA",
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
        tabName = "pathGSEA",
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
