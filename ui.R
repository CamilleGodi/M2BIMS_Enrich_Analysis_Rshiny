# Rshiny application created by 
# Victor BAILLEUL 
# Camille GODI
# Benjamin MARSAC
# Komlan Dieu-Donné TOTO
# Affiliation : Université de Rouen Normandie

# This application facilitates functional enrichment analysis from
# differential expression results

source("./utils/global.R")

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
        c("", "Arabidopsis thaliana", "Escherichia coli (K12)", "Homo sapiens", "Mus musculus", "Saccharomyces cerevisiae"),
        selected = "",
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
          h4("1/ Select a CSV (or CSV2) file. It must have the following columns and no others : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'."),
          h4("2/ Select the scientific name of the organism from which the data originates. Tip : you can type to search in the box."),
          h4("3/ (Optionnal) Explore your data through the 'Whole data inspection' tab"),
          h4("4/ Perform desired analysis through the appropriate tab."),
          br(),
          h4("Note : a plot can be downloaded by doing right-click > 'save image as'.")
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
          p("This app is the result of a group work in second year of Bioinformatics Master's Degree, 'BIMS', year 2023-2024s, Université de Rouen Normandie ( URN )."),
          br(),
          strong("Acknowledgments :"),
          p("We thank Solène Pety and Hélène Dauchel for their guidance and advices.")
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
          column(
            width = 6,
            box(
              title = "GENE ONTOLOGY SETTINGS",
              status = "warning",
              solidHeader = TRUE,
              width = 12,  
              fluidRow(
                # GENE ONTOLOGY SETTINGS
                column(
                  width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      checkboxGroupInput("goAnnotationORA", "Select a GO Annotation:",
                                         choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"),
                                         selected = "BP")))))),
            
            # GO Level Selection
            box(
              title = "GO Level Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              column(
                width = 12,
                fluidRow(
                  column(
                    width = 12,
                    radioButtons("levelSelectionORA", "Select a GO Level:", 
                                 choices = c("One-level GO ORA" = "OneLevelORA", "All-level GO ORA" = "AllLevelORA"))),
                  # Conditional Panel for OneLevelORA
                  conditionalPanel(
                    condition = "input.levelSelectionORA == 'OneLevelORA'",
                    fluidRow(
                      column(
                        width = 12,
                        fluidRow(
                          column(
                            width = 12,
                            sliderInput("levelSliderORA", "Select a GO level:", min = 1, max = 7, value = 1)))))))))),
          
          # DEG Selection
          column(
            width = 6,
            box(
              title = "DEG profile Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              column(
                width = 12,
                fluidRow(
                  column(
                    width = 12,
                    radioButtons("DEGSelectionORAGo", "Select the correct DEG processing:", 
                                 choices = c("Over expressed DEG only" = "OverDEG", "Under expressed DEG only" = "UnderDEG", "Both" = "BothDEG")))))),
            
            # Parameters
            box(
              title = "Parameters",
              status = "warning",
              solidHeader = TRUE,
              width = 12, 
              fluidRow(
                column(
                  width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      sliderInput("PValueORA", "Select a P-Value:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      sliderInput("QValueORA", "Select a Q-Value:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      sliderInput("adjustedPValueCutoffORA", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05)))))))),
      
        # ORA GO OUTPUTS
        fluidRow(
          box(width = 10,
              plotOutput("ORAgoDotPlot", height="600px"),
             ),
          box(width = 10,
              plotOutput("ORAgoCNETPlot", height="800px"),
          ),
          box(width = 10,
              plotOutput("ORAgoTreePlot", height="800px"),
          ),
          box(width = 10,
              plotOutput("ORAgoEmapPlot", height="800px"),
          ),
          
          box(width =12,
              DTOutput("results_ora_go_preview_table", height = "1000px")
          )
        )
      ),
      
      ### GO Term Enrichment : GSEA ###
      tabItem(
        tabName = "GO_term_GSEA_subtab",
        h2("GO Term Enrichment : GSEA", style = "text-align: center"),
        
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "GENE ONTOLOGY SETTINGS",
              status = "warning",
              solidHeader = TRUE,
              width = 12,    
              fluidRow(
                # GENE ONTOLOGY SETTINGS
                column(
                  width = 8,
                  fluidRow(
                    column(
                      width = 12,
                      checkboxGroupInput("goAnnotationGSEA", "Select a GO Annotation:",
                                         choices = c("Biological Process" = "BP", "Molecular Function" = "MF", "Cellular Component" = "CC"),
                                         selected = "BP")))))),
            
            # GO Level Selection
            box(
              title = "GO Level Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,  
              column(
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
                      fluidRow(
                        column(
                          width = 12,
                          sliderInput("levelSliderGSEA", "Select a GO level:", min = 1, max = 7, value = 1))))))))),
          
          # DEG Selection
          
          column(
            width = 6,
            box(
              title = "DEG profile Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              column(
                width = 8,
                fluidRow(
                  column(
                    width = 12,
                    radioButtons("DEGSelectionORAPathways", "Select the correct DEG processing:", 
                                 choices = c("Over expressed DEG only" = "OverDEG", "Under expressed DEG only" = "UnderDEG", "Both" = "BothDEG")))))),
            
            # Parameters
            
            box(
              title = "Parameters",
              status = "warning",
              solidHeader = TRUE,
              width = 12, 
              fluidRow(
                column(
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
                      sliderInput("adjustedPValueCutoffGSEA", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      checkboxInput("metricAbsoluteValGSEA", "Absolute value", value = TRUE))))))))),
      
      ### Pathway Enrichment : ORA ###
      tabItem(
        tabName = "pathway_ORA_subtab",
        h2("Pathway Enrichment : ORA", style = "text-align: center"),
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Databases",
              status = "warning",
              solidHeader = TRUE,
              width = 12,    
              fluidRow(
                column(
                  width = 8,
                  fluidRow(
                    column(
                      width = 12,
                      radioButtons("DBSelectionORA", "Select a Database", 
                                   choices = c("KEGG" = "kegg", "REACTOME" = "reactome"), selected = "kegg"))))))),
          
          # DEG Selection
          column(
            width = 6,
            box(
              title = "DEG profile Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              column(
                width = 12,
                fluidRow(
                  column(
                    width = 12,
                    radioButtons("DEGSelectionORA", "Select the correct DEG processing:", 
                                 choices = c("Over expressed DEG only" = "OverDEGORA", "Under expressed DEG only" = "UnderDEGORA", "Both" = "BothDEG")))))))),
        
        #parameters
        fluidRow(
          fluidRow(
            column(
              width = 12,
              fluidRow(
                column(
                  width = 12,
                  box(
                    title = "Parameters",
                    status = "warning",
                    solidHeader = TRUE,
                    width = 12,
                    fluidRow(
                      column(
                        width = 12,
                        sliderInput("PValueORAPathways", "Select a P-Value:", min = 0, max = 1, value = 0.05)),
                      column(
                        width = 12,
                        sliderInput("QValueORAPathways", "Select a Q-Value:", min = 0, max = 1, value = 0.05)),
                      column(
                        width = 12,
                        sliderInput("adjustedPValueCutoffORAPathways", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05))))))))),
      
      # ORA KEGG OUTPUTS
      fluidRow(
        box(width = 10,
            plotOutput("ORAPathwaysDotPlot", height="600px"),
        ),
        box(width = 10,
            plotOutput("ORAPathwaysCNETPlot", height="800px"),
        ),
        box(width = 10,
            plotOutput("ORAPathwaysTreePlot", height="800px"),
        ),
        box(width = 10,
            plotOutput("ORAPathwaysEmapPlot", height="800px"),
        ),
        
        box(width =12,
            DTOutput("results_ora_pathways_preview_table", height = "1000px")
        )
      )
    ),
    
      ### Pathway Enrichment : GSEA ###
      tabItem(
        tabName = "pathway_GSEA_subtab",
        h2("Pathway Enrichment : GSEA", style = "text-align: center"),
        
        fluidRow(
          column(
            width = 6,
            box(
              title = "Databases",
              status = "warning",
              solidHeader = TRUE,
              width = 12,    
              fluidRow(
                column(
                  width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      radioButtons("DBSelectionGSEA", "Select a Database", 
                                   choices = c("KEGG" = "kegg", "REACTOME" = "reactome"), selected = "kegg"))))))),
          
          
          column(
            width = 6,
            box(
              title = "DEG profile Selection",
              status = "warning",
              solidHeader = TRUE,
              width = 12,
              column(
                width = 12,
                fluidRow(
                  column(
                    width = 12,
                    radioButtons("DEGSelectionGSEAPathways", "Select the correct DEG processing:", 
                                 choices = c("Over expressed DEG only" = "OverDEGGSEA", "Under expressed DEG only" = "UnderDEGGSEA", "Both" = "BothDEG")))))))),
        
        # Parameters
        fluidRow(
          fluidRow(
            column(
              width = 12,
              column(
                width = 12,
                box(
                  title = "Parameters",
                  status = "warning",
                  solidHeader = TRUE,
                  width = 12,
                  fluidRow(
                    column(
                      width = 12,
                      sliderInput("PValueCutoffGSEAPathways", "Select a P-Value:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      sliderInput("QValueGSEAPathways", "Select a Q-Value:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      sliderInput("adjustedPValueCutoffGSEAPathways", "Select an adjusted P-Value Cutoff:", min = 0, max = 1, value = 0.05)),
                    column(
                      width = 12,
                      checkboxInput("metricAbsoluteValGSEAPathways", "Absolute value", value = TRUE))))))))))))


################################################################################
################################################################################
################################################################################