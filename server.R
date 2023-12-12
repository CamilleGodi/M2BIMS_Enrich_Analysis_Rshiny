# Rshiny application produced by Victor BAILLEUL (Université de Rouen Normandie)
# Last updated: 24/11/2023

# This application aimed to facilitate functional enrichment analysis from
# differential expression results

# Contact me at: victor.bailleul@univ-rouen.fr

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)       # Utilisé pour l'affichage de tableaux
library(ggplot2) 
library(plotly)   # Utilisé pour le plot interactif

################################################################################
################################################################################
################################################################################

createEnsemblLink <- function(organism, gene_id) {
  sprintf("https://www.ensembl.org/%s/Gene/Summary?g=%s", organism, gene_id)
}

createEnsemblHTMLlink <- function(organism, gene_id) {
  sprintf('<a href="https://www.ensembl.org/%s/Gene/Summary?g=%s">%s</a>', organism, gene_id, gene_id)
}

################################################################################
################################################################################
################################################################################


function(input, output, session) {
  
  ### Stockage des données du .csv dans un objet "reactive" (si .csv valide) ###
  reactiveDataExpDiff <- reactive({
    req(input$inputFile)
    inFile <- input$inputFile
    ext <- tools::file_ext(inFile$datapath) # Récupération de l'extension
    data <- read.csv(inFile$datapath, sep = ";", header = TRUE, stringsAsFactors = FALSE) # Lecture du csv
    
    # Verification de l'extension et des colonnes du csv
    required_columns <- c('GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj')
    
    if (ext != "csv") {
      showNotification("Error: The uploaded file doesn't have the extension (.csv)", type = "error", duration = NULL)
    } else if (!all(required_columns %in% colnames(data))) {
      showNotification("Error: The uploaded file does not contain the required columns or use the wrong separator (only ';' accepted)", type = "error", duration = NULL)
    }
    
    shiny::validate(
      need(ext == "csv", "Incorrect file type"),
      need(sum(!colnames(data) %in% required_columns) == 0, "Incorrect columns in file")
    )
    
    return(data)
  })
  
  ### Vérifier le format du .csv dès l'upload du fichier ###
  observeEvent(input$inputFile, {
    reactiveDataExpDiff()
  })
  
  ### Gestion de l'affichage du tableau filtré [Whole data inspection] ###
  output$dataPreview <- renderDT({
    data <- reactiveDataExpDiff()
    
    # On filtre les valeurs par rapport aux valeurs des sliders
    cutoff_logFC <- input$logFC
    cutoff_padj <- input$pValueCutoff
    
    data$highlight <- ifelse(
      data$log2FC <= -cutoff_logFC & data$padj <= cutoff_padj, 'Underexpressed',
      ifelse(data$log2FC >= cutoff_logFC & data$padj <= cutoff_padj, 'Overexpressed', 'Not Highlighted')
    )
    filtered_data <- data[data$highlight != 'Not Highlighted', ]
    
    # Affichage de l'aperçu
    DT::datatable(filtered_data, options = list(pageLength = 5))
  })
  
  ### Gestion du téléchargement du tableau filtré [Whole data inspection] ###
  output$downloadFilteredTable <- downloadHandler(
    filename = function() {
      paste("filtered-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      data <- reactiveDataExpDiff()
      
      # On filtre les valeurs par rapport aux valeurs des sliders
      cutoff_logFC <- input$logFC
      cutoff_padj <- input$pValueCutoff
      
      data$highlight <- ifelse(
        data$log2FC <= -cutoff_logFC & data$padj <= cutoff_padj, 'Underexpressed',
        ifelse(data$log2FC >= cutoff_logFC & data$padj <= cutoff_padj, 'Overexpressed', 'Not Highlighted')
      )
      
      filtered_data <- data[data$highlight != 'Not Highlighted', ]
      
      # Write the filtered data to the specified file
      write.csv(filtered_data, file, row.names = FALSE)
    }
  )
  
  ### Gestion du volcano plot [Whole data inspection] ###
  output$volcanoPlot <- renderPlotly({
    data <- reactiveDataExpDiff()
    
    # Récupération des cutoff via le slider
    cutoff_logFC <- input$logFC
    cutoff_padj <- input$pValueCutoff
    
    # Ajout d'une colonne sur le tableau triant les points en trois classes
    # en fonction des cutoffs sur la padj et le log2FC choisi
    data$highlight <- ifelse(
      data$log2FC <= -cutoff_logFC & data$padj <= cutoff_padj, 'Underexpressed',
      ifelse(data$log2FC >= cutoff_logFC & data$padj <= cutoff_padj, 'Overexpressed', 'Not Highlighted')
    )
    
    # Creation du plot
    p <- ggplot(data, aes(x = log2FC, y = -log10(padj), color = highlight)) +
      geom_point(alpha = 0.8) +
      scale_color_manual(values = c('Underexpressed' = 'blue', 'Overexpressed' = 'red', 'Not Highlighted' = 'grey')) +
      labs(title = 'Volcano Plot', x = 'Log2 Fold Change', y = '-Log10 p-adjusted') +
      theme_minimal()
    
    # Interactivité du plot
    ggplotly(p, tooltip = c("text"), dynamicTicks = TRUE)
  })
}


################################################################################
################################################################################
################################################################################

