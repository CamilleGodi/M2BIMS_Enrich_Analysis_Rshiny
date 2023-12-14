# Rshiny application produced by Victor BAILLEUL (Université de Rouen Normandie)
# Last updated: 24/11/2023

# This application aimed to facilitate functional enrichment analysis from
# differential expression results

# Contact me at: victor.bailleul@univ-rouen.fr

source("global.R")
source("./utils/inspection.R")

################################################################################
################################################################################
################################################################################

function(input, output, session) {
  
  ### Storage of the .csv input file date in a "reactive" object (if the .csv is valid) ###
  reactiveDataExpDiff <- reactive({
    req(input$inputFile)
    inFile <- input$inputFile
    ext <- tools::file_ext(inFile$datapath) # Récupération de l'extension
    data <- data.table::fread(inFile$datapath, header = TRUE) # csv reading
    
    # Check the input file extension and the column names
    required_columns <- c('GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj')
    
    if (ext != "csv") {
      shinyalertWrapper(title = "Error: The uploaded file doesn't have the right extension (.csv)",
                        message = "",
                        type = "error")
    } else if (!all(required_columns %in% colnames(data))) {
      shinyalertWrapper(title = "Error: Incorrect columns in file",
                        message = "Expected columns : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'",
                        type = "error")
    }
    
    shiny::validate(
      need(ext == "csv", "Incorrect file type"),
      need(sum(!colnames(data) %in% required_columns) == 0, "Incorrect columns in file")
    )
    
    return(data)
  })
  
  ### Check the input file extension (.csv) as soon as the file is uploaded ###
  observeEvent(input$inputFile, {
    reactiveDataExpDiff()
  })
  
  ### Management of the preview of the filtered data table [Whole data inspection] ###
  output$dataPreview <- DT::renderDT({
    filtered_data <- filter_dt(reactiveDataExpDiff(),
                      fc_cutoff = as.numeric(input$fc_cutoff),
                      padj_cutoff = as.numeric(input$padj_cutoff)
                      )
    
    
    # Create Ensembl link for genes with an Ensembl ID
    filtered_data$ID <- createEnsemblHTMLlink(input$selectOrganism, filtered_data$ID)
    
    # Preview of the filtered data table ( "escape = FALSE" allows HTML formatting )
    DT::datatable(filtered_data, options = list(scrollX = TRUE, pageLength = 25), escape = FALSE)
  })
  
  ### Management of the download of the filtered data table [Whole data inspection] ###
  output$downloadFilteredTable <- downloadHandler(
    filename = function() {
      paste("filtered-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      filtered_data <- filter_dt(reactiveDataExpDiff(),
                        fc_cutoff = as.numeric(input$fc_cutoff),
                        padj_cutoff = as.numeric(input$padj_cutoff)
                        )
      
      # Create Ensembl link for genes with an Ensembl ID
      filtered_data$EnsemblLink <- createEnsemblLink(input$selectOrganism, filtered_data$ID)
      
      
      
      # Write the filtered data to the specified file
      write.csv(filtered_data, file, row.names = FALSE)
    }
  )
  
  
  ### Management of the volcano plot [Whole data inspection] ###
  
  ranges_volca_plot <- reactiveValues(x = NULL, y = NULL)

  output$volcanoPlot <- renderPlotly({
    if(!is.null(reactiveDataExpDiff())){
      draw_volcano(reactiveDataExpDiff(),
                   fc_cutoff = as.numeric(input$fc_cutoff),
                   padj_cutoff = as.numeric(input$padj_cutoff),
                   xlim = ranges_volca_plot$x,
                   ylim = ranges_volca_plot$y,
                   lines = TRUE)
    }
    
  })
}


################################################################################
################################################################################
################################################################################

