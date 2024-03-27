# Rshiny application created by 
# Victor BAILLEUL 
# Camille GODI
# Benjamin MARSAC
# Komlan Dieu-Donné TOTO
# Affiliation : Université de Rouen Normandie

# This application facilitates functional enrichment analysis from
# differential expression results

source("./utils/global.R")
source("./utils/inspection.R")

################################################################################
################################################################################
################################################################################

function(input, output, session) {
  ### Storage of the .csv input file date in a "reactive" object (if the .csv is valid) ###
  reactive_data_expr_diff <- reactive({
    req(input$input_file)
    input_file <- input$input_file
    ext <- tools::file_ext(input_file$datapath) # Récupération de l'extension
    data <- data.table::fread(input_file$datapath, header = TRUE) # csv reading
    
    # Check the input file extension and the column names
    required_columns <- c('GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj')
    
    if (ext != "csv") {
      shinyalert_wrapper(title = "Error: The uploaded file doesn't have the right extension (.csv)",
                        message = "",
                        type = "error")
    } else if (!all(required_columns %in% colnames(data)) | length(colnames(data)) > 6) {
      shinyalert_wrapper(title = "Error: Incorrect columns in file",
                        message = "Expected columns : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj' (and no other ones)",
                        type = "error")
    }
    
    shiny::validate(
      need(ext == "csv", "Incorrect file type"),
      need(sum(!colnames(data) %in% required_columns) == 0, "Incorrect columns in file. Expected columns : 'GeneName', 'ID', 'baseMean', 'log2FC', 'pval', 'padj'  (and no other ones)")
    )
    
    return(data)
  })
  
  ### Check the input file extension (.csv) as soon as the file is uploaded ###
  observeEvent(input$input_file, {
    reactive_data_expr_diff()
  })
  
  organism_library_go <- reactive({organism_conversion_table[input$select_organism, "annotation_db"]})

  reactive_data_new_ids <- reactive({prepare_pipe(reactive_data_expr_diff(),
                                                  organism_db = organism_library_go(),
                                                  from = "ENSEMBL")})

  filtered_data <- reactive({filter_dt(reactive_data_new_ids(),
                                       fc_cutoff = as.numeric(input$fc_cutoff),
                                       padj_cutoff = as.numeric(input$padj_cutoff))
  })
  
  ### Management of the preview of the filtered data table [Whole data inspection] ###
  output$data_preview_table <- DT::renderDT({
    filtered_data <- show_filtered_df(filtered_data())
    
    # Create Ensembl link for genes with an Ensembl ID
    filtered_data$ID <- create_Ensembl_html_link(input$select_organism, filtered_data$ID)
    
    # Preview of the filtered data table ( "escape = FALSE" allows HTML formatting )
    DT::datatable(filtered_data, options = list(scrollX = TRUE, pageLength = 25), escape = FALSE)
  })
  
  ### Management of the download of the filtered data table [Whole data inspection] ###
  output$download_filtered_table <- downloadHandler(
    filename = function() {
      paste("filtered-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the filtered data to the specified file
      write.csv(show_filtered_df(filtered_data()), file, row.names = FALSE)
    }
  )
  
  
  ### Management of the volcano plot [Whole data inspection] ###
  ranges_volca_plot <- reactiveValues(x = NULL, y = NULL)

  output$volcano_plot <- renderPlotly( {
    if(!is.null(reactive_data_expr_diff())){
      draw_volcano(filtered_data(),
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

