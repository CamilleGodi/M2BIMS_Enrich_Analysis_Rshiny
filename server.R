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
source("./utils/ORA.R")
source("./utils/enrichment.R")
source("./utils/GSEA.R")

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
  
  organism_library_go <- reactive({ 
    if (!is.null(input$select_organism)) {
      organism_conversion_table[input$select_organism, "annotation_db"]
    } else {
      shinyalert_wrapper("Please choose the organism before attempting an analysis", type = "warning")
      NULL
    }})
  
  organism_kegg_code <- reactive({ 
    if (!is.null(input$select_organism)) {
      organism_conversion_table[input$select_organism, "kegg_name"]}
  })
  
  organism_reactome_name <- reactive({ 
    if (!is.null(input$select_organism)) {
      organism_conversion_table[input$select_organism, "reactome_name"]}
  })
  
  reactive_data_new_ids <- reactive({
    if (!is.null(organism_library_go())) {
      prepare_pipe(reactive_data_expr_diff(),
                  organism_db = organism_library_go(),
                  from = "ENSEMBL")
      }
    })
        
  
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
  
  
  ### Prepare "universe" for ORA/GSEA
  organism_universe<- reactive({ 
    if(!is.null(organism_library_go())){
      prepare_universe(filtered_data(), organism_library_go(), "ENSEMBL")
    }
  })
  
  
  ### ORA GO results ###
  results_ora_go <- reactive({ 
    if(!is.null(filtered_data())){
      res_ora_go <- do_ora_go_terms(
        filtered_data(),
        organism_library_go(),
        organism_universe()
      )
      return(res_ora_go)
    }
  })
  
  results_ora_go_filtered <- reactive({ 
    if(!is.null(results_ora_go())){
      
      ora_go_after_filter_ontologies <- filter_go_enrich_results(results_ora_go(), ontology = input$goAnnotationORA)
 
      res_ora_go_filtered <- filter_table_enrich_results(ora_go_after_filter_ontologies, 
                                  p_value_cutoff = input$PValueORA, 
                                  p_adj_cutoff   = input$adjustedPValueCutoffORA, 
                                  q_value_cutoff = input$QValueORA)
      return(res_ora_go_filtered)
    }
  })
  
  output$results_ora_go_preview_table <- DT::renderDT({
    preview_table <- results_ora_go_filtered() %>% enrich_pagination(alpha_cutoff = input$adjustedPValueCutoffORA)
    # Preview of the filtered data table ( "escape = FALSE" allows HTML formatting )
    DT::datatable(preview_table, options = list(scrollX = TRUE, pageLength = 25), escape = FALSE)
  })
  
  output$ORAgoDotPlot <- renderPlot({
    if(!is.null(results_ora_go_filtered())){
      results_ora_go_filtered() %>% draw_dotplot(
        show_category = 30, 
        title = paste("ORA - GO termes -", paste(input$goAnnotationORA, collapse="&"), "- Dot plot")
       )
    }
  })
  
  output$ORAgoCNETPlot <- renderPlot({
    if(!is.null(results_ora_go_filtered())){
      results_ora_go_filtered() %>% draw_cnetplot(
          category_label = 1,
          gene_list = results_ora_go()@result$geneID,
          category_color = "red",
          node_label = "category",
          title = paste("ORA - GO termes -", paste(input$goAnnotationORA, collapse="&"), "- CNET plot")
        )
    }
  })
  
  output$ORAgoTreePlot <- renderPlot({
    if(!is.null(results_ora_go_filtered())){
      results_ora_go_filtered() %>% draw_treeplot(
        gradient_col = c("red", "blue"),
        show_category = 30,
        n_cluster = 10,
        label_words_n = 4,
        h_clust_method = "ward.D2",
        title = paste("ORA - GO termes -", paste(input$goAnnotationORA, collapse="&"), "- Tree plot")
      )
    }
  })
  
  output$ORAgoEmapPlot <- renderPlot({
    if(!is.null(results_ora_go_filtered())){
      results_ora_go_filtered() %>% draw_emapplot(
        show_category = 10,
        category_label = 1,
        title = paste("ORA - GO termes -", paste(input$goAnnotationORA, collapse="&"), "- EMAP plot")
      )
    }
  })
  
  output$ORAgoBarPlot <- renderPlot({
    if(!is.null(results_ora_go_filtered())){
      results_ora_go_filtered() %>% barplot(
        title = paste("ORA - GO termes -", paste(input$goAnnotationORA, collapse="&"), "- Barplot")
      )
    }
  })
  
  ### Management of the download of the filtered data table [ORA GO TERMS] ###
  output$download_ora_go <- downloadHandler(
    filename = function() {
      paste("ora-go-", paste(input$goAnnotationORA, collapse="-"), "-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the filtered data to the specified file
      write.csv(results_ora_go(), file, row.names = FALSE)
    }
  )
  
  output$download_ora_go_filtered <- downloadHandler(
    filename = function() {
      paste("ora-go-", paste(input$goAnnotationORA, collapse="-"), "-filtered-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the filtered data to the specified file
      write.csv(results_ora_go_filtered(), file, row.names = FALSE)
    }
  )
  
  
  ### ORA KEGG/REACTOME results ###
  results_ora_pathways <- reactive({ 
    if (!is.null(filtered_data()) & input$DBSelectionORA == "kegg") {
      res_ora_pathways <- do_ora_kegg(
        filtered_data(),
        organism_library_go(),
        organism_universe(),
        organism_kegg_code()
      )
    } else if (!is.null(filtered_data()) & input$DBSelectionORA == "reactome") {
      if (is.na(organism_reactome_name()) ) {
        shinyalert_wrapper(title = "Error: organism doesn't have Reactome data available. Please use KEGG instead.",
                           message = "",
                           type = "error")
      } else {
        res_ora_pathways <- do_ora_reactome(
          filtered_data(),
          organism_reactome_name(),
          organism_universe()
        )
        print(res_ora_pathways)
      }
      
      return(res_ora_pathways)
    }
  })

  
  results_ora_pathways_filtered <- reactive({ 
    if(!is.null(results_ora_pathways())){
      res_ora_pathways_filtered <- filter_table_enrich_results(results_ora_pathways(), 
                                                         p_value_cutoff = input$PValueORAPathways, 
                                                         p_adj_cutoff   = input$adjustedPValueCutoffORAPathways, 
                                                         q_value_cutoff = input$QValueORAPathways
                                                         )
      return(res_ora_pathways_filtered)
    }
  })
  
  output$results_ora_pathways_preview_table <- DT::renderDT({
    preview_table <- results_ora_pathways_filtered() %>% enrich_pagination(alpha_cutoff = input$adjustedPValueCutoffORAPathways)
    # Preview of the filtered data table ( "escape = FALSE" allows HTML formatting )
    DT::datatable(preview_table, options = list(scrollX = TRUE, pageLength = 25), escape = FALSE)
  })
  
  output$ORAPathwaysDotPlot <- renderPlot({
    if(!is.null(results_ora_pathways_filtered())){
      results_ora_pathways_filtered() %>% draw_dotplot(
        show_category = 30, 
        title = paste("ORA -", input$DBSelectionORA, "pathways - Dot plot")
      )
    }
  })
  
  output$ORAPathwaysCNETPlot <- renderPlot({
    if(!is.null(results_ora_pathways_filtered())){
      results_ora_pathways_filtered() %>% draw_cnetplot(
        category_label = 1,
        gene_list = results_ora_pathways_filtered()@result$geneID,
        category_color = "red",
        node_label = "category",
        title = paste("ORA -", input$DBSelectionORA, "pathways - CNET plot")
      )
    }
  })
  
  output$ORAPathwaysTreePlot <- renderPlot({
    if(!is.null(results_ora_pathways_filtered())){
      results_ora_pathways_filtered() %>% draw_treeplot(
        gradient_col = c("red", "blue"),
        show_category = 30,
        n_cluster = 10,
        label_words_n = 4,
        h_clust_method = "ward.D2",
        title = paste("ORA -", input$DBSelectionORA, "pathways - Tree plot")
      )
    }
  })
  
  output$ORAPathwaysEmapPlot <- renderPlot({
    if(!is.null(results_ora_pathways_filtered())){
      results_ora_pathways_filtered() %>% draw_emapplot(
        show_category = 10,
        category_label = 1,
        title = paste("ORA -", input$DBSelectionORA, "pathways - EMAP plot")
      )
    }
  })
  
  output$ORAPathwaysBarPlot <- renderPlot({
    if(!is.null(results_ora_pathways_filtered())){
      results_ora_pathways_filtered() %>% barplot(
        title = paste("ORA -", input$DBSelectionORA, "pathways - Barplot")
      )
    }
  })
  
  ### Management of the download of the filtered data table [ORA GO TERMS] ###
  output$download_ora_pathway <- downloadHandler(
    filename = function() {
      paste("ora-", input$DBSelectionORA, "-pathways-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the filtered data to the specified file
      write.csv(results_ora_pathways(), file, row.names = FALSE)
    }
  )
  
  output$download_ora_pathway_filtered <- downloadHandler(
    filename = function() {
      paste("ora-", input$DBSelectionORA, "-pathways-filtered-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Write the filtered data to the specified file
      write.csv(results_ora_pathways_filtered(), file, row.names = FALSE)
    }
  )
  
  
  ### GSEA GO results ###
  results_gsea_go <- reactive({ 
    if(!is.null(filtered_data())){
      res_gsea_go <- do_gsea_go_terms(
        filtered_data(),
        organism_library_go(),
        reactive_abs = input$metricAbsoluteValGSEA
      )
      return(res_gsea_go)
    }
  })
  
  results_gsea_go_filtered <- reactive({ 
    if(!is.null(results_gsea_go())){
      
      res_gsea_go_filtered <- filter_table_enrich_results(results_gsea_go(), 
                                                          p_value_cutoff = input$PValueCutoffGSEA, 
                                                          p_adj_cutoff   = input$adjustedPValueCutoffGSEA, 
                                                          q_value_cutoff = input$QValueGSEA)
      
      gsea_go_after_filter_ontologies <- filter_go_enrich_results(res_gsea_go_filtered, ontology = input$goAnnotationGSEA)
      

      return(gsea_go_after_filter_ontologies)
    }
  })
  
  output$GSEAgoCNETPlot <- renderPlot({
    if(!is.null(results_gsea_go_filtered())){
      results_gsea_go_filtered() %>% draw_cnetplot(
        category_label = 1,
        gene_list = results_gsea_go()@result$geneID,
        category_color = "red",
        node_label = "category",
        title = paste("GSEA - GO termes -", paste(input$goAnnotationGSEA, collapse="&"), "- CNET plot")
      )
    }
  })
  
  
}


################################################################################
################################################################################
################################################################################

