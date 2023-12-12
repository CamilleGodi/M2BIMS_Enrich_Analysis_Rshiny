#' @title server.R
#' @author Benjamin Marsac
#' @description Fichier de serveur pour interface Rshiny
#' @date 2023-2024

source("./global.R")

server = function(input, output) {
  ranges_MA_plot = reactiveValues(x = NULL, y = NULL)
  ranges_volca_plot = reactiveValues(x = NULL, y = NULL)
  
  
  df_products_upload = reactive({
    in_file = input$select_csv

    if (is.null(in_file))
      return(NULL)
    
    # on check si le fichier est bien au format csv (sinon erreur précisant le format), 
    # le if est juste là pour avoir le jalon rouge d'erreur en plus de l'information dans le box du DT
    ext = tools::file_ext(in_file$datapath)
    req(in_file)
    if(ext != "csv"){                       
      shiny::showNotification("Message text",
                              action ="Your file is not in the required format",
                              type = "error",
                              duration = 10)
      validate(need(ext == "csv", "Please upload a csv file"))
    }
    # on charge le fichier (fread trouve tout seul le séparateur normalement)
    
    df = data.table::fread(in_file$datapath, header = TRUE)
    
    required_colnames = c("GeneName",	"ID",	"baseMean",	"log2FC",	"pval",	"padj")
    
    if(sum(!colnames(df) %in% required_colnames) > 0 | length(colnames(df)) < 6){
      shiny::showNotification("Message text",
                              action ="You must have 6 columns named : GeneName ID baseMean log2FC pval padj",
                              type = "error",
                              duration = 10)
      validate(need(sum(!colnames(df) %in% required_colnames) > 0,
                    "Please check your column's name"))
      validate(need(length(colnames(df)) < 6,
                    "Please the number of column"))
    }
    return(df)
  })
  
  # Affichage du datatable scrollable avec par défaut 25gènes/pages, on n'affiche que 3 décimales car trop c'est moche (mais l'info n'est pas perdu)
  output$table = DT::renderDataTable({
    if(!is.null(df_products_upload())){
      DT::datatable(df_products_upload()[abs(df_products_upload()$log2FC) > as.numeric(input$slider_of_log2_foldchange_for_volcanoplot) &
                                                  df_products_upload()$padj < as.numeric(input$slider_of_pvalue_for_plots),],options = list(scrollX = TRUE,pageLength = 25))
    }
  })

  output$draw_volcanoplot = renderPlot({
      if(!is.null(df_products_upload())){
        draw_volcano(df_products_upload(),
                     fc_cutoff = as.numeric(input$slider_of_log2_foldchange_for_volcanoplot),
                     alpha = as.numeric(input$slider_of_pvalue_for_plots),
                     xlim = ranges_volca_plot$x,
                     ylim = ranges_volca_plot$y)
      }
    })

  output$draw_maplot = renderPlot({
    if(!is.null(df_products_upload())){
      draw_MA(df_products_upload(),alpha =as.numeric(input$slider_of_pvalue_for_plots),xlim = ranges_MA_plot$x, ylim = ranges_MA_plot$y)
    }
  })

  observeEvent(input$maplot_double_click, {
    brush <- input$maplot_brush
    if (!is.null(brush)) {
      ranges_MA_plot$x <- c(brush$xmin, brush$xmax)
      ranges_MA_plot$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_MA_plot$x <- NULL
      ranges_MA_plot$y <- NULL
    }
  })
  
  observeEvent(input$volcanoplot_double_click, {
    brush <- input$volcanoplot_brush
    if (!is.null(brush)) {
      ranges_volca_plot$x <- c(brush$xmin, brush$xmax)
      ranges_volca_plot$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges_volca_plot$x <- NULL
      ranges_volca_plot$y <- NULL
    }
  })
  
  output$download <- downloadHandler(
    filename = function(){"thename.csv"}, 
    content = function(fname){
      write.csv(df_products_upload()[abs(df_products_upload()$log2FC) > as.numeric(input$slider_of_log2_foldchange_for_volcanoplot) &
                                       df_products_upload()$padj < as.numeric(input$slider_of_pvalue_for_plots),], fname)
    }
  )
  
}

