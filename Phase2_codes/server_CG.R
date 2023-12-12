################################################################################
####################### M2.1 BIMS - PROJET R SHINY #############################
##################### Camille GODI - Promo 2022-2025 ###########################
################################################################################

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinyalert)
library(data.table)
library(ggplot2)
library(ggrepel)

################################################################################

createEnsemblLink <- function(organism, gene_id) {
  sprintf("https://www.ensembl.org/%s/Gene/Summary?g=%s", organism, gene_id)
}

createEnsemblHTMLlink <- function(organism, gene_id) {
  sprintf('<a href="https://www.ensembl.org/%s/Gene/Summary?g=%s">%s</a>', organism, gene_id, gene_id)
}

checkColumnNames <- function(colNames) {
  # expected column names (format DeSEQ2) : 
  # geneName (added a bit after DeSEQ2 analysis), geneID (same), baseMean, log2FoldChange, pvalue, padj
  for ( expectedColumn in list("geneName", "geneID", "baseMean", "log2FoldChange", "pvalue", "padj") ) {
    if ( !( expectedColumn %in% colNames) ) {
      return(FALSE)
    }
  } 
  return(TRUE)
}


InteractiveVolcanoPlot <- function(data, p_value, FCthreshold, dynamic_ranges){
  
  colours <- c("hotpink", "blue") ; names(colours) = c("upregulated", "downregulated")
  
  # Create diff expression categories -> colour later
  data$Differentially_Expressed_Genes <- "no"
  data$Differentially_Expressed_Genes[data$log2FoldChange > FCthreshold & data$padj < p_value] <- "upregulated"
  data$Differentially_Expressed_Genes[data$log2FoldChange < - FCthreshold & data$padj < p_value] <- "downregulated"
  
  # Add label for DE genes
  data$DElabel <- NA
  data$DElabel[data$Differentially_Expressed_Genes != "no"] <- data$geneName[data$Differentially_Expressed_Genes != "no"]
  
  vplot <- ggplot(data = data, aes(x = log2FoldChange,  
                                   y = -log10(padj), 
                                   col = Differentially_Expressed_Genes,
                                   label = DElabel)) +
    ggtitle(sprintf("Volcano plot of the provided data ( thresholds : padj = %g, log2FC = %g )", p_value, FCthreshold) ) +
    theme(legend.position = "bottom", legend.text = element_text(size=10) ) +
    labs(color = NULL) +
    coord_cartesian(xlim = dynamic_ranges$x, ylim = dynamic_ranges$y, expand = FALSE) +
    geom_point(size = 0.5) +
    geom_text_repel(size = 4, force = 1, 
                    max.overlaps = 40,
                    min.segment.length = 0.2, show.legend = FALSE) +
    scale_colour_manual(values = colours) +
    geom_vline(xintercept = c(-FCthreshold, FCthreshold), col = "black", size = 0.2) +
    geom_hline(yintercept = -log10(p_value), col = "red", size = 0.2)
  vplot
}


################################################################################

function(input, output) {
  
  # Allow for file input after organism only (WIP: not working yet)
  observeEvent(input$organism != "" , enable("csv_upload") )
  
  # If uploaded file is not a .csv, create an error pop-up ...
  observeEvent({ input$csv_upload } ,
      {if ( !grepl(".csv", input$csv_upload$name, fixed=TRUE) ) {
          shinyalert(
          title = "Please upload a csv file",
          text = "",
          size = "xs", 
          closeOnEsc = TRUE,
          closeOnClickOutside = TRUE,
          type = "error",
          showConfirmButton = FALSE,
          showCancelButton = TRUE,
          cancelButtonText = "Cancel",
          animation = TRUE
          )

      # ... Else, load it and create a renderDataTable object for the "Table" box if column names are alright
      } else {  
        
        csv_contents <- fread(input$csv_upload$datapath, header =TRUE)
        csv_contents_copy <- csv_contents          # copy for downloadable table
        
        if ( !checkColumnNames(colnames(csv_contents)) ){
          shinyalert(
            title = "Please upload a csv file with the right column names",
            text = "expected column names ( derived from DeSEQ2 format ) :\ngeneName, geneID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj",
            size = "xs", 
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            type = "error",
            showConfirmButton = FALSE,
            showCancelButton = TRUE,
            cancelButtonText = "Cancel",
            animation = TRUE
          )
        } else {
        
          # Get index of the adjusted p-value column, to order automatically by it
          index_padj <- (grep("^padj$", colnames(csv_contents)) - 1 )
          
          # Add Ensembl link for each genes, in place of the ID for the working table, in a new column for the downloadable table
          csv_contents$geneID <- createEnsemblHTMLlink(input$organism, csv_contents$geneID)
          csv_contents_copy$geneEnsemblLinks <- createEnsemblLink(input$organism, csv_contents_copy$geneID)
          
          # Enable download of the table
          enable("download_table")
          output$download_table <- downloadHandler(
            filename = function() { "DE_genes_table.csv" },
            content = function(file) {
              write.csv(csv_contents_copy, file = file, row.names = FALSE)
            })
          
          # Create data (table) visualisation 
          render_data_options <- list( pageLength = 10, 
                                       ordering = TRUE, 
                                       order = list(index_padj, "asc")
                                       )
          output$input_csv_table <- renderDataTable( csv_contents, options = render_data_options, escape = FALSE )
          
          # Interactive ranges for zoom
          volcano_plot_ranges <- reactiveValues(x = NULL, y = NULL)
          
          # Create volcano Plot + downloader
          volcano_plot <- reactive(InteractiveVolcanoPlot(csv_contents, input$p_value_slider, input$log2_slider, volcano_plot_ranges ))
          output$volcano_plot <- renderPlot(volcano_plot())
          enable("download_volcano")
          output$download_volcano <- downloadHandler(
            filename = function() { "volcano_plot.png" },
            content = function(file) {
              png( file, width = 800, height = 600 )
              plot( volcano_plot() )
              dev.off()
            })
          
          # When a double-click happens, check if there's a brush on the plot.
          # If so, zoom to the brush bounds; if not, reset the zoom.
          observeEvent(input$volcano_plot_dblclick, {
            brush <- input$volcano_plot_brush
            if (!is.null(brush)) {
              volcano_plot_ranges$x <- c(brush$xmin, brush$xmax)
              volcano_plot_ranges$y <- c(brush$ymin, brush$ymax)
              
            } else {
              volcano_plot_ranges$x <- NULL
              volcano_plot_ranges$y <- NULL
            }
          })
        }
      }
    }
  )
}
