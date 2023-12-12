# A shiny application for RNA-seq data analysis : server code
# Author: Komlan Dieu-Donné TOTO
#email: komlan-dieu-donne.toto@univ-rouen.fr
# Date: 24/11/2023

# Infos: To run this application, you need to install the following packages:
# shiny, shinydashboard, shinythemes, ggplot2, plotly, DT, data.table




library(data.table)
library(shiny)
library(shinyFiles)
library(plotly)
library(DT)
library(ggplot2)


# The generateVolcanoPlot function generates a volcano plot based on the provided data and user-defined thresholds.
# takes: data The input data frame containing the necessary variables for the plot.
#        Pvalue_cutoff_from_input The p-value cutoff threshold.
#        log2_FoldChange_cutoff_from_input The log2 fold change cutoff threshold.
# return A plotly volcano plot.

generateVolcanoPlot <- function(data, Pvalue_cutoff_from_input, log2_FoldChange_cutoff_from_input) {
   # I filter the data to keep only the points that are under or equal to the p-value or above the log2FC thresholds
  filtered_data <- data[(abs(data$pval) < Pvalue_cutoff_from_input | abs(data$log2FC) > log2_FoldChange_cutoff_from_input),] #abs allows to take the absolute value of the log2FC
  #  I create a vector of colors for the points
  color <- ifelse(data$pval < Pvalue_cutoff_from_input & abs(data$log2FC) > log2_FoldChange_cutoff_from_input, "blue", "red") # if the point is under the thresholds, it is blue, otherwise it is red
  color[!(data$pval < Pvalue_cutoff_from_input | abs(data$log2FC) > log2_FoldChange_cutoff_from_input)] <- "gray" # if the point is not under the p-value or above the log2FC thresholds, it is gray
  size <- rep(4, nrow(data))  # here I create a vector of size for the points
  size[data$pval < Pvalue_cutoff_from_input & abs(data$log2FC) > log2_FoldChange_cutoff_from_input] <- 4  # if the point is under the thresholds, its size is 4
  plot_ly( # i use the plot_ly function from the plotly package to create the plot
    filtered_data,
    x = ~log2FC,
    y = ~-log10(pval),
    type = 'scatter',
    mode = 'markers',
    marker = list(
      color = color,
      size = size,
      opacity = 1,   # points are not transparent
      line = list(width = 0, color = 'rgba(0,0,0,0)')  # points have no border, so they are not surrounded by a black line
    )
  ) %>%
    layout( # I use the layout function to add a title and axis titles
      title = "-log(pval) vs log2FC",
      xaxis = list(title = "Log2 Fold Change"),
      yaxis = list(title = "-Log10(P-value)"),
      hovermode = "closest"
    ) %>%
    config(displayModeBar = FALSE)   # I remove the toolbar that is displayed by default
}

server <- function(input, output) {
  # charging the data
  data <- reactive({
    req(input$file) # if no file is selected, the app will not crash
    file_data <- fread(input$file$datapath, sep = ";", header = TRUE) 
    
    #  I check if the file has the right columns names, otherwise the app will crash
    required_columns <- c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
    if (!all(required_columns %in% colnames(file_data))) {
      stop("Le fichier n'a pas le bon format de colonnes. Veuillez bien formater votre fichier et recommencez.")
    }
    
    return(file_data)
  })
  
  # generating the volcano plot once the data is loaded
  output$volcanoPlot <- renderPlotly({
    req(input$Pvalue_cutoff_from_input, input$log2_FoldChange_cutoff_from_input)
    generateVolcanoPlot(data(), input$Pvalue_cutoff_from_input, input$log2_FoldChange_cutoff_from_input)
  })
  
  #  creating a download button to download the volcano plot
  output$downloadVolcanoPlot <- downloadHandler(
    filename = function() {
      "volcano_plot.html"
    },
    content = function(file) {
      g <- generateVolcanoPlot(data(), input$Pvalue_cutoff_from_input, input$log2_FoldChange_cutoff_from_input)
      htmlwidgets::saveWidget(g, file, selfcontained = TRUE)
    }
  )
  
  #  I create a reactive function to filter the data according to the user's thresholds choice
  filtered_data_table <- reactive({
    data_subset <- data()
    pvalue_threshold <- input$Pvalue_cutoff_from_input
    log2foldchange_threshold <- input$log2_FoldChange_cutoff_from_input
    data_subset[data_subset$pval > pvalue_threshold | abs(data_subset$log2FC) < log2foldchange_threshold,] <- NA # I replace the points that are not under the thresholds by NA
    data_subset
  })
  
  #  creating the table
  output$table <- DT::renderDataTable({
    filtered_data_table()
  })
}
