################################################################################

### To create Ensembl URL and clickable link (HTML) on the preview table

#' @description
#' Checks if a given ID is in Ensemble ID format
#' @param gene_id : ID to check (string)
#' @return TRUE if gene_id is an Ensembl ID, FALSE otherwise.
#'
#' @example is_Ensembl_ID("ENSGXXXXXXXXXX")
#' 
is_Ensembl_ID <- function( gene_id ) {
  status <- ifelse ( str_detect( gene_id, "^ENS[:upper:]+[:digit:]{11}$"), TRUE, FALSE )
  return(status)
}

#' @description
#' Create an Ensembl link to website
#' @param organism : (string) target organism
#' @param gene_id  : (string) given Ensembl ID
#' @return (string) Raw link to the Ensembl page corresponding to the given gene_id in the given organism 
#'
#' @example create_Ensembl_link("Mus_musculus", "ENSGXXXXXXXXXX")
#' 
create_Ensembl_link <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( is_Ensembl_ID( gene_id ),
                   sprintf("https://www.ensembl.org/%s/Gene/Summary?g=%s", organism, gene_id),
                   NA )
  return(links)
}

#' @description
#' Create an Ensembl link to website as a HTML hyperlink, to include in preview data tables
#' @param organism : (string) target organism
#' @param gene_id  : (string) given Ensembl ID
#' @return (string) HTML formatted hyperlink to the Ensembl page corresponding to the given gene_id in the given organism
#'
#' @example create_Ensembl_html_link("Mus_musculus", "ENSGXXXXXXXXXX")
#' 
create_Ensembl_html_link <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( is_Ensembl_ID( gene_id ),
                   sprintf('<a href="https://www.ensembl.org/%s/Gene/Summary?g=%s" target="_blank">%s</a>', organism, gene_id, gene_id),
                   gene_id )
  return(links)
}

################################################################################

### FILTER DATATABLE

#' @description
#' Filters data table according to given adjusted p-value and Fold Change cut-offs
#' @param deseq_data  : (data table) differential expression data table
#' @param fc_cutoff   : (float) Fold Change cut-off (default : 1)
#' @param padj_cutoff : (float) adjusted p-value cut-off (default : 0.05)
#' @return (data table) differential expression data table without non-differential-expressed genes. new columns : "log2FC" and "diff_expressed" ("UP", "DOWN")
#'
#' @example filter_dt(deseq_datatable, 2, 0.01)
#' 
filter_dt <- function(deseq_data, fc_cutoff = 1, padj_cutoff = 0.05) {
  deseq_data$diff_expressed <- ifelse(
    deseq_data$log2FC <= - fc_cutoff & deseq_data$padj <= padj_cutoff, "DOWN",
    ifelse(deseq_data$log2FC >= fc_cutoff & deseq_data$padj <= padj_cutoff, "UP", "NO_DE")
  )
  return(deseq_data)
}

show_filtered_df = function(deseq_data){
  filtered_data <- deseq_data[deseq_data$diff_expressed != 'NO_DE', ]
  return(filtered_data)
}


### VOLCANO PLOT
#' @description
#' Draws a ggplot2-based volcano plot of given differential expression data table, with adjusted p-value and Fold Change cut-offs
#' @param deseq_data  : (data table) differential expression data table
#' @param title       : (string) volcano plot title (default : "Volcano plot")
#' @param xlab        : (string) label of the x axis (default : "Log2(FoldChange)")
#' @param ylab        : (string) label of the x axis (default : "-Log10(p-value adjusted)")
#' @param xlim        : (vector c(float,float) ) limits of the x axis (default : NA)
#' @param ylim        : (vector c(float,float) ) limits of the y axis (default : NA)
#' @param fc_cutoff   : (float) Fold Change cut-off (default : 1)
#' @param padj_cutoff : (float) adjusted p-value cut-off (default : 0.05)
#' @param lines       : (bool) whether lines should be drawn to represent the FC and adjuste p-value cut_offs
#' @return (ggplot2 plot) volcano plot of given differentially expression data table
#'
#' @example draw_volcano(deseq_datatable, fc_cutoff = 2, padj_cutoff = 0.01, lines = TRUE)
#' 
draw_volcano <- function(deseq_data,
                        title = "Volcano plot",
                        xlab = "Log2(FoldChange)",
                        ylab = "-Log10(p-value adjusted)",
                        xlim = NA,
                        ylim = NA,
                        fc_cutoff = 1,
                        padj_cutoff = 0.05,
                        lines = FALSE) {
  
  deseq_data$diff_expressed <- factor(deseq_data$diff_expressed, levels = c("UP", "DOWN", "NO_DE"))
  fig <- ggplot2::ggplot(data = deseq_data,
                        ggplot2::aes(
                          x = log2FC,
                          y = -log10(padj),
                          col = diff_expressed )
                        ) +
    ggplot2::scale_color_manual(values = c('DOWN' = 'blue', 'UP' = 'red', 'NO_DE' = 'grey')) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Differentially \n expressed genes \n (DE)")) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_bw()

  if(lines == TRUE){
    fig = fig +
      ggplot2::geom_hline(yintercept = -log10(padj_cutoff), col = "black") +
      ggplot2::geom_vline(xintercept = -fc_cutoff, col = "black") +
      ggplot2::geom_vline(xintercept = fc_cutoff, col = "black")
  }
  
  if(is.numeric(xlim)){
    fig = fig + ggplot2::xlim(min(xlim), max(xlim))
  }
  
  if(is.numeric(ylim)){
    fig = fig + ggplot2::ylim(min(ylim), max(ylim))
  }
  
  
  fig <- ggplotly(fig)
  
  return(fig)

}