################################################################################

### To create Ensembl URL and clickable link (HTML) on the preview table

# Checks the gene ID format corresponds to Ensembl
is_Ensembl_ID <- function( gene_id ) {
  status <- ifelse ( str_detect( gene_id, "^ENS[:upper:]+[:digit:]{11}$"), TRUE, FALSE )
  return(status)
}

create_Ensembl_link <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( is_Ensembl_ID( gene_id ),
                   sprintf("https://www.ensembl.org/%s/Gene/Summary?g=%s", organism, gene_id),
                   NA )
  return(links)
}

create_Ensembl_html_link <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( is_Ensembl_ID( gene_id ),
                   sprintf('<a href="https://www.ensembl.org/%s/Gene/Summary?g=%s" target="_blank">%s</a>', organism, gene_id, gene_id),
                   gene_id )
  return(links)
}

################################################################################

### FILTER DATATABLE


filter_dt <- function(deseq_data, fc_cutoff = 1, padj_cutoff = 0.05) {
  deseq_data$diff_expressed <- ifelse(
    deseq_data$log2FC <= - fc_cutoff & deseq_data$padj <= padj_cutoff, "DOWN",
    ifelse(deseq_data$log2FC >= fc_cutoff & deseq_data$padj <= padj_cutoff, "UP", "NO_DE")
  )
  filtered_data <- deseq_data[deseq_data$diff_expressed != 'NO_DE', ]
  return(filtered_data)
}




### VOLCANO PLOT

draw_volcano <- function(deseq_data,
                        title = "Volcanoplot",
                        xlab = "Log2(FoldChange)",
                        ylab = "-Log10(p-value adjusted)",
                        xlim = NA,
                        ylim = NA,
                        fc_cutoff = 1,
                        padj_cutoff = 0.05,
                        lines = FALSE) {
  deseq_data <- deseq_data[!is.na(deseq_data$padj), ]       # retirer les valeurs n'ayant pas passé le filtre indépendant si resultat de DESeq2
  deseq_data$diff_expressed <- ifelse(
    deseq_data$log2FC <= -fc_cutoff & deseq_data$padj <= padj_cutoff, "DOWN",
    ifelse(deseq_data$log2FC >= fc_cutoff & deseq_data$padj <= padj_cutoff, "UP", "NO_DE")
  )
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