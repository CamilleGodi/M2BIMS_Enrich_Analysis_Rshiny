################################################################################

### To create Ensembl URL and clickable link (HTML) on the preview table

# Checks the gene ID format corresponds to Ensembl
isEnsemblID <- function( gene_id ) {
  status <- ifelse ( str_detect( gene_id, "^ENS[:upper:]+[:digit:]{11}$"), TRUE, FALSE )
  return(status)
}

createEnsemblLink <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( isEnsemblID( gene_id ),
                   sprintf("https://www.ensembl.org/%s/Gene/Summary?g=%s", organism, gene_id),
                   NA )
  return(links)
}

createEnsemblHTMLlink <- function(organism, gene_id) {
  organism <- sub(" ", "_", organism)
  links <- ifelse( isEnsemblID( gene_id ),
                   sprintf('<a href="https://www.ensembl.org/%s/Gene/Summary?g=%s" target="_blank">%s</a>', organism, gene_id, gene_id),
                   gene_id )
  return(links)
}

################################################################################

draw_volcano = function(res,
                        title = "Volcanoplot",
                        xlab = "Log2(FoldChange)",
                        ylab = "-Log10(p-valeurs)",
                        xlim = NA,
                        ylim = NA,
                        fc_cutoff = 1,
                        alpha = 0.05,
                        lines = FALSE) {
  df = as.data.frame(res)
  df = df[!is.na(df$padj),]       # retirer les valeurs n'ayant pas passé le filtre indépendant si resultat de DESeq2
  df$diffexpressed = "NO_DE"
  df$diffexpressed[df$log2FC < - fc_cutoff & df$padj < alpha] = "DOWN"
  df$diffexpressed[df$log2FC > fc_cutoff & df$padj < alpha] = "UP"
  
  colors <- c()
  for (i in sort(unique(df$diffexpressed)))
  {
    if (i == "DOWN")
    { colors <- c(colors, "cyan4")
    } else if (i == "NO_DE")
    { colors <- c(colors, "gray70")
    } else if (i == "UP")
    { colors <- c(colors, "orangered1")
    }
  }
  fig = ggplot2::ggplot(data = df,
                        ggplot2::aes(
                          x = df[,"log2FC"],
                          y = -log10(df[,"pval"]),
                          col = diffexpressed
                        )) +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::geom_point() +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::guides(col = ggplot2::guide_legend(title = "Différentiellement \n exprimés")) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = colors)
  
  if(lines == TRUE){
    fig = fig +
      ggplot2::geom_hline(yintercept = -log10(alpha), col = "black") +
      ggplot2::geom_vline(xintercept = -fc_cutoff, col = "black") +
      ggplot2::geom_vline(xintercept = fc_cutoff, col = "black")
  }
  
  if(is.numeric(xlim)){
    fig = fig + ggplot2::xlim(min(xlim),max(xlim))
  }
  
  if(is.numeric(ylim)){
    fig = fig + ggplot2::ylim(min(ylim),max(ylim))
  }
  
  
  return(ggplotly(fig, tooltip = c("text"), dynamicTicks = TRUE))
}