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
