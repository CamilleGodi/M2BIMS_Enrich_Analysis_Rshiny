# Download libraries that aren't installed, load them. Source : https://stackoverflow.com/a/44660688
load_libs <- function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if(length(need)>0) { 
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

load_libs( "Rcpp", "RcppEigen")
load_libs("BiocManager", "DT", "ggplot2", "shiny", "shinydashboard", "shinyalert", "plotly", "tidyverse")


################################################################################


# Download organisms BioCondutoR annotation databases that aren't installed yet. Adapted from source : https://stackoverflow.com/a/44660688
load_libs_biocmanager <- function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs, require, character.only = TRUE))
  need <- libs[req == FALSE]
  if(length(need)>0) { 
    BiocManager::install(need)
    lapply(need, require, character.only = TRUE)
  }
}

load_libs_biocmanager("clusterProfiler", "ReactomePA", "org.At.tair.db", "org.EcK12.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Sc.sgd.db")


################################################################################


### NOTIFICATIONS
shinyalert_wrapper <- function(title, message = "", type) {
  shinyalert(
    title = title,
    text = message,
    size = "xs", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    type = "error",
    showConfirmButton = FALSE,
    showCancelButton = TRUE,
    cancelButtonText = "Cancel",
    animation = TRUE
  )
}


################################################################################


### Organism names/codes conversion table
## Examples for Mus musculus :
# Get kegg name     : organism_conversion_table["Mus musculus", "kegg_name"]
# Get annotation db : organism_conversion_table["Mus musculus", "annotation_db"]
# Get reactome name : organism_conversion_table["Mus musculus", "reactome_name"]

add_organism_in_conversion_table <- function(conversion_table, species_name, annotation_db, kegg_name, reactome_name) {
  conversion_table <- rbind(conversion_table, c(annotation_db, kegg_name, reactome_name))
  rownames(conversion_table)[nrow(conversion_table)] <- species_name       # Updates last row's name

  return(conversion_table)
}


# Init conversion table with human
organism_conversion_table <- data.frame(annotation_db = "org.Hs.eg.db", kegg_name = "hsa", reactome_name = "human")
rownames(organism_conversion_table) <- "Homo sapiens"

# Add organisms
organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Arabidopsis thaliana", "org.At.tair.db", "ath", NA)
organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Escherichia coli (K12)", "org.EcK12.eg.db", "ecoc", NA)
organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Mus musculus", "org.Mm.eg.db", "mmu", "mouse")
organism_conversion_table <- add_organism_in_conversion_table(organism_conversion_table, "Saccharomyces cerevisiae", "org.Sc.sgd.db", "sce", "yeast")


################################################################################


### SCRAPPED : Ensembl REST API

# load_libs("httr", "jsonlite", "xml2") # for Ensembl REST API

get_ensembl_organisms_list <- function() {
  server <- "https://rest.ensembl.org"
  ext <- "/info/species?"
  
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  res <- content(r)$species
  
  ensembl_organisms <- c()
  for ( element in res ) {
    ensembl_organisms = c(ensembl_organisms, element$name)
  }
  ensembl_organisms <- sort(ensembl_organisms)
  ensembl_organisms <- sub("_", " ", ensembl_organisms)
  
  return(ensembl_organisms)
}


################################################################################