# Download libraries that aren't installed, load them. Source : https://stackoverflow.com/a/44660688
load_libs <-function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs,require, character.only = TRUE))
  need <- libs[req == FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need, require, character.only = TRUE)
  }
}

load_libs("shiny", "shinydashboard", "shinyalert", "plotly", "DT", "tidyverse")
load_libs("httr", "jsonlite", "xml2") # for Ensembl REST API

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


### Ensembl REST API
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