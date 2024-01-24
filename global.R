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