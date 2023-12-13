if(!require("shiny")){
  install.packages("shiny")
}

if(!require("shinydashboard")){
  install.packages("shinydashboard")
}

if(!require("shinyalert")){
  install.packages("shinyalert")
}

if(!require("tidyverse")){
  install.packages("tidyverse")
}
if(!require("DT")){
  install.packages("DT")
}
if(!require("tidyverse")){
  install.packages("tidyverse")
}
if(!require("plotly")){
  install.packages("plotly")
}

library(shiny)
library(shinydashboard)
library(shinyalert)
library(plotly)
library(DT)
library(tidyverse)
library(plotly)


################################################################################

### NOTIFICATIONS

shinyalertWrapper <- function(title, message = "", type) {
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