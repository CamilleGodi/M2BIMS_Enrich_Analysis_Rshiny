#' @title global.R
#' @author Benjamin Marsac
#' @description Fichier d'interface pour interface Rshiny
#' @date 2023-2024

source("./global.R")

# création du header

header = dashboardHeader(
  title = "My sub-optimal dashboard"
)



# création des items à ajouter à la sidebar
sidebar = dashboardSidebar(
  sidebarMenu(
    menuItem( icon = icon("home"),"Home", tabName = "home"),
    fileInput("select_csv","Select a csv file", accept = ".csv"),
    selectInput("dataset", "Dataset", choices = c("Homo sapiens", "Mus musculus")),
    menuItem(icon = icon("database"),"Whole Data Inspection", tabName = "whole_data_inspection"),
    menuItem(icon = icon("network-wired"),"GO term Enrichment", tabName = "go_term_enrichment"),
    menuItem(icon = icon("vial"),"Pathway Enrichment", tabName = "pathway_enrichment"),
    br(),
    menuItem(icon = icon("question"),"About", tabName = "about")
    
    
  )
)

#création du body "home"
home_body = tabItem(tabName = "home",
                    fluidRow(
                      column(12,
                             span(h2(strong("Home"), align = "center", style = 'font-size:30px;color:#020202;')))
                    ),
                    fluidRow(
                             span(h4(strong("Bienvenue dans mon application Rshiny. Cette application a vocation à faire le traitement de données issues d'analyse transcriptomique afin d'en faire l'analyse d'enrichissement."), align = "justify", style = 'font-size:15px;color:#020202;')),
                             br(),
                             span(h4(strong("Une analyse d'enrichissement permet à déterminer les voies métaboliques qui sont les plus susceptible d'être affectées à l'issue  d'une analyse différentielle. "), align = "justify", style = 'font-size:15px;color:#020202;'))
                    )
                    
)

whole_data_inspection_body = tabItem(tabName = "whole_data_inspection",
                                    #les deux boxs pour les futurs plots
                                    fluidRow(
                                      column(12,
                                             span(h2(strong("Whole data inspection"), align = "center", style = 'font-size:30px;color:#020202;')))
                                    ),
                                    fluidRow(
                                      column(width = 12,
                                             box(label = "maplot_box",
                                                 plotOutput('draw_maplot',dblclick = "maplot_double_click", brush = "maplot_brush"),
                                                 title = "MA plot",
                                                 status = "primary",
                                                 height = 470,
                                                 solidHeader = TRUE),
                                             box(label = "volcanoplot_box",
                                                 plotOutput('draw_volcanoplot',dblclick = "volcanoplot_double_click", brush = "volcanoplot_brush"),
                                                 title = "Volcanoplot",
                                                 status = "primary",
                                                 height = 470,
                                                 solidHeader = TRUE)
                                      )
                                    ),
                                    fluidRow(
                                      #les slides pour choisir les paramètres de p value et de log2foldchange, défaut sur les paramètres généraux
                                      column(5,
                                             sliderInput(label = "Set your adjusted p-value threshold",
                                                         inputId = "slider_of_pvalue_for_plots",
                                                         min = 0.001,
                                                         max = 0.2,
                                                         value = 0.05,
                                                         step = 0.001,
                                                         width = 500)),
                                      column(5,  
                                             sliderInput(label = "Set your log2(foldchange) threshold (absolute)",
                                                         inputId = "slider_of_log2_foldchange_for_volcanoplot",
                                                         min = 0,
                                                         max = 5,
                                                         value = 1,
                                                         step = 0.1,
                                                         width = 500),
                                      ),
                                      column(1,downloadButton('download',"Download the data"))
                                    ),
                                    fluidPage(
                                      box(label = "datatable_box",
                                          DT::dataTableOutput('table'),
                                          title = "Datatable",
                                          align = "center",
                                          width = 12,
                                          status = "primary",
                                          solidHeader = TRUE)
                                    )
)

# création du body
boardbody = dashboardBody(
  tabItems(home_body,
           whole_data_inspection_body
  )
)
# appel des différentes par dans une variable ui
ui = dashboardPage(
  header,sidebar,boardbody
)

