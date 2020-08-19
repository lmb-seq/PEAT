library(shinydashboard)
library(shiny)
library(DT)
library(shinyjs)

dashboardPage(
  dashboardHeader(title = "PEAT"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Select Input", tabName = "cufflinks_gene_level", icon = icon("folder-open"), startExpanded = TRUE,selected=TRUE,
               uiOutput('input_selection')
               #fluidRow(column(10, verbatimTextOutput("value"))),
               #fluidRow(column(10, verbatimTextOutput("value2")))
               #fluidRow(column(10, verbatimTextOutput("value3"))),
               #fluidRow(column(10, verbatimTextOutput("value4")))
               ),
      menuItem("Filter Input", tabName = "cufflinks_gene_level", icon = icon("filter"), #startExpanded = TRUE,
               #fluidRow(column(10, verbatimTextOutput("value"))),
               #fluidRow(column(10, verbatimTextOutput("value2"))),
               #fluidRow(column(10, verbatimTextOutput("value3"))),
               #fluidRow(column(10, verbatimTextOutput("value4"))),
               uiOutput('filter_selection')
               #fluidRow(column(10, verbatimTextOutput("value")))
      ),
      menuItem("Add annotations", tabName = "cufflinks_gene_level", icon = icon("plus"),
               uiOutput("annotation_selection")
               #fluidRow(column(10, verbatimTextOutput("value")))
      ),
      # menuItem("Analyse Selection", tabName = "cufflinks_gene_level", icon = icon("book"),
      #          uiOutput("analyse_selection")
      # ),
      menuItem("Gene lists", tabName = "cufflinks_gene_level", icon = icon("briefcase"),
               #fluidRow(column(10, verbatimTextOutput("value"))),
               #fluidRow(column(10, verbatimTextOutput("value2"))),
               #fluidRow(column(10, verbatimTextOutput("value3"))),
               #fluidRow(column(10, verbatimTextOutput("value4"))),
               uiOutput("gene_lists")
      ),
      menuItem("Table columns", tabName = "cufflinks_gene_level", icon = icon("ellipsis-v"),
               uiOutput("column_selection"),
               br()
      ),
      menuItem("Download data", tabName = "cufflinks_gene_level", icon = icon("download"),
               br(),
               #fluidRow(column(10, verbatimTextOutput("value"))),
               uiOutput("download_menu"),
               br()
      )
    ),
    
    br(),
    #Put things in here that you want to happen no matter what menu option is selected,
    #or if none are
    uiOutput("user_update"), #needed to initialise and update user_input
    br(),
    uiOutput("bottom_left_text"),
    fluidRow(column(3,offset = 1, "vb2.26")),
    br()
    
  ),
  dashboardBody(
    useShinyjs(),
    fluidRow(uiOutput("tab_panel"))
  )
)



