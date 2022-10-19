# Shiny application for 182 site Teton run
# Heili Lowman
# Created: September 9, 2021
# Revised for the second Teton run: February 5, 2022
# Revised again for third Teton run: October 18, 2022

#### Setup ####

# Load packages.
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(here)

# All data will be pulled from this project.

# Dataset of available sites/years.
dat <- readRDS("df_182rivers_sitesyrsgpp.rds") 

sy_dat <- dat %>%
  select(site_name, year)

# Dataset of site names.
names_dat <- dat %>%
  select(site_name, long_name) %>%
  distinct() %>% # removes duplicates
  rename("NWIS ID" = "site_name",
         "Long Site Name" = "long_name")

# Model diagnostics dataset
model_dat <- readRDS("teton_182rivers_model_diags_101522.rds")

#### UI ####

# User interface:
ui <- fluidPage(
  titlePanel("River Resilience Modeling"), # title
  
  # Set of 3 tabs
  tabsetPanel(type = "tabs",
              
              ##### Tab 1: Covariate Figures #####
              tabPanel(h4("Covariate Figures"),
                       br(),
                       p("This part of the application allows you to sort through the available stream sites (n = 182) and display the original covariate data used to fit the model."),
                       br(),
                       p("You may use the drop down menu to select your site of interest, and the corresponding covariate figures should populate below."),
                       br(),
                       column(width = 12,
                              
                              column(width = 3,
                                     selectInput("select_site", label = h3("Select stream site:"), # site dropdown
                                     choices = unique(sy_dat$site_name)))),
                              
                              hr(),
                       column(width = 12,
                       imageOutput("covplot"))),
              
              ##### Tab 2: Site Listing #####
              tabPanel(h4("Site Names and NWIS IDs"),
                       br(),
                       p("This part of the application provides a table to convert between NWIS identification numbers and long names of stream sites."),
                       fluidRow(
                         column(width=12,
                                dataTableOutput('names'))
                       )),
              
              ##### Tab 3: Table Display of Model Output #####
              tabPanel(h4("Summarized Model Results & Diagnostics"),
                       br(),
                       p("This part of the application allows you to view the model results of fitting the model to the dataset."),
                       br(),
                       p("You may toggle through the column headers to sort in an ascending/descending manner, or you may used the 'Search' bar to search for a particular site."),
                       br(),
                       fluidRow(
                         column(width = 12,
                                dataTableOutput('table')
                         )
                       )
                       ),
              
              ##### Tab 4: S vs. C Figures #####
              tabPanel(h4("S vs. C Figures for All Iterations"),
                       br(),
                       p("This part of the application allows you to sort through the available stream sites and display the sensitivy of the persitence curve (s) and critical disturbance threshold (c) values for every iteration of the model run."),
                       br(),
                       p("You may use the drop down menu to select your site of interest, and the corresponding figures should populate below."),
                       br(),
                       column(width = 12,
                              
                              column(width = 3,
                                     selectInput("select_site2", label = h3("Select stream site:"), # site dropdown
                                                 choices = unique(sy_dat$site_name)))),
                       hr(),
                       column(width = 12,
                              imageOutput("scplot")))
              
              )
)

#### Server ####

# Server:
server <- function(input, output){
  
  ##### Tab 1 #####
  
  # Spit out the appropriate figure - tab 1
  
  output$covplot <- renderImage({
    
    filename <- normalizePath(file.path("site_covariate_plots",
                                        paste(input$select_site, 
                                              'covar.jpg',
                                              sep = "")))
    
    # Return a list containing the filename and alt text
    list(src = filename,
         width = "600",
         height = "600")
    
  }, deleteFile = FALSE)
  
  ##### Tab 2 #####
  
  # Display site names
  
  output$names <- renderDataTable(names_dat)
  
  ##### Tab 3 #####
  
  # Display model diagnostics
  
  output$table <- renderDataTable(model_dat)
  
  ##### Tab 4 #####
  
  # Spit out the appropriate figure - tab 4
  
  output$scplot <- renderImage({
    
    filename <- normalizePath(file.path("site_sc_plots",
                                        paste(input$select_site2, 
                                              'sc.jpg',
                                              sep = "")))
    
    # Return a list containing the filename and alt text
    list(src = filename,
         width = "600",
         height = "600")
    
  }, deleteFile = FALSE)
  
}

# Combine the user interface and server:
shinyApp(ui = ui, server = server)

# Additional Shiny resources:
# https://shiny.rstudio.com/gallery/widget-gallery.html
# http://shinyapps.dreamrs.fr/shinyWidgets/
# https://shiny.rstudio.com/gallery/#user-showcase

#### Deployment ####

# To deploy, use the following code:
# deployApp(here("code/shiny"))
# ... so that I don't deploy the whole project.

# End of script.
