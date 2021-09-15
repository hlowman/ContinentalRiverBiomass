# Shiny application for 34 site Teton run
# Heili Lowman
# September 9, 2021

#### Setup ####

# Load packages.
library(tidyverse)
library(shiny)
library(shinyWidgets)
library(here)

# All data will be pulled from this project.

# Dataset of available sites/years.
dat <- readRDS("teton_34rivers_sitesyrsgpp.rds") %>%
  select(site_name, yearf)

# Model diagnostics dataset
model_dat <- readRDS("teton_34rivers_model_diagnostics_090821.rds")

#### UI ####

# User interface:
ui <- fluidPage(
  titlePanel("Continental Metabolism Modeling"), # title
  
  # Set of 2 tabs
  tabsetPanel(type = "tabs",
              position = c("fixed=top"),
              
              # Tab 1: Covariate Figures
              tabPanel(h4("Sampling Site Covariate Figures"),
                       br(),
                       p("This part of the application allows you to sort through the available stream sites and display the original covariate data used to fit the Ricker model."),
                       br(),
                       p("You may use the drop down menu to select your site of interest, and the corresponding covariate figures should populate below, with one figure per year of data available."),
                       br(),
                       column(width = 12,
                              
                              column(width = 3,
                                     selectInput("select_site", label = h3("Select stream site:"), # site dropdown
                                     choices = unique(dat$site_name))),
                              
                              column(width = 3,
                                     htmlOutput("secondSelection"))), # year dropdown
                       hr(),
                       column(width = 12,
                       imageOutput("covplot"))),
              
              # Tab 2: Table Display of Model Output
              tabPanel(h4("Summarized Model Results & Diagnostics"),
                       br(),
                       p("This part of the application allows you to view the model results of fitting the Ricker model to the dataset."),
                       br(),
                       p("You may toggle through the column headers to sort in an ascending/descending manner, or you may used the 'Search' bar to search for a particular site."),
                       br(),
                       fluidRow(
                         column(width = 12,
                                dataTableOutput('table')
                         )
                       )
                       )
              )
)

#### Server ####

# Server:
server <- function(input, output){
  
  # Create dependent dropdown menu:
  
  output$secondSelection <- renderUI({
    
    selected_site <- input$select_site # assign chosen site to "selected_site"
    
    dat_new <- dat %>% # take original dataset
      filter(site_name %in% selected_site) %>% # filter by chosen site
      mutate(yearf_new = yearf) # new subset of years
    
    selectInput(inputId ="select_year", 
                label = h3("Select year:"), # site dropdown
                choices = unique(dat_new$yearf_new))})
  
  # Spit out the appropriate figure
  
  output$covplot <- renderImage({
    
    filename <- normalizePath(file.path("site_covariate_plots",
                                        paste(input$select_site, 
                                              input$select_year, 
                                              'covar.jpg', sep='_')))
    
    # Return a list containing the filename and alt text
    list(src = filename,
         width = "600",
         height = "600")
    
  }, deleteFile = FALSE)
  
  # Display model diagnostics
  
  output$table <- renderDataTable(model_dat)
  
}

# Combine the user interface and server:
shinyApp(ui = ui, server = server)

# Additional Shiny resources:
# https://shiny.rstudio.com/gallery/widget-gallery.html
# http://shinyapps.dreamrs.fr/shinyWidgets/
# https://shiny.rstudio.com/gallery/#user-showcase

# To deploy, use the following code:
# deployApp(here("code/teton_34sites/shiny"))
# ... so that I don't deploy the whole project.

# End of script.
