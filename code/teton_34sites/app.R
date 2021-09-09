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
dat <- readRDS(here("data_working/teton_34rivers_sitesyrsgpp.rds")) %>%
  select(site_name, yearf)

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
                       p("This application allows you to sort through the available stream sites and display the original covariate data used to fit the Ricker model. You may use the drop down menu to select your site of interest, and the corresponding covariate figures should populate below, with one figure per year of data available."),
                       br(),
                       column(width = 12,
                              
                              column(width = 3,
                                     selectInput("select_site", label = h3("Select stream site:"), # site dropdown
                                     choices = unique(dat$site_name))),
                              
                              column(width = 3,
                                     htmlOutput("secondSelection"))), # year dropdown
                       hr(),
                       plotOutput("covplot")),
              
              # Tab 2: Table Display of Model Output
              tabPanel(h4("Summarized Model Results & Diagnostics")))
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
    
    filename <- normalizePath(here("figures/teton_34sites/site_covariate_plots",
                                        paste(input$select_site, 
                                              input$select_year, 
                                              'covar.jpg', sep='_')))
    
    # Return a list containing the filename and alt text
    list(src = filename)
    
  }, deleteFile = FALSE)
  
}

# Combine the user interface and server:
shinyApp(ui = ui, server = server)

# Additional Shiny resources:
# https://shiny.rstudio.com/gallery/widget-gallery.html
# http://shinyapps.dreamrs.fr/shinyWidgets/
# https://shiny.rstudio.com/gallery/#user-showcase

# End of script.