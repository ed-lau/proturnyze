##################### HOUSEKEEPING ###########################
library(shiny); library(shinythemes); library(networkD3); library(plotly)
#source("about_page.R"); source("news_page.R"); source("doc_page.R"); 
source("ui_elements.R"); source("ui_browse.R"); source("ui_enrich.R")
source("ui_analyze.R"); source("ui_instructions.R"); source("ui_chat.R")
source("ui_refit.R")
##################### END HOUSEKEEPING ###########################


# Define UI for dataset viewer application
shinyUI(
         navbarPage(
                theme = shinytheme("cosmo"),
                title="Proturnyze",
                
                # This is an entire tab panel as a function in ui_browse.R
                #instructions_page(),
                refit_page(),
                browse_page(),
                enrich_page(),
                chat_page()
))