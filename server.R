##################### HOUSEKEEPING ###########################
library(shiny); library(dplyr); library(ggplot2); library(reshape2); library(heatmap3)

#library(biomaRt); (Migrating to using AWS RDS for all data)

library(R.utils); library(networkD3); library(tidyr); library(caret); library(DT); require(stringr)
#library(mclust); library(ipred); library(plotly); 
library(RMySQL)

source("server_elements.R", local=T)


options(shiny.maxRequestSize=10*1024^2) 
################################################################



shinyServer(function(input, output, session) {
        
        
        ###
        ### Housekeeping
       
        # Loading a local annotation file
        annotations <- reactive({ readRDS("data/annot_10090.Rds") 
               
        })
        
        # Server-side functions for fitting models
        source("server_models.R", local=T)
        
        # Server-side functions for chatting functions
        source("server_interactive.R", local=T)
        
        # Server-side functions for entire-file refitting functions
        #source("server_refit.R", local=T)
        
        # Server-side functions for enrichment functions
        #source("server_enrich.R", local=T)

        # Server-side functions for chatting functions
        source("server_chat.R", local=T)
})