##################### HOUSEKEEPING ###########################
library(shiny); library(dplyr); library(ggplot2); library(reshape2); library(heatmap3)

#library(biomaRt); (Migrating to using AWS RDS for all data)

library(R.utils); library(networkD3); library(tidyr); library(caret)
library(mclust); library(ipred); library(plotly); library(DT); require(stringr)
source("server_elements.R", local=T)


options(shiny.maxRequestSize=10*1024^2) 
################################################################



shinyServer(function(input, output, session) {
        
        
        
        ###
        ### Update databases (Migrating to AWS)
        ###
        # withProgress(message = 'Updating GO database...', value = 0.5, {
        #         last_update <- read.table("data/lastupdate.txt", as.is=T)
        #         incProgress(0.1)
        #         # If last update was a week ago, do:
        #         if (Sys.Date() - as.Date(last_update[1,1]) >= 20){
        #                 #all_human_go <- read.table("http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=-1&tax=9606&db=UniProtKB&col=proteinID,goID,goName,aspect", fill=T, header=T, as.is=T, sep="\t")
        #                 #all_human_go <- all_human_go[,1:4]
        #                 #names(all_human_go) <- c("uniprot","go_id","go_name","aspect")
        #                 #all_human_go <- all_human_go %>% distinct(uniprot, go_id)
        #                 #write.table(all_human_go, "data/all_human_go.txt", row.names=F, quote=F, sep="\t")
        #                 
        #                 all_mouse_go <- read.table("http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=-1&tax=10090&db=UniProtKB&col=proteinID,goID,goName,aspect", fill=T, header=T, as.is=T, sep="\t")
        #                 all_mouse_go <- all_mouse_go[,1:4]
        #                 names(all_mouse_go) <- c("uniprot","go_id","go_name","aspect")
        #                 all_mouse_go <- all_mouse_go %>% distinct(uniprot, go_id)
        #                 write.table(all_mouse_go, "data/all_mouse_go.txt", row.names=F, quote=F, sep="\t")
        #                 incProgress(0.4)
        #                 
        #                 write.table(Sys.Date(),"data/lastupdate.txt", row.names=F, quote=F, col.names=F, sep="\t")
        #                 
        #         } else {incProgress(0.4)}
        #         incProgress(0.4)
        # })
        
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