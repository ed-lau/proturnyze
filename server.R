##################### HOUSEKEEPING ###########################
library(shiny); library(dplyr); library(ggplot2); library(reshape2); library(heatmap3)
library(biomaRt); library(R.utils); library(networkD3); library(tidyr); library(caret)
library(mclust); library(ipred); library(plotly); library(DT)
source("server_elements.R", local=T)

options(shiny.maxRequestSize=10*1024^2) 
################################################################



shinyServer(function(input, output, session) {
        
        source("server_chat.R", local=T)
        
        ######
        ######
        ###### GO and Data Retrieval Functions
        ######
        ######
        
        withProgress(message = 'Updating GO database...', value = 0.5, {
                last_update <- read.table("data/lastupdate.txt", as.is=T)
                incProgress(0.1)
                # If last update was a week ago, do:
                if (Sys.Date() - as.Date(last_update[1,1]) >= 20){
                        #all_human_go <- read.table("http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=-1&tax=9606&db=UniProtKB&col=proteinID,goID,goName,aspect", fill=T, header=T, as.is=T, sep="\t")
                        #all_human_go <- all_human_go[,1:4]
                        #names(all_human_go) <- c("uniprot","go_id","go_name","aspect")
                        #all_human_go <- all_human_go %>% distinct(uniprot, go_id)
                        #write.table(all_human_go, "data/all_human_go.txt", row.names=F, quote=F, sep="\t")
        
                        all_mouse_go <- read.table("http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=-1&tax=10090&db=UniProtKB&col=proteinID,goID,goName,aspect", fill=T, header=T, as.is=T, sep="\t")
                        all_mouse_go <- all_mouse_go[,1:4]
                        names(all_mouse_go) <- c("uniprot","go_id","go_name","aspect")
                        all_mouse_go <- all_mouse_go %>% distinct(uniprot, go_id)
                        write.table(all_mouse_go, "data/all_mouse_go.txt", row.names=F, quote=F, sep="\t")
                        incProgress(0.4)
                        
                        write.table(Sys.Date(),"data/lastupdate.txt", row.names=F, quote=F, col.names=F, sep="\t")
                        
                } else {incProgress(0.4)}
                incProgress(0.4)
        })
        
        
#         getAllHumanGO <- reactive({
#                 # Loading a file on all Gene Ontology terms associated with all human (add mouse later) genes
#                 read.table("data/all_human_go.txt", sep ="\t", header=T, as.is=T)
#         })
        
        getAllMouseGO <- reactive({
                # Loading a file on all Gene Ontology terms associated with all human (add mouse later) genes
                read.table("data/all_mouse_go.txt", sep ="\t", fill=T, header=T, as.is=T, quote="")
        })
        
        output$GOUpdateDate <- renderText({
                readLines("data/lastupdate.txt")
        })        
        
        
        
  
  # Return the requested dataset
  datasetInput = reactive({

    input$proteinValue 
    proteinValue <- filter(proteinValue, Strain %in% input$strain)
    proteinValue <- filter(proteinValue, Factor %in% input$factor)
    proteinValue <- filter(proteinValue, go_name %in% input$compartment)
    #proteinValue <- filter(proteinValue, proteinValue %in% input$everthing)
    proteinValue
    })
  

    updateSelectizeInput(session, 'strain', choices = levels(proteinValue$Strain), selected=levels(proteinValue$Strain), server=T)
    updateSelectizeInput(session, 'factor', choices = levels(proteinValue$Factor), selected=levels(proteinValue$Factor), server=T)
    updateSelectizeInput(session, 'compartment', choices = levels(proteinValue$go_name), selected=levels(proteinValue$go_name), server=T)
   # updateSelectizeInput(session, 'everything', choices= proteinValue, selected=NULL, server = TRUE)

  # based on variable values
  
  # Generate a summary of the dataset
  output$summary <- renderPrint({
          datasetInput() %>% dplyr::select(-PN, -go_name, -Uniprot) %>% summary()
  })
  
  # Render Data Table
  output$view <- renderDataTable({
    datasetInput() %>% distinct(Uniprot, Strain, Factor) %>% dplyr::select(-go_name) %>% 
                  mutate(Uniprot = makeUniLink(Uniprot), k=round(k, 3), t=round(t, 3), p=round(p, 3))
  }, escape=F, options = list(lengthMenu = c(10, 50, 100), pageLength = 10))
  
  # Make Uniprot Link
  makeUniLink <- function(val) {
          paste('<a href="http://www.uniprot.org/uniprot/', val,'" target="_blank" class="btn btn-primary"> ',val,'</a>', sep = "")
  }
  
  output$viewExplorePlot <- renderPlotly({
  
          plotData <- datasetInput() %>% dplyr::select(-GN, -PN)
          
          if (input$pickExplorePlotType == "scat"){        
                  g <- ggplot(data = plotData, aes(y = k, x = p))
                  g <- g + facet_wrap(~Factor)
                  g <- g + geom_density2d(alpha = 1)
                  g <- g + theme_bw(base_size = 10) #+ coord_cartesian(ylim = c(-3,3))
                  g <- g + theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, hjust = 0)) + labs(x="log2 Abundance", y = "log2 k")
                  g
          
          } else if (input$pickExplorePlotType == "box"){
                  g <- ggplot(data = plotData, aes(y = k, x = Strain))
                  g <- g + facet_wrap(~Factor)
                  g <- g + geom_boxplot(alpha = 0.25)
                  g <- g + theme_bw(base_size = 10) #+ coord_cartesian(ylim = c(-3,3))
                  g <- g + theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, hjust = 0)) + labs(x="Mouse Strain", y = "log2 k")
                  g
          } else if (input$pickExplorePlotType == "cc"){
                  g <- ggplot(data = plotData, aes(y = k, x = go_name))
                  g <- g + facet_wrap(~Factor)
                  g <- g + geom_violin(alpha = 0.25)
                  g <- g + theme_bw(base_size = 10) #+ coord_cartesian(ylim = c(-3,3))
                  g <- g + theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, hjust = 0)) + labs(x="Compartment", y = "log2 k")
                  g
          }
          ggplotly(g)
  })
  
  ######
  ######
  ###### Data Analysis Demo Functions
  ######
  ######
  
  filterDemoData <- function(){
          
          # Dependency on the Run Button
          #input$startClusterButton
          
          # Get the checkbox inputs
          datatypeChoice <- input$datatypeChoice
          strainChoice <- input$strainChoice
          factorChoice <- input$factorChoice
          
          datatypeColumns <- rep(F,ncol(demoData))
          if(1 %in% datatypeChoice){ datatypeColumns[14:25] <- T}
          if(2 %in% datatypeChoice){ datatypeColumns[2:13] <- T}
          if(3 %in% datatypeChoice){ datatypeColumns[26:37] <- T}
          
          strainColumns <- rep(F,ncol(demoData))
          if(1 %in% strainChoice){ strainColumns[c(2,8,14,20,26,32)] <- T}
          if(2 %in% strainChoice){ strainColumns[c(3,9,15,21,27,33)] <- T}
          if(3 %in% strainChoice){ strainColumns[c(4,10,16,22,28,34)] <- T}
          if(4 %in% strainChoice){ strainColumns[c(5,11,17,23,29,35)] <- T}
          if(5 %in% strainChoice){ strainColumns[c(6,12,18,24,30,36)] <- T}
          if(6 %in% strainChoice){ strainColumns[c(7,13,19,25,31,37)] <- T}

          factorColumns <- rep(F,ncol(demoData))
          if(1 %in% factorChoice){ factorColumns[c(2:7,14:19,26:31)] <- T}
          if(2 %in% factorChoice){ factorColumns[c(8:13,20:25,32:37)] <- T}
          
          
          allColumns <- datatypeColumns & strainColumns & factorColumns
          allColumns[1] <- T
          return(allColumns)
          
  }
  
  output$filterDemoDataFeedback <- renderText({
          if(sum(filterDemoData()) >= 4){
                  return(
                          paste("You have selected ",
                                as.character(sum(filterDemoData())-1),
                                " parameters for analysis.",
                                #"The following columns will be analyzed:",
                                #as.character(colnames(demoData)[filterDemoData()]),
                                " Please proceed to the next step."
                                )
                        ) 
          } else return("Please select at least 3 columns for data analysis.")
  })
  
  subsetDemoDataByFilter <- reactive({
          # Load the Demo Data again, and filter by columns selected
          datafile <- demoData[filterDemoData()]
          row.names(datafile) <- datafile$uniprot
          datafile <- datafile[-1]
          # Filter only by complete cases - if no imputation selected.
          datafile <- datafile[complete.cases(datafile),]
          datafile
  })
  
  clusterDemoData <- reactive({
          #if (is.null(imputeTable())){return(NULL)}
          
          # Dependency on the Run Button
          input$startClusterButton
          datafile <- isolate(subsetDemoDataByFilter())
          
          withProgress(message = "Clustering selected data...", value = 0.9, {
                  
                  # Shrink data to 4 principal components if more than 4 parameters selected - to save clustering time.
                  if(ncol(datafile) > 4){
                          datafile2 <- preProcess(datafile, method = c("pca"), pcaComp = 4)
                          datafile3 <- predict(datafile2, datafile)
                          row.names(datafile3) <- row.names(datafile)
                  } else {datafile3 <- datafile}
                  
                  clus <- Mclust(datafile3)
          })
          return(clus)
  })
  
  output$renderDemoClusters <- renderPlot({
          # Dependency on the Run Button
          input$startClusterButton
          # Read the Cluster results
          clus <- isolate(clusterDemoData())
          # Plot the MCluster results
          plot(clus, what = "classification")
          })
  
  numberOfClusters <- reactive({
          clus <- isolate(clusterDemoData())
          clus_tree <- as.data.frame(clus$classification) 
          colnames(clus_tree) <- "cluster"
          clus_tree <- clus_tree %>% distinct(cluster)
          clus_tree$cluster %>% max()
                })
  
  # Update the UI input for picking the cluster to display
  observe({
          # Create a list of n items where n is the number of clusters from clustering
          s_options <- list()
          for(x in 1:numberOfClusters()){
          s_options[[sprintf("Cluster %d", x)]] <- x
          }
          updateSelectInput(session, "pickCluster", label="Choose cluster to display", choices = s_options, selected = 1)
  })
  
  output$renderDemoSingleCluster <- renderPlotly({
          
          # Dependency on the Run Button
          input$startClusterButton
          
          clus <- isolate(clusterDemoData())
          clus_tree <- as.data.frame(clus$classification) 
          colnames(clus_tree) <- "cluster"
          clus_tree$uniprot <- rownames(clus_tree)
          
          # Load the Demo Data again, and filter by columns selected
          datafile <- isolate(subsetDemoDataByFilter())
          
          datafile_plot <- mutate(datafile, uniprot = rownames(datafile)) %>% melt()
          datafile_plot <- datafile_plot %>% left_join(clus_tree) %>% filter(cluster == input$pickCluster) #)
          #print(glimpse(datafile_plot))
          
          box_plot <- ggplot(data = datafile_plot, aes(x = variable, y = value))
          box_plot <- box_plot + geom_boxplot(alpha = 0.25)
          box_plot <- box_plot + theme_bw(base_size = 10) #+ coord_cartesian(ylim = c(-3,3))
          box_plot <- box_plot + theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, hjust = 0)) + labs(x="Sample", y = "Value")
          box_plot <- ggplotly(box_plot)
          
          line_plot <- ggplot(data = datafile_plot, aes(x = variable, y = value, group = uniprot))
          line_plot <- line_plot + geom_line(alpha = 0.25)
          line_plot <- line_plot + theme_bw(base_size = 10) #+ coord_cartesian(ylim = c(-3,3))
          line_plot <- line_plot + theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 90, hjust = 0)) + labs(x="Sample", y = "Value")
          line_plot <- ggplotly(line_plot)
          
          if(input$pickPlotType == "boxplot"){
                  box_plot      # Note to self; plotly doesn't work with the explicit return() function
          }else{ 
                  line_plot     # Note to self; plotly doesn't work with the explicit return() function
          }
  })
  
  
  prepClusterDownload <- reactive({
          
          clus <- clusterDemoData()
          clus_tree <- as.data.frame(clus$classification) 
          colnames(clus_tree) <- "cluster"
          clus_tree$uniprot <- rownames(clus_tree)
          # Filter by the selected cluster
          clus_tree <- clus_tree %>% filter(cluster == input$pickCluster)
          
          # Load the filtered Demo Data
          datafile <- isolate(subsetDemoDataByFilter())
          row.names(datafile) %>% unlist() %>% as.data.frame() %>% filter(. %in% clus_tree$uniprot)
  })
    
  prepBackgroundDownload <- reactive({
          datafile <- isolate(subsetDemoDataByFilter())
          row.names(datafile) %>% unlist() %>% as.data.frame()
  })      
          
  output$downloadCluster <- downloadHandler(
          filename = function(){paste("cluster", input$pickCluster, ".txt", sep="")},
          content = function(file){
                  write.table(prepClusterDownload(), file, quote=F, row.names=F, col.names=F)
  })
  
  output$downloadBackground <- downloadHandler(
          filename = function(){paste("background", ".txt", sep="")},
          content = function(file){
                  write.table(prepBackgroundDownload(), file, quote=F, row.names=F, col.names=F)
          })
  
  
  
  ######
  ######
  ###### GeneList Analysis Functions
  ######
  ######

          # Associate the gene list and background with Gene Ontologies
  getGOList <- reactive({
          withProgress(message = 'Retrieving GO Annotations...', value = 0, {
                  
                  all_go <- getAllMouseGO()
                  
                  genelist <- processLists() %>% filter(type == "genelist") %>% dplyr::select(uniprot)
                  background <- processLists() %>% filter(type == "background") %>% dplyr::select(uniprot)
                  
                  #### ALT: Get the GO terms from biomaRt instead (this could take up to 3 minutes)
                  #ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host = "www.ensembl.org")
                  #all_go <- getBM(attributes=c("uniprot_swissprot","go_id","name_1006","namespace_1003"), filters = "uniprot_swissprot", values=background, mart= ensembl)
                  #all_go <- mutate(all_go, uniprot = uniprot_swissprot)
                  
                  genelist_go <- all_go %>% filter(uniprot %in% unlist(genelist)) %>% distinct(uniprot,go_id) %>% mutate(type = "genelist")
                  background_go <- all_go %>% filter(uniprot %in% unlist(background)) %>% distinct(uniprot,go_id) %>% mutate(type = "background")
                  
          })
          bind_rows(genelist_go, background_go)
  })
  
          # Get all reviewed mouse UniProt ID from UniProt
  getMouseUniprot <- reactive({
          withProgress(message = "Fetching Uniprot entries...", value = 0.9, {
                  allUniprot <- read.table("http://www.uniprot.org/uniprot/?query=organism:10090+reviewed:yes&format=tab&columns=id,database(OMA)", header=T, fill=T, as.is=T) %>% dplyr::select(uniprot = Entry, OMA = Cross.reference)
          })
          allUniprot
  })
  
          
          # Handle the user input values in the Gene List text box. Map them to the complete human UniProt list, get homology if necessary.
  userInputGeneList <- reactive({
          input$goButton
          genelist <- isolate(input$genelist)
          genelist <- strsplit(genelist,"[[:space:]+[:punct:]+]") %>% as.data.frame() %>% distinct() 
          names(genelist) <- "uniprot"
          genelist <- genelist %>% filter(uniprot != "") %>% filter(nchar(as.character(uniprot)) == 6)
          genelist
  })
          
          # Handle the user input values in the Background List text box. Map them to the complete human UniProt list, get homology if necessary
  userInputBackground <- reactive({
          input$goButton
          background <- isolate(input$background)
          background <- strsplit(background,"[[:space:]+[:punct:]+]") %>% as.data.frame() %>% distinct()
          names(background) <- "uniprot"
          background <- background %>% filter(uniprot != "") %>% filter(nchar(as.character(uniprot)) == 6)
          background
  })  
          
          # Get the species taxonomy
  userInputTaxonomy <- reactive({
          genelist <- userInputGeneList()
          gene <- genelist$uniprot[1]
          getTax(gene)
  })
          
  getTaxonomyUniprot <- reactive({
          withProgress(message = 'Fetching Orthologs...', value = 0.5, {
                  tax <- userInputTaxonomy()
                  taxUniprot <- getTaxUniprot(tax)
                  incProgress(0.5)
                  #genelist2 <- taxUniprot2 %>% dplyr::select(uniprot)
          })
          taxUniprot
  })
  
          
  # Clean up the lists and match them to UniProt mouse
  processLists <- reactive ({
                  
          mouseUniprot <- getMouseUniprot()  
          genelist <- userInputGeneList() 
          background <- userInputBackground()
          
          tax <- userInputTaxonomy()
          
          # If the input list is human, then simply filter the list of Uniprot IDs with the downloaded human files.
          if (tax == "10090") {
                  genelist <- genelist %>% filter(uniprot %in% mouseUniprot$uniprot)
                  if (nrow(background) > 0){
                          background <- background %>% filter(uniprot %in% mouseUniprot$uniprot) %>% bind_rows(genelist) %>% distinct()
                  } else {background <- mouseUniprot %>% dplyr::select(uniprot)}
                  
                  # If the input list is not human, get all Uniprot IDs from the associated taxonomy, then filter.
          } else if (tax != "10090"){
                  genelist <- getTaxonomyUniprot() %>% filter(tax_uniprot %in% genelist$uniprot) %>% left_join(mouseUniprot) %>% filter(uniprot != "NA") %>% dplyr::select(uniprot)
                  if (nrow(background) > 0){
                          background <- getTaxonomyUniprot() %>% filter(tax_uniprot %in% background$uniprot) %>% left_join(mouseUniprot) %>% filter(uniprot != "NA") %>% dplyr::select(uniprot) %>% bind_rows(genelist) %>% distinct()
                  } else {background <- getTaxonomyUniprot() %>% left_join(mouseUniprot) %>% filter(uniprot != "NA") %>% dplyr::select(uniprot)}
          }
          
          # Output a combined list
          genelist$type <- "genelist"
          background$type <- "background"
          bind_rows(genelist, background)
          })
          
  output$genelist <- renderDataTable({
                  processLists() 
          }, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
          
          makeUniLink <- function(val) {
                  paste('<a href="http://www.uniprot.org/uniprot/', val,'" target="_blank" class="btn btn-primary"> ',val,'</a>', sep = "")
          }
          
  output$downloadData <- downloadHandler(
                  filename = function() { paste(userInputTopic(), '.csv', sep='') },
                  content = function(file) {
                          write.csv(head(results(),userInputDlCount()), file)})
          
  output$nGenelist <- renderText({
                  tax_text <- read.table(paste("http://www.uniprot.org/uniprot/?query=organism:",userInputTaxonomy(),"&format=tab&limit=1&columns=organism", sep=""), header=T, fill=T, as.is=T, sep="\t")
                  tax_text <- paste(as.character(tax_text[1,1])," (",userInputTaxonomy(),")",sep="")
                  nGenelist <- processLists() %>% filter(type=="genelist") %>% nrow()
                  nBackground <- processLists() %>% filter(type=="background") %>% nrow()
                  if(nGenelist > 0){
                          paste("You entered ", nGenelist," recognized UniProt entries in the gene list
                                and ", nBackground, " recognized UniProt entries in the background. The assumed
                                taxonomy for the Uniprot entries is ", tax_text, ". Proceed to the next step.", sep="")
                  }
                  })
          
          
  enrichTerms <- reactive({
                  genelist_go <- getGOList() %>% filter(type == "genelist") %>% dplyr::select(uniprot, go_id, go_name, aspect)
                  background_go <- getGOList() %>% filter(type == "background") %>% dplyr::select(uniprot, go_id, go_name, aspect)
                  
                  output <- enrichFunction(genelist_go, background_go, pfilter = 0.05)
                  output
          })
          
  makeGOMatrix <- reactive({
                  enriched <- enrichTerms()$go_id
                  sig_GO_terms <- getGOList() %>% filter(type == "genelist") %>% filter(go_id %in% enriched) %>% dplyr::select(go_id) %>% distinct()
                  
                  GO_matrix <- as.data.frame(matrix(data = NA, nrow = nrow(sig_GO_terms), ncol= nrow(sig_GO_terms)))
                  rownames(GO_matrix) <- sig_GO_terms$go_id
                  colnames(GO_matrix) <- sig_GO_terms$go_id
                  
                  for(i in 1: nrow(GO_matrix)){
                          for(j in 1:ncol(GO_matrix)){
                                  GO_term_1_uniprot <- getGOList() %>% filter(go_id == rownames(GO_matrix)[i])
                                  GO_term_2_uniprot <- getGOList() %>% filter(go_id == colnames(GO_matrix)[j])
                                  overlap <- GO_term_1_uniprot %>% filter(uniprot %in% GO_term_2_uniprot$uniprot) %>% nrow()
                                  GO_matrix[i,j] <- overlap
                          }
                  }
                  
                  GO_matrix_edges <- GO_matrix %>% mutate(key = rownames(GO_matrix)) %>% tidyr::gather()
                  colnames(GO_matrix_edges) <- c("GO_term_1","GO_term_2","overlap")
                  GO_matrix_edges <- filter(GO_matrix_edges, overlap > 0)
                  GO_matrix_edges
          })
          
  output$go_network <- renderSimpleNetwork({
                  #Node <- c(as.character(makeGOMatrix()$GO_term_1), as.character(makeGOMatrix()$GO_term_2)) %>% as.data.frame() %>% distinct()
                  #Node$Group <- "Default"
                  simpleNetwork(Data = makeGOMatrix())
                  #forceNetwork(Links = GO_matrix_edges, Source = 1, Target = 2, Nodes = GO_matrix_edges, Node)
          })
          
  makeGoLink <- function(val) {
                  paste('<a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=', val,'" target="_blank" class="btn btn-primary"> ',val,'</a>', sep = "")
          }
          
          output$enrich_out <- renderDataTable({
                  enrichTerms() %>% filter(aspect == input$go_type) %>%
                          mutate(est = round(est, 1), conf1 = round(conf1,1), p = signif(p,2), BH = signif(BH,2), conf2 = round(conf2, 1)) %>%
                          mutate(go_id = makeGoLink(go_id)) %>%
                          dplyr::select(go_id, go_name, aspect, proteins, F = est, Lo = conf1, Hi = conf2, P = p, BH) %>% dplyr::arrange(BH)
          },escape = FALSE)
          
  source("server_refit.R", local=T)
})