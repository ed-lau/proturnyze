###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###


enrich_page <- function(){
        
        tabPanel("Enrich",
                 tagline(),
                 sidebarLayout(
                         
                        sidebarPanel(
                                         h3("Functional Enrichment"),
                                         tags$hr(),
                                         p("This page contains a helper function that enriches for Gene Ontology terms in
                                           any protein list against a particular background (or the entire Uniprot, if left blank)."),
                                         br(),
                                         p("This function is currently optimized to handle mouse proteins for demo purposes. If proteins from any other
                                           taxonomies are put in, the program will attempt to map the Uniprot ID to the mouse homolog.")
                                         )
                         , 
                         mainPanel(
                                 tabsetPanel(type ="tabs",
                                             
                                             tabPanel("Step 1: Start",
                                                      br(),
                                                      h3("Step 1A: Enter Uniprot ID for the foreground gene list:"),
                                                      p("Example: ",code("A2AAJ9, A2ADY9, A2AGT5, A2AIL4, A2AMM0, A2AN08, A2APY7, A2AQ25")),
                                                      textInput("genelist", label=NULL, width="100%", value = "A2AAJ9, A2ADY9, A2AGT5, A2AIL4, A2AMM0, A2AN08, A2APY7, A2AQ25"),
                                                      tags$hr(),
                                                      h3("Step 1B: Enter Uniprot ID for the background gene list:"),
                                                      p("Leave as blank to perform enrichment analysis against the whole proteome."),
                                                      textInput("background", label=NULL, width="100%"),
                                                      tags$hr(),
                                                      h3("Step 1C: Select Gene Ontology aspect to use:"),
                                                      radioButtons("go_type", "Gene Ontology Aspect:",
                                                                   c("Biological Process" = "Process",
                                                                   "Cellular Component" = "Component",
                                                                   "Molecular Function" = "Function")),
                                                      tags$hr(),
                                                      h3("Step 1D: Press Go:"),
                                                      actionButton("goButton","Upload Gene Lists", width="100%", icon("play")),
                                                      h6("Note: This may take up to a minute.", align="center"),
                                                      br(),
                                                      h5(textOutput("nGenelist")),
                                                      br(),
                                                      dataTableOutput("genelist"),
                                                      br()
                                             ),
                                             tabPanel("Step 2: Output",
                                                      br(),
                                                      h4("Annotate Results"),
                                                      br(),
                                                      p("This table shows the Gene Ontology Biological Processes terms that are significantly associated with
                                                        the proteins that are significantly represented in the query (Benjamini-Hochberg corrected P <
                                                        0.05.)"),
                                                      h6("go_id: Gene Ontology term ID."),
                                                      h6("GO: Gene Ontology term."),
                                                      h6("proteins: Number of proteins associated with GO term,"),
                                                      h6("F: Fold enrichment of GO term."),
                                                      h6("Lo: Lower bound of fold enrichment range."),
                                                      h6("Hi: Upper bound of fold enrichment range."),
                                                      h6("P: P value of Fisher's exact test of enrichment."),
                                                      h6("BH: Benjamini-Hochberg corrected P value of enrichment."),
                                                      dataTableOutput("enrich_out")
                                                      ),
                                             tabPanel("Step 3: Visualize",
                                                      br(),
                                                      h4("Network Results"),
                                                      br(),
                                                      p("This graph displays the IntAct interacting network of the top 10 proteins"),
                                                      br(),
                                                      h5("Interactome", align="center"),
                                                      #plotOutput("cluster_plot")
                                                      simpleNetworkOutput("go_network")
                                             ),
                                             tabPanel("Step 4: Download",
                                                      br(),
                                                      h4("Download Data"),
                                                      br(),
                                                      p("Select the number of top proteins to download"),
                                                      br(),
                                                      sliderInput("dl_count", label = h4("Number of proteins"), min = 10, max = 260, value = 100),
                                                      br(),
                                                      p("Press the button below to download the results as a csv file."),
                                                      #downloadButton('downloadData', 'Download'),
                                                      br()
                                             )
                                 )
                         )
                 )
        )
        
        

        
}
