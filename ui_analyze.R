###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###

# Merged into Browse tab.

analyze_page <- function(){

        
        tabPanel("Analyze",
                 tagline(),
                 sidebarLayout(
                         sidebarPanel(
                                 h3("Dataset Analysis"),
                                 tags$hr(),
                        
                                 br(),
                                 p(" Here you can explore one way to analyze the data you saw in the previous tab. To save time, we have tranformed 
                                        and standardized the data,
                                     so that the three data types (protein turnover, protein abundance, transcript abundance) may be analyzed together. 
                                   One way to do so is through clustering analysis. Clustering is a method to group similar proteins together based on 
                                   their behaviors in transcript and protein parameters."),
                                 br(),
                                 p("To start, simply go to the ", code("Step 1: Input"), " tab and select your desired parameters. After you are satifisied with 
                                the selections, go to", code("Step 2: Output"), " and press the Run button to begin clustering."),
                                 br(),
                                 p("The method here uses a model-based clustering method implemented in the ", code("Mclust"), "package in R. For Hacker Track 
                                        participants planning to analyze the data by other means,
                                   click here to ", tags$a(href = 'trinity.txt', target="blank","download the dataset"), "as a tab-separated text file."),
                                 br()
                                 ), 
                         mainPanel(
                                 tabsetPanel(type ="tabs",
                                             
                                      tabPanel("Step 1: Input",

                                               checkboxGroupInput("datatypeChoice", label = h3("Step 1A: Select Data Type"), 
                                                                 choices = list("Protein turnover" = 1, "Protein abundance" = 2, "Transcript abundance" = 3),
                                                                 selected = c(1)),
                                               p("Select the omics data type(s) to be included in the data analysis."),
                                               tags$hr(),
                                               checkboxGroupInput("strainChoice", label = h3("Step 1B: Select Mouse Strain"), 
                                                                  choices = list("A/J" = 1, "BALB/cJ" = 2, "C57BL/6/J" = 3, "CE/J" = 4, "DBA/2J" = 5, "FVB/NJ" = 6),
                                                                  selected = c(1,2,3,4,5,6)),
                                               p("Select the mouse strain(s) to be included in the data analysis."),
                                               tags$hr(),
                                               checkboxGroupInput("factorChoice", label = h3("Step 1C: Select Disease Factor"), 
                                                                  choices = list("Normal" = 1, "Hypertrophy" = 2),
                                                                  selected = c(1,2)),
                                               p("Select the disease factor(s) to be included in the data analysis."),
                                               tags$hr(),
                                               h3("Step 1D: What is next?"),
                                               textOutput("filterDemoDataFeedback"),
                                               br(),
                                               br(),
                                               br()
                                               
                                      ),
                                      tabPanel("Step 2: Output",
                                               h3("Step 2A: Press Run!"),
                                               br(),
                                               actionButton("startClusterButton", "Run!", icon = icon("play", lib = "glyphicon"), width = "100%"),
                                               br(),
                                               p("This will start clustering and may take a minute."),
                                               tags$hr(),
                                               h3("Step 2B: View Clustering Result"),
                                               plotOutput("renderDemoClusters"),
                                               tags$hr(),
                                               h3("Step 2C: Visualize Individual Clusters"),
                                               tags$hr(),
                                               selectInput("pickCluster", label="Choose cluster to display:", choices=c(1), selected = 1),
                                               selectInput("pickPlotType", label="Choose plot type:", choices=c("Box plot" = "boxplot", "Line plot" = "lineplot"), selected = "boxplot"),
                                               plotlyOutput("renderDemoSingleCluster"),
                                               h3("Step 2D: Download Cluster Protein List"),
                                               p("Clicking on the two buttons below will download the list of proteins present in the selected cluster above, as well as the 
                                                 background list of all analyzed proteins, respectively. The protein lists can be analyzed for their biological significance 
                                                 in the ", code("Enrich"), " tab. See below for more details."),
                                               downloadButton("downloadCluster", label = "Download Proteins in Cluster"),
                                               downloadButton("downloadBackground", label = "Download Background Proteins"),
                                               tags$hr(),
                                               h3("Step 2E: What is next?"),
                                               br(),
                                               p("Here are some questions we can ask with the clustered data. E.g., Do the clusters co-segregate with certain biological 
                                                 pathways? Are proteins involved in a particular biological function preferentially increased in their turnover during 
                                                 cardiac hypertrophy?"),
                                               br(),
                                               p("This can be answered by doing a functional analysis of the clusters we got from the analysis on this page, using sources
                                                 of biological annotations such as Gene Ontology or Reactome. Under the ", code("Enrich"), " tab you will find a simple tool that performs 
                                                 one such analysis, allowing us to compare a protein list against a background and look for over-representation of biological
                                                 functions."),
                                               br(),
                                               p("Hacker Track Participants: If you are already familiar with the process, here are some other links where you may visit to explore other analysis routes."),
                                               p(tags$a(href = 'http://reactome.org', target="blank","Reactome pathway analysis of gene lists")),
                                               p(tags$a(href = "http://dx.doi.org/10.7303/syn2289125", target="blank", "Synapse project for the raw protein turnover/expression data")),
                                               p(tags$a(href = "http://heartproteome.org/copa/ProteinTurnover.aspx", target="blank", "Instructions for analyzing the Synapse data using R"))
                                               )
                                 )
                         )
                 )
        )
        
     
}
