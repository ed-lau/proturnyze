###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###


instructions_page <- function(){
        
        
        tabPanel("Instructions",
                 tagline(),
                 sidebarLayout(
                         sidebarPanel(
                                 h3("Welcome"),
                                 tags$hr(),
                                 p("Proturynze is an online extension of the protein turnover analysis software ProTurn. This 
                                   customized version is built to demonstrate interactive data analysis."),
                                 br(),
                                 tags$hr(),
                                 h4("About"),
                                 p("Server Version:"),
                                 p(code("1.0.01")),
                                 br(),
                                 p("Gene Ontology Updated:"),
                                 p(code(textOutput("GOUpdateDate"))),
                                 br()
                                 ), 
                         mainPanel(
                                 h3("Step 1: Instructions"),
                                 p("Here you can view a cardiac multi-omics dataset and interact with it. During this Hackathon we are also hosting a
                                  machine learning practicum where you will learn more about application of simple machine learning methods.
                                    Here you can interact with a relatively small dataset and perform simple data-driven discovery."),
                                 br(),
                                 p("There are many routes to analyzing these data and generating new hypotheses from the process. For example, 
                                   our investigators have learned from this dataset the following about the heart:"),
                                 p("- The half-life of a protein may be negatively correlated with its steady-state abundance in the heart."),
                                 p("- Proteins within the mitochondria have different half-lives that span almost two orders of magnitude."),
                                 p("- Members of cardiac biochemical pathways appear to have clustered protein turnover rates."),
                                 br(),
                                 p("For beginners: Our goal here is to try to formulate one such similar hypothesis, or extend one of the above, as guided by the data.
                                   At the end of the section, you can present your one-sentence claim or hypothesis about the data to the group"),
                                 br(),
                                 p("For experienced participants: The raw dataset can be found under the ", code("Analyze"), "tab. Show us what you can do with the data! Further
                                   data can be found on our repository on ", tags$a(href = "http://dx.doi.org/10.7303/syn2289125", target="blank", "Synapse project"), " where you can download the raw data and upload updated files."),
                                 tags$hr(),
                                 h3("Step 2: What is next?"),
                                 p("Whenever you are ready, feel free to go to the",code("Explore"), "tab. There you can follow further instructions
                                  to browse the dataset and also generate some simple exploratory graphs."),
                                 br(),
                                 p("If you have any questions, please don't hesitate to go to the dedicated chat room in the ", code("Chat"), " tab and ask them."),
                                 tags$hr(),
                                 h3("Dataset Citation"),
                                 p(strong("Protein turnover and abundance data are from the following studies:")),
                                 p("Lam et al., Protein kinetic signature of the remodeling heart following isoproterenol stimulation. Journal of Clinical Investigation; 124(4): 1734 (2014)", 
                                   tags$a(href = "http://www.jci.org/articles/view/73787", target="blank", "Link")),
                                 br(),
                                 p("Lau et al., A large dataset of protein dynamics in the mammalian heart proteome. Scientific Data, 3: 160015 (2016).",  
                                   tags$a(href = "http://www.nature.com/articles/sdata201615", target="blank", "Link")),
                                 br(),
                                 p(strong("Transcript abundance data are from the following dataset:")),
                                 p("Rau et al., Transcriptomes of the hybrid mouse diversity panel subjected to isoproterenol challenge. NCBI GEO: GSE48670 (2015).", 
                                   tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48670", target="blank", "Link"))
                         )
                                 )
                 )
        
        
}
