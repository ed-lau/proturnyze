###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###
# 
# sidebarLayout(
#         
#         sidebarPanel(
#                 
#         ),
#         
        
refit_page <- function(){

tabPanel("Refit",
         tagline(),

                         tabsetPanel(type ="tabs",
                                     
                                     tabPanel("Step 1: Input Parameters",
                                              br(),
                                              wellPanel(
                                                      h3("Start Here"),
                                                      tags$hr(),
                                                      p("This page displays ProTurn result or refits the ProTurn outputs to a custom kinetic model."),
                                                      br(),
                                                      p("Click here to download sample files for ", tags$a(href = 'hl.out', target="blank","hl.out"), " and ", tags$a(href = 'hl-data.out', target="blank","hl-data.out"),
                                                        "for refitting. The full set of raw data can be found on our", tags$a(href = "http://dx.doi.org/10.7303/syn2289125", target="blank", "Synapse project for the raw protein turnover/expression data"), "."),
                                                      br(),
                                                      p(strong("Note that the refitting process can take a long time. If the server stops responding, simply refresh your browser.")),
                                                      br()
                                                      
                                              ),
                                              br(),
                                              h4("Refit Data"),
                                              h3("Step 1A: Upload ProTurn output files:"),
                                              fluidRow(
                                                      column(4, fileInput('file1', 'Upload hl.out')),
                                                      column(4, fileInput('file2', 'Upload hl-data.out')),
                                                      column(4, fileInput('file3', "Upload an annotation file (optional)"))
                                              ),
                                           
                                              tags$hr(),
                                              h3("Step 1B: Choose whether to refit"),
                                              br(),
                                              fluidRow(
                                                      column(4,
                                                                radioButtons("refit_choice", "Graph or fit? ", 
                                                                c("Graphing Only" = "Graph",
                                                                "Refit Data" = "Fit"),
                                                                selected = "Graph")
                                                      ),
                                                      column(8,
                                                             p(strong("Instructions:")),
                                                             p( "Choose Graph Only if you would just like to visualize an existing ProTurn output (Note: This
                                                               output must come from NS fitting option."),
                                                             p("Choose Refit Data if you would like to refit each indivudal peptide time-series using one of
                                                               the selected models. (Note: For performance reasons, this is currently limited to 500 peptides.)")
                                                             )
                                              ),
                                              hr(),
                                              radioButtons("fitModel", "Choose kinetic model: ", 
                                                           c("Quasi-two-compartment model (Lam et al.)" = "NS",
                                                             "Multi-parameter model (Modified from Kim et al. 2012)" = "SS",
                                                             "Two-compartment model (Modified from Guan et al. 2014)" = "CC"),
                                                        
                                                           selected = "NS"),
                                              hr(),
                                              radioButtons("fitMethod", "Fitting algorithm: ", 
                                                           c("Simplex (Nelder-Mead)" = "Nelder-Mead",
                                                             "Gradient/Quasi-Newton (Broyden-Fletcher-Goldfarb-Shanno))" = "BFGS",
                                                             "Secant/Bisection (Brent-Dekker)" = "Brent"),
                                                           selected = "Nelder-Mead"),
                                              hr(),
                                              actionButton("startRefitButton", "Graph data", icon = icon("play", lib = "glyphicon"), width = "100%"),
                                              textOutput("refitStatus"),
                                              br()
                                     ),
                                     tabPanel("Step 2: Peptide-Centric View",
                                              
                                              h3("Filter input by R2 and DP:"),
                                              fluidRow(
                                                      column(4,
                                                             sliderInput("result_R2filter", "Filter results by coefficient of determination:", min = 0, max = 1, step = 0.01, value = c(0.9,1))),
                                                      column(4,
                                                             sliderInput("result_DPfilter", "Filter results by data points", min = 2, max = 20, step = 1, value = c(2,20)) ),
                                                      column(4,
                                                             sliderInput("result_lenfilter", "Filter results by peptide length", min = 5, max = 30, step = 1, value = c(7,25)))
                                              ),
                                              hr(),
                                              wellPanel(
                                                         h3("Select Proteins"),
                                                         br(),
                                                         DT::dataTableOutput("refitProteinSummary"),
                                                         br()
                                                         ),
                                              br(),
                                              textOutput("refitProteinText"),
                                              hr(),
                                              wellPanel(
                                                      h3("Select Peptides"),
                                                      br(),
                                                      DT::dataTableOutput("refitPeptideSummary"),
                                                      br()
                                              ),
                                              
                                              h3("Plot"),
                                              br(),
                                              fluidRow(
                                                 column(6,
                                                        plotOutput("displayPeptide")),
                                                 column(6,
                                                        plotOutput("displayResidual"))
                                               ),
                                              br()
                                      
                                     ),
                                    tabPanel("Step 3: Protein-Centric View",
                                             h3("Filter input by R2 and DP:"),
                                             fluidRow(
                                                     column(4,
                                                            sliderInput("result_R2filter2", "Filter results by coefficient of determination:", min = 0, max = 1, step = 0.01, value = c(0.9,1))),
                                                     column(4,
                                                            sliderInput("result_DPfilter2", "Filter results by data points", min = 2, max = 20, step = 1, value = c(2,20)) ),
                                                     column(4,
                                                            sliderInput("result_lenfilter2", "Filter results by peptide length", min = 5, max = 30, step = 1, value = c(7,25)))
                                             ),
                                             wellPanel(
                                                     h3("Select Proteins"),
                                                     br(),
                                                     DT::dataTableOutput("refitProteinSummary2"),
                                                     br()
                                             ),
                                             textOutput("refitProteinText2"),
                                             hr(),
                                             wellPanel(
                                                     h3("Select Peptides"),
                                                     br(),
                                                     DT::dataTableOutput("refitPeptideSummary2"),
                                                     br()
                                             )
                                             )
                         )

)
}
