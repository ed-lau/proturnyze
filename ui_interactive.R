
interactive_page <- function() {

        tabPanel("Interactive Fitting",
                tagline(),
                wellPanel(
                        h3("Step 1: Choose Organ and Label"),
                        p("See documentation for details (coming soon)."),
                        fluidRow(
                                column(4,
                                        radioButtons("livOrgan", "Choose an organ to display: ", 
                                        c("Heart" = "h",
                                        "Kidney" = "k",
                                        "Liver" = "l",
                                        "Muscle" = "m"),
                                        selected = "h")),
                                column(4,
                                        radioButtons("livLabel", "Choose a labeling method: ", 
                                        c("Amino acid" = "aa"),
                                        selected = "aa" )),
                                column(4,
                                       radioButtons("livFitMethod", "Choose an optimization algorithm: ", 
                                                    c("Simplex (Nelder-Mead)" = "Nelder-Mead",
                                                      "Gradient/Quasi-Newton (Broyden-Fletcher-Goldfarb-Shanno)" = "BFGS",
                                                      "Secant/Bisection (Brent-Dekker)" = "Brent"),
                                                    selected = "Nelder-Mead"))
                        )
                
                ),
                hr(),
                wellPanel(h3("Step 2: Select the peptide to be analyzed"),
                        fluidRow(
                                column(4,
                                       br(),
                                       sliderInput("livFilter_num_k", "Filter peptides by the number of lysines", min = 0, max = 5, step = 1, value = c(1,2)),
                                       hr(),
                                       sliderInput("livFilter_DP", "Filter peptides by the number of data points", min = 2, max = 20, step = 1, value = c(5,15)),
                                       hr(),
                                       sliderInput("livFilter_len", "Filter peptides by sequence length", min = 5, max = 40, step = 1, value = c(7,25)),
                                       hr(),
                                       sliderInput("livFilter_ess", "Filter peptides by sequence amino acid essentiality", min = 0, max = 1, step = 0.01, value = c(0,1)),
                                       br()
                                ),
                                column(8,
                                       DT::dataTableOutput("livDisplayPeptides"),
                                        textOutput("livDisplayPeptideText")
                                )
                        )
                ),
                hr(),
                wellPanel(
                        
                        h3("Step 3: Optimize for kinetic model and visualize results"),
                
                
                        fluidRow(
                                column(4,
                                       br(),
                                       selectInput("livPlotType", label="Choose a plot to display", 
                                                   choices = list("Best-fit curve"=1, "Residual plot"=2), 
                                                   selected=1),
                                       hr(),
                                       sliderInput("livFitting_kp", "Select precursor turnover rate:", min = 0.01, max = 2, step = 0.01, value = 0.3),
                                       hr(),
                                       sliderInput("livFitting_pss", "Select plateau precursor enrichment:", min = 0.01, max = 0.5, step = 0.01, value = 0.4),
                                       br()
                                ),
                                column(4,
                                       plotOutput("livDisplayPeptide_SS1")),
                                column(4,
                                       plotOutput("livDisplayPeptide_CC1"))
                              )
                        )
                )
                

}