
interactive_page <- function() {

        tabPanel("Analysis",
                tagline(),
                wellPanel(
                        h3("Step 1: Choose Organ and Label"),
                        p("Choose a pre-loaded dataset (Heart, Kidney, Liver, Muscle) or upload a custom dataset. Uploaded datasets
                          need to follow specific format requirements - see ", code("Documentation"), " page for details."),
                        fluidRow(
                                column(4,
                                        radioButtons("livOrgan", "Choose a pre-loaded dataset:", 
                                        c("Heart" = "h",
                                        "Kidney" = "k",
                                        "Liver" = "l",
                                        "Muscle" = "m"),
                                        selected = "h")),
                                
                                column(4,
                                       p(tags$b("Alternatively, upload a custom dataset:")),
                                       fileInput('file1', 'Peptide'),
                                       fileInput('file2', 'Isotopomer')),
                                
                                column(4,
                                       radioButtons("livLabel", "Choose metabolic labeling method: ", 
                                                    c("Amino acid (13C6 Lysine)" = "aa",
                                                      "Heavy water (2H)" = "hw"),
                                                    selected = "aa" ))
                                )
                
                ),
                hr(),
                tabsetPanel(type ="tabs",
                    tabPanel(h4("Option A: Peptide-Centric View"),
                        br(),
                        wellPanel(h3("Step 2: Select the peptide to be analyzed"),
                                  p("Use the dials on the left to filter peptides from the dataset. Then select a peptide on the table on the right
                                    to display its isotope labeling profile."),
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
                                p("Display the best-fitted curve or residual plot for the selected peptide based on two kinetic models. 
                                  Use the numeric input fields to adjust precursor enrichment rates and plateau."),
                                fluidRow(
                                        column(4,
                                               br(),
                                               selectInput("livPlotType", label="Choose a plot to display", 
                                                           choices = list("Best-fit curve"=1, "Residual plot"=2), 
                                                           selected=1),
                                               br(),
                                               checkboxInput("livToggleResidual", label="Toggle residual drop line", value = T),
                                               hr(),
                                               
                                               # Changing to numerical input per request
                                               #sliderInput("livFitting_kp", "Select precursor turnover rate:", min = 0.01, max = 2, step = 0.01, value = 0.3),
                                               numericInput("livFitting_kp", label = "Select precursor enrichment rate:", value = 0.3, min=0.01, max=10, step=0.01),
                                               hr(),
                                               # Changing to numerical input per request
                                               #sliderInput("livFitting_pss", "Select plateau precursor enrichment:", min = 0.01, max = 0.5, step = 0.01, value = 0.4),
                                               numericInput("livFitting_pss", label = "Select precursor enrichment plateau:", value=0.4, min=0.01, max=1, step=0.01),
                                               hr(),
                                               radioButtons("livFitMethod", "Select optimization algorithm: ", 
                                                            c("Simplex (Nelder-Mead)" = "Nelder-Mead",
                                                              "Gradient/Quasi-Newton (Broyden-Fletcher-Goldfarb-Shanno)" = "BFGS",
                                                              "Secant/Bisection (Brent-Dekker)" = "Brent"),
                                                            selected = "Nelder-Mead"),
                                               br()
                                        ),
                                        column(4,
                                               plotOutput("livDisplayPeptide_SS1")),
                                        column(4,
                                               plotOutput("livDisplayPeptide_CC1"))
                                      )
                                )
                        ),
                        tabPanel(h4("Option B: Protein-Centric View"),
                                br(),
                                wellPanel(h3("Step 2: Select the protein to be analyzed"),
                                          fluidRow(
                                                  column(4,
                                                         br(),
                                                         sliderInput("livFilter_num_pep", "Filter proteins by the number of peptides", min = 1, max = 1000, step = 1, value = c(1,1000)),
                                                         br()
                                                  ),
                                                  column(8,
                                                         DT::dataTableOutput("livDisplayProteins"),
                                                         textOutput("livDisplayProteinText")
                                                  )
                                          )
                                ),
                                
                                wellPanel(h3("Step 3: Select one or more peptides within a protein to be analyzed"),
                                          fluidRow(
                                                  column(4,
                                                         br(),
                                                         # For the peptides within the selected protein, how many data points are they found in.
                                                         sliderInput("livFilter_pro_pep_DP", "Filter peptides within selected protein 
                                                                     by DP", min = 2, max = 20, step = 1, value = c(5,15)),
                                                         br()
                                                  ),
                                                  column(8,
                                                         DT::dataTableOutput("livDisplayProteinPeptides"),
                                                         textOutput("livDisplayProteinPeptideText"),
                                                         br()
                                                  )
                                          )
                                ),
                                
                                wellPanel(h3("Step 4: Optimize for kinetic model and visualize results"),
                                          fluidRow(
                                                  column(4,
                                                         br(),
                                                         p("Coming soon"),
                                                         #sliderInput("livFilter_pep_R2_CC", "Filter peptides by R2", min = 1, max = 1000, step = 1, value = c(1,1000)),
                                                         br()
                                                  ),
                                                  column(4,
                                                          #plotOutput("livDisplayProteinPeptide_SS1")
                                                          br()
                                                  ),
                                                  column(4,
                                                         plotOutput("livDisplayProteinPeptide_CC1"))
                                                  )
                                          
                                ),
                                hr()
                        )
                )
        )

}