###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###


documentation_page <- function(){
        
        
        tabPanel("Documentation",
                 tagline(),
                 sidebarLayout(
                         sidebarPanel(
                                 h3("Welcome"),
                                 tags$hr(),
                                 p("IsoVista allows users to visualize stable isotope enriched protein turnover experiment data."),
                                 br(),
                                 tags$hr(),
                                 h4("About"),
                                 p("Server Version:"),
                                 p(code("1.1.1")),
                                 br(),
                                 p("Annotation Updated:"),
                                 p(code(textOutput("annotUpdateDate"))),
                                 br(),
                                 p("Contact:"),
                                 p("Edward Lau - edward.lau@me.com")
                                 ), 
                         mainPanel(
                                 h3("Introduction"),
                                 p("IsoVista allows users to visualize stable isotope enriched protein turnover experiment data. Follow
                                   the instructions on the ", code("Analysis > Peptide-Centric Fitting"), " tab to get started. Protein-Centric 
                                   Fitting function will be implemented in the near future. "),
                                 hr(),
                                 h3("Data format requirements"),
                                 p("IsoVista comes with preloaded amino acid labeling data files. You may override the preloaded dataset
                                   with custom data files. To do so, two tab-separated unquoted text files are required."),
                                 br(),
                                 p("The Peptide table must contain the following columns with the exact headers:"),
                                 p(tags$b("ID: "), " a unique integer ID for the data entry (e.g., 1)"),
                                 p(tags$b("Uniprot: "), "a string containing the Uniprot accession (e.g., P68033)"),
                                 p(tags$b("Peptide: "), "a string containing the peptide sequence (e.g., DLTDYLMK)"),
                                 p(tags$b("DP: "), "an integer denoting the number of time points from which the peptide was quantified (e.g., 11)"),
                                 p(tags$b("z: "), "an integer denoting the charge of the peptide (e.g., 2)"),
                                 br(),
                                 p("The isotopomer table must contain the following columnes with the exact headers:"),
                                 p(tags$b("ID: "), "an integer ID that maps the isotopomer entry to the ID of a peptide in the peptide table"),
                                 p(tags$b("t: "), "a floating point number denoting the time point (in days) from which the isotopomer was measured"),
                                 p(tags$b("A0: "), "a floating point number denoting the RIA of this particular time point. For amino acid labeling experiment this is currently hard-coded
                                   to represent heavy/(light+heavy) ratio."),
                                 br(),
                                 p("Download sample ", tags$a(href = "small-pep.txt", target="blank", "peptide"), 
                                   " and ", tags$a(href = "small-iso.txt", target="blank", "isotopomer"), "data here."),
                                 hr(),
                                 h3("Kinetic Models"),
                                 p("The steady-state model follows classic first-order kinetics, with slight
                                   modifications in implementation for heavy water labeling noted in Kim et al. Mol Cell Proteomics. 2012."),
                                 uiOutput("SS_Model_Mathjax"),
                                 p("The two-compartment model is adapted
                                   from Guan et al. Anal Chem. 2012, with slight modifications noted in Lam et al. J Clin Invest. 2014 for heavy water labeling."),
                                 uiOutput("CC_Model_Mathjax"),
                                 p("Note on the the maximum attainable isotopomer fractional abundance: Amax is equal to the plateau precursor enrichment in the
                                   amino acid labeling method. In the heavy water method it is calculated as the combinatorial probability of the plateau precursor
                                   label to the power of the number of deuterium accessible and irreversible labeling sites on the peptide sequence."),
                                 hr(),
                                 h3("Peptide Essentiality"),
                                 p("Esssentiality is calculated as the average of synthetic flux of each amino in the sequence
                                   based on Reeds J Nutr. 2000."),
                                 hr(),
                                 h3("References"),
                                 p(strong("Additional information on kinetic models for heavy water labeling:")),
                                 p("Lam MP et al., Protein kinetic signature of the remodeling heart following isoproterenol stimulation. Journal of Clinical Investigation; 124(4): 1734 (2014)",
                                   tags$a(href = "http://www.jci.org/articles/view/73787", target="blank", "Link")),
                                 br(),
                                 p(strong("Additional datasets for heavy water labeling:")),
                                 p("Lau E et al., A large dataset of protein dynamics in the mammalian heart proteome. Scientific Data, 3: 160015 (2016).",
                                   tags$a(href = "http://www.nature.com/articles/sdata201615", target="blank", "Link")),
                                 br(),
                                 p(strong("Additional information on compartment modeling:")),
                                 p("Guan S et al., Compartment modeling for mammalian protein turnover studies by stable isotope metabolic labeling. Anal Chem; 84(9): 4012 (2012)",
                                   tags$a(href = "http://www.ncbi.nlm.nih.gov/pubmed/22444387", target="blank", "Link")),
                                 br(),
                                 p(strong("Additional information on amino acid essentiality:")),
                                 p("Reeds PJ, Dispensable and indispensable amino acids for humans. J Nutr; 130(7): 1835 (2000)",
                                   tags$a(href = "http://www.ncbi.nlm.nih.gov/pubmed/10867060", target="blank", "Link")),
                                 br(),
                                 hr()
                         )
                                 )
                 )
        
        
}
