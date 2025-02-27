
######
######
###### Equations for Documentation
######
######

output$SS_Model_Mathjax <- renderUI({
        withMathJax(
                helpText('First-order kinetic model
                         $$ a + (A_{max} - a) \\cdot \\left( 1 - e^{-k \\cdot t} \\right) $$'))
})

output$CC_Model_Mathjax <- renderUI({
        withMathJax(
                helpText('Simple two-compartment kinetic model
                         $$ a + (A_{max} - a) \\cdot \\left( 1 - \\frac{ e^{-k \\cdot t} \\cdot k_p - 
                         e^{-k_p \\cdot t} \\cdot k}{k_p - k} \\right)$$'))
})



######
######
######  Peptide-Centric Interactive Refitting Functions
######
######


livRefitSinglePeptideFunction <- function(hl, 
                                          hldata, 
                                          #row_id, 
                                          fitMethod="Nelder-Mead", 
                                          fitModel="SteadyStateOneParameter", 
                                          labelMethod="aa",
                                          pss,
                                          kp) {
        
        
        SS1_Refit <- function(var){
                Refitting_Predicted <- sapply(ds$t, function(x) SSModel(x, var, a, Amax))
                Refitting_SS <- sum((ds$A0-Refitting_Predicted)^2)
                return(Refitting_SS)
        }
        
        # var is kdeg. ksyn will have to be supplied; a is beginning and Amax is plateau
        CC1_Refit <- function(var){
                Refitting_Predicted <- sapply(ds$t, function(x) CCModel(x, var, kp, a, Amax))
                Refitting_SS <- sum((ds$A0-Refitting_Predicted)^2)
                return(Refitting_SS)
        }
        

        #
        #
        # This should now loop through all the peptide IDs in the sliced peptide table
        # And do fitting for each and every one ID
        #
        #
        
        # Get a list of all the IDs present in the hl file
        ids = hl %>% dplyr::distinct(ID) %>% arrange(ID) %>% unlist()
   
        print(ids)     
        # How many IDs are present?
        num_ids = length(ids)
        
        withProgress (message = 'Refitting Data...', value = 0.5, {
   
        output_for_all_id <- list()
        
        for (i in 1:length(ids)){
                incProgress(0.5/num_ids)
                        
                row_id <- ids[i]
                
                # Slicing hl and hldata.out by the ID being considered
                hl. <- hl %>% filter(ID == row_id) %>% slice(1)
                ds <- hldata %>% filter(ID == row_id)
                
                # Get the correct kp, Amax, and a value
                if(labelMethod == "aa"){
                        Amax <- pss
                        a <- 0
                } else if (labelMethod == "hw") {
                        seq <- hl.$Peptide[1] %>% unlist() %>% as.character()
                        N <- calcLabelingSite(seq)[1]
                        a <- calc_thr_a(seq) 
                        Amax <- a* ((1 - pss)^N) 
                }
                
                # LOOP VERSION
                if(fitModel == "SteadyStateOneParameter"){

                        if(fitMethod == "Brent"){
                                Optimize <- optim(0.249, SS1_Refit, method = "Brent", lower = 0, upper = 4, control = list(maxit=15))
                        } else {
                                Optimize <- optim(0.249, SS1_Refit, method = fitMethod, control = list(maxit=100))
                        }
                        
                        k <- Optimize$par
                        SS <- Optimize$value
                        SE <- (Optimize$value/(hl.$DP-1))^0.5
                        dk <- sapply(ds$t, function(x) SSModel_dk(x, Optimize$par, a, Amax, SE)) %>% abs() %>% min()
                        
                }
                
                if(fitModel == "TwoCompartmentOneParameter"){

                        if(fitMethod == "Brent"){
                                Optimize <- optim(0.249, CC1_Refit, method = "Brent", lower = 0, upper = 4, control = list(maxit=15))
                        } else {
                                Optimize <- optim(0.249, CC1_Refit, method = fitMethod, control = list(maxit=100))
                        }

                        k <- Optimize$par
                        SS <- Optimize$value
                        SE <- (Optimize$value/(hl.$DP-1))^0.5
                        dk <- sapply(ds$t, function(x) CCModel_dk(x, Optimize$par, kp, a, Amax, SE)) %>% abs() %>% min()
                        
                }
                
                
                R2 <- 1- (SS/(sum((ds$A0 - mean(ds$A0))^2)))
                
                output_for_this_id <- list(round(k, 3),
                                           round(dk, 3),
                                           round(R2, 3))
                
                # Merge all the outputs from each ID into a nested list.
                output_for_all_id[[length(output_for_all_id)+1]] <- output_for_this_id

        }
        
        
        })
       
        output_for_all_id
        
        
}



######
######
###### Interactive Refit SERVER reactive functions
######
######

livRead_hl <- reactive ({
        
       
        withProgress(message = 'Reading Data...', value = 0.2, {
                
        
        organ <- input$livOrgan 
        label_type <- input$livLabel 
        file_path <- paste0("data/", label_type, "_", organ, ".txt")
            
        # Override with user uploaded file    
        if (!is.null(input$file1) & !is.null(input$file2)) {file_path <- input$file1$datapath}
                
        hl <- read.table(file_path, header=T, fill=T, as.is=T, sep="\t", quote="")
        
        incProgress(0.5)
        
        # Row-wise is hella slow. Use Mutate instead?
        
        # hl <- hl %>% group_by(ID) %>% do({
        #         concat = paste0(.$Peptide, "[", .$z, "+]") 
        #         #N = calcLabelingSite(.$Peptide)[1]
        #         #a = calc_thr_a(.$Peptide)
        #         #Amax = a * (1 - .$pss)^N
        #         ess = calc_essentiality(.$Peptide)
        #         len = nchar(gsub("\\(.*\\)","", .cal$Peptide))
        #         num_k = stringr::str_count(.$Peptide, "K")
        #         data.frame(., concat, len, num_k, ess)
        #         
        # }) %>% ungroup()
        
        
        hl <- hl %>% ungroup() %>% dplyr::mutate(concat = paste0(Peptide, "[", z, "+]"),
                                   #N = calcLabelingSite(Peptide)[1],
                                   #a = calc_thr_a(Peptide),
                                   #Amax = a * (1 - pss)^N,
                                   len = nchar(gsub("\\(.*\\)","", Peptide)),
                                   ess = sapply(Peptide, calc_essentiality) %>% round(2),
                                   num_k = stringr::str_count(Peptide, "K"))
        
                                })
        
        
        hl 
        
})  
        

livRead_hldata <- reactive({

        organ <- input$livOrgan
        label_type <- input$livLabel
        file_path <- paste0("data/", label_type, "_", organ, "_data.txt")
        
        # Override with user uploaded file
        if (!is.null(input$file1) & !is.null(input$file2)) {file_path <- input$file2$datapath}
        hldata <- read.table(file_path, header=T, fill=T, as.is=T, sep="\t", quote="")
        hldata
        
})   



# Display the peptides from the chosen dataset
livFilterPeptides <- reactive({
        
        hl <- livRead_hl()
        if(is.null(hl)) return(NULL)

        hl <- hl %>% dplyr::filter(num_k >= input$livFilter_num_k[1],
                                   num_k <= input$livFilter_num_k[2],
                                   len >= input$livFilter_len[1],
                                   len <= input$livFilter_len[2],
                                   DP >= input$livFilter_DP[1],
                                   DP <= input$livFilter_DP[2],
                                   ess >= input$livFilter_ess[1],
                                   ess <= input$livFilter_ess[2])

        hl
        
})



# Display the peptides on a data table
output$livDisplayPeptides <- DT::renderDataTable({
        
        if(!is.null(livFilterPeptides())) livFilterPeptides() %>% dplyr::select(ID, Uniprot, Peptide, z, DP, ess)
        
        }, selection="single",  options = list(pageLength=8,
                                               lengthChange=F,
                                               filter="top")
        )   


# Display the peptide table's summary statistics
output$livDisplayPeptideText <- renderText({
        if(!is.null(livFilterPeptides())){
                filtered_data <- livFilterPeptides()
                paste0("After filtering, ", nrow(filtered_data), " peptides from ", nrow(distinct(filtered_data, Uniprot)), " proteins remain.")
        }
        
})

# The core function to display a fitted peptide or the residual. All plots should call this from wrappers.
# This is not a reactive so all variables must be passed as arguments

livPlotPeptideFunction <- function (hl, 
                                    hldata, 
                                    labelMethod, 
                                    plotModel="TwoCompartment", 
                                    plotType, 
                                    toggleResidual, 
                                    pss, 
                                    kp, 
                                    fit,
                                    isProtein=F,
                                    R2filter=c(0,1)) {
        
        
        # Get a list of all the IDs present in the hl file
        ids = hl %>% dplyr::distinct(ID) %>% arrange(ID) %>% unlist()
        
        print(ids)
        print(fit)
        # How many IDs are present?
        num_ids = length(ids)
        
        
        # Get the correct kp, Amax, and a value
        if(labelMethod == "aa"){
                Amax <- pss
                a <- 0
        } else if (labelMethod == "hw") {
                seq <- hl$Peptide[1] %>% unlist() %>% as.character()
                N <- calcLabelingSite(seq)[1]
                a <- calc_thr_a(seq) 
                Amax <- a* ((1 - pss)^N) 
        }
        
        
        
        # Get residuals from hldata.
        # Match each row's ID to the proper ID in the fit result list
        # then use that k to get model and residual
        
        hldata$predicted <- 0
        hldata$residual <- 0
        for (i in 1:length(ids)) {
                
                # Get the optimized k
                k = fit[i] %>% unlist() %>% (function(x) x[1])
                
                for (j in 1:nrow(hldata)){
                        if (hldata$ID[j] == ids[i]){
                                if (plotModel == "SteadyState"){
                                        hldata$predicted[j] <- SSModel(hldata$t[j], k, a, Amax)
                                } else if (plotModel == "TwoCompartment"){
                                        hldata$predicted[j] <- CCModel(hldata$t[j], k, kp, a, Amax)
                                } 
                                hldata$residual[j] <- hldata$A0[j] -  hldata$predicted[j]
                        }
                }
        }
        
               
        # Prepare GGplot 
        g <- ggplot(data = hldata, aes(x = t, y = A0))
        

        withProgress (message = 'Plotting Graph...', value = 0.5, {
        # Draw each of the fitted peptide
        for (i in 1:length(ids)){
                incProgress(0.5/length(ids))
                
                # Get the optimized k, dk, and R2
                k = fit[i] %>% unlist() %>% (function(x) x[1])
                dk= fit[i] %>% unlist() %>% (function(x) x[2])
                R2= fit[i] %>% unlist() %>% (function(x) x[3])
                
                
                if ((R2 > R2filter[1] & R2 < R2filter[2]) | isProtein==F){
                g <- g + geom_point(data = hldata %>% filter(ID %in% ids[i]), col="black") 
                } else {
                g <- g + geom_point(data = hldata %>% filter(ID %in% ids[i]), col="grey")       
                }
                
                if (plotModel == "SteadyState"){
                        
                        # If this is a protein plot, only add the line when the R2 is within filter
                        if ((R2 > R2filter[1] & R2 < R2filter[2]) | isProtein==F){
                        # The "regular" stat_function syntax (fun= , args=, etc.) is required if we want to loop
                        g <- g + stat_function(fun=SSModel,
                                               args=list(k=k, 
                                                         a=a, 
                                                         Amax=Amax),
                                               alpha=0.33
                                               )
                        }
                        print(paste0("added ggplot for ", i, " of ", num_ids))

                        
                        fitTitle = paste0("Protein-Level First-Order Kinetics Curve \n of ", num_ids, " Input Peptide(s)")
                        # Only add fitting error lines if displaying one peptide
                        if(isProtein==F & num_ids == 1){
                                g <- g + stat_function(fun = function(x) SSModel(x, k+dk, a, Amax), color = "red")
                                g <- g + stat_function(fun = function(x) SSModel(x, k^2/(k+dk), a, Amax), color="red")
                                fitTitle <- paste0("First-Order Kinetics Model \n", "Peptide: ", hl$concat, "\n k: ", k, " R2: ", R2)
                                residualTitle = "First-Order Kinetics Model - Residual Plot"
                        }
                        
                } 
                
                if (plotModel == "TwoCompartment") {
                        
                        # If this is a protein plot, only add the line when the R2 is within filter
                        if ((R2 > R2filter[1] & R2 < R2filter[2]) | isProtein==F){
                        # The "regular" stat_function syntax (fun= , args=, etc.) is required if we want to loop
                        g <- g + stat_function(fun = CCModel,
                                               args=list(kdeg=k, 
                                                         ksyn=kp, 
                                                         a=a, 
                                                         Amax=Amax),
                                               alpha=0.33
                                               )
                        }
                        
                        print(paste0("added ggplot for ", i, " of ", num_ids))

                        
                        fitTitle = paste0("Protein-Level 2-Compartment Kinetics Curve \n of ", num_ids, " Input Peptide(s)")
                        # Only add fitting error lines if displaying one peptide
                        if(isProtein == F & num_ids == 1){
                                g <- g + stat_function(fun = function(x) CCModel(x, k+dk, kp, a, Amax), color = "red")
                                g <- g + stat_function(fun = function(x) CCModel(x, k^2/(k+dk), kp, a, Amax), color="red")
                                
                                fitTitle <- paste0("Two Compartment Kinetics Model \n", "Peptide: ", hl$concat, "\n k: ", k, " R2: ", R2)
                                residualTitle = "Two Compartment Kinetics Model - Residual Plot"
                        }
                }
        
                #Dropping a line from the data point to the curve - only do this if there is one peptide selected or else it gets too busy.
                if(toggleResidual == T & num_ids == 1){
                        g <- g + geom_segment(aes(x=t, xend=t, y=predicted, yend=A0), col="grey")
                }
                
        }
        })
        
        #Find x limit
        gx <- max(ggplot_build(g)$data[[1]]$x)
        gy <- max(ggplot_build(g)$data[[1]]$y)
        
        g <- g + xlim(0, gx) + ggtitle(fitTitle)
        
        # Plot residual plot instead of best-fit curve
        if (plotType == 2){
                g <- ggplot(data = hldata, aes(x = t, y = residual))
                g <- g + geom_point() + ggtitle(residualTitle)
                g <- g + geom_hline(yintercept=0, col="black")
                g <- g + geom_hline(yintercept= c(abs(max(hldata$A0)-min(hldata$A0))/2,-abs(max(hldata$A0)-min(hldata$A0))/2), col="red")
        }
        
        g <- g + theme(aspect.ratio=1,    
                       #panel.grid.minor = element_blank(), 
                       #panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       plot.background = element_blank()) 
        #g <- ggplotly(g)
        g
}
        
        


# Display individual peptides from the refitting results (STEADY STATE GRAPH)
output$livDisplayPeptide_SS1 <- renderPlot({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(!nrow(hl)) return(NULL)

        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- isolate(livRead_hldata()) %>% filter(ID %in% hl$ID)

        # Get optimized k
        fit <- livRefitSinglePeptideFunction(hl=hl, 
                                      hldata=hldata, 
                                      #row_id=row_id, 
                                      fitMethod=input$livFitMethod,
                                      fitModel="SteadyStateOneParameter",
                                      labelMethod=isolate(input$livLabel),
                                      pss=input$livFitting_pss,
                                      kp=input$livFitting_kp)
        
        # Display GG Plot using the fitted values
        livPlotPeptideFunction(hl=hl,
                           hldata=hldata,
                           labelMethod=input$livLabel,
                           plotModel="SteadyState",
                           plotType=input$livPlotType,
                           toggleResidual=input$livToggleResidual,
                           pss=input$livFitting_pss,
                           kp=input$livFitting_kp,
                           fit=fit,
                           isProtein=F,
                           R2filter=c(0,1))
        
}, bg="transparent")

# Display individual peptides from the refitting results (Two-Compartment Graph)
output$livDisplayPeptide_CC1 <- renderPlot({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(!nrow(hl)) return(NULL)

        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- livRead_hldata() %>% filter(ID %in% hl$ID)

        # Get the fitted values
        # Get optimized k
        fit <- livRefitSinglePeptideFunction(hl=hl, 
                                             hldata=hldata, 
                                             #row_id=row_id, 
                                             fitMethod=input$livFitMethod,
                                             fitModel="TwoCompartmentOneParameter",
                                             labelMethod=isolate(input$livLabel),
                                             pss=input$livFitting_pss,
                                             kp=input$livFitting_kp)

        # Display GG Plot using the fitted values
        livPlotPeptideFunction(hl=hl,
                           hldata=hldata,
                           labelMethod=input$livLabel,
                           plotModel="TwoCompartment",
                           plotType=input$livPlotType,
                           toggleResidual=input$livToggleResidual,
                           pss=input$livFitting_pss,
                           kp=input$livFitting_kp,
                           fit=fit,
                           isProtein=F,
                           R2filter=c(0,1))
        
}, bg="transparent")



######
######
###### Protein Centric View functions below
######
######

# Display the protein-centric view from the chosen dataset
livFilterProteins <- reactive({
        
        hl <- livRead_hl()
        if(is.null(hl)) return(NULL)
        
        
        hl <- hl %>% group_by(Uniprot) %>% summarize(num_pep = n()) %>% ungroup()
        
        
        hl <- hl %>% dplyr::filter(num_pep >= input$livFilter_num_pep[1],
                                      num_pep <= input$livFilter_num_pep[2])
        
        hl
        
})

# Display the proteins from the dataset on a data table
output$livDisplayProteins <- DT::renderDataTable({
        
        if(!is.null(livFilterProteins())) {
                
                withProgress (message = 'Retrieving Annotations...', value = 0.5, {
                annot <- annotations()
                incProgress(0.3)
                livFilterProteins() %>% left_join(annot) %>% dplyr::select(Uniprot, GN, num_pep)
                })
                }
}, selection="single",  options = list(pageLength=10,
                                       lengthChange=F,
                                       filter="top")
)   


# Display the summary statistics of the chosen protein
output$livDisplayProteinText <- renderText({
        if(!is.null(livFilterProteins())){
                filtered_data <- livFilterProteins()
                paste0("After filtering, ", nrow(filtered_data), " proteins remain.")
        }
        
})

# Display the peptides from the chosen PROTEIN
livFilterProteinPeptides <- reactive({
        
        # Get the filtered protein table
        filtered_proteins <- livFilterProteins()
        if(is.null(filtered_proteins)) return(NULL)
        if(nrow(filtered_proteins) == 0) return(NULL)

        # Get the selected protein's Uniprot ID to use it to filter the peptides
        s <- input$livDisplayProteins_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        
        selected_uniprot <- filtered_proteins %>% slice(s) %>% dplyr::select(Uniprot) %>% unlist %>% as.character()
        
        # Read the original, unfildered dataset
        hl <- livRead_hl()
        if(is.null(hl)) return(NULL)
        
        hl <- hl %>% dplyr::filter(Uniprot %in% selected_uniprot,
                                   DP >= input$livFilter_pro_pep_DP[1],
                                   DP <= input$livFilter_pro_pep_DP[2],
                                   
                                   num_k >= input$livFilter_pro_pep_num_k[1],
                                   num_k <= input$livFilter_pro_pep_num_k[2]
                                   
                                   )
        
        hl
        
})

# Display the peptides from the chosen protein on a data table
output$livDisplayProteinPeptides <- DT::renderDataTable({
        
        if(!is.null(livFilterProteinPeptides())) livFilterProteinPeptides() %>% dplyr::select(ID, Uniprot, Peptide, z, DP, ess)
        
}, selection="multiple",  options = list(pageLength=10,
                                       lengthChange=T,
                                       filter="top")
)   


# Display summary statistics from the peptides of the chosen protein
output$livDisplayProteinPeptideText <- renderText({
        if(!is.null(livFilterProteinPeptides())){
                filtered_data <- livFilterProteinPeptides()
                paste0("After filtering, ", nrow(filtered_data), " qualifying peptides exist for the chosen protein.")
        }
        
})



# Fit all selected peptides form selected protein with Steady State
livFitProteinPeptide_SS1 <- reactive({
        
        input$livProteinPeptideGoButton
        
        # Get the dataset
        hl <- livFilterProteinPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl)==0) return(NULL)
        
        # Get the selected row
        s <- isolate(input$livDisplayProteinPeptides_rows_selected) %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- livRead_hldata() %>% filter(ID %in% hl$ID)
        
        # Get the fitted values
        # Get optimized k
        livRefitSinglePeptideFunction(hl=hl, 
                                     hldata=hldata, 
                                     #row_id=row_id, 
                                     fitMethod=input$livFitMethod,
                                     fitModel="SteadyStateOneParameter",
                                     labelMethod=isolate(input$livLabel),
                                     pss=input$livFitting_pss,
                                     kp=input$livFitting_kp)

        })

# Display all selected peptides from the selected protein
output$livDisplayProteinPeptide_SS1 <- renderPlot({
        
        if (is.null(livFitProteinPeptide_SS1())) return(NULL)
        
        
        input$livProteinPeptideGoButton
        
        # Get the dataset
        hl <- livFilterProteinPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl)==0) return(NULL)
        
        # Get the selected row
        s <- isolate(input$livDisplayProteinPeptides_rows_selected) %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- livRead_hldata() %>% filter(ID %in% hl$ID)
        
        
        # Display GG Plot using the fitted values
        livPlotPeptideFunction(hl=hl,
                               hldata=hldata,
                               labelMethod=input$livLabel,
                               plotModel="SteadyState",
                               plotType=input$livPlotType,
                               toggleResidual=input$livToggleResidual,
                               pss=input$livFitting_pss,
                               kp=input$livFitting_kp,
                               fit=livFitProteinPeptide_SS1(),
                               isProtein=T,
                               R2filter=input$livFilter_pro_pep_R2)
        
}, bg="transparent")



# Fit all selected peptides from the selected proteins with Two-Compartment
livFitProteinPeptide_CC1 <- reactive({
        
        input$livProteinPeptideGoButton
        
        # Get the dataset
        hl <- livFilterProteinPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl)==0) return(NULL)
        
        # Get the selected row
        s <- isolate(input$livDisplayProteinPeptides_rows_selected) %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- livRead_hldata() %>% filter(ID %in% hl$ID)
        
        # Get the fitted values
        # Get optimized k
        livRefitSinglePeptideFunction(hl=hl, 
                                             hldata=hldata, 
                                             #row_id=row_id, 
                                             fitMethod=input$livFitMethod,
                                             fitModel="TwoCompartmentOneParameter",
                                             labelMethod=isolate(input$livLabel),
                                             pss=input$livFitting_pss,
                                             kp=input$livFitting_kp)
        
})


# Display all selected peptides from the selected protein
output$livDisplayProteinPeptide_CC1 <- renderPlot({

        
        if (is.null(livFitProteinPeptide_CC1()))  return(NULL)
        
        input$livProteinPeptideGoButton
        
        # Get the dataset
        hl <- livFilterProteinPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl)==0) return(NULL)
        
        # Get the selected row
        s <- isolate(input$livDisplayProteinPeptides_rows_selected) %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the hl-data RIA values then filter by selected ID
        hldata <- livRead_hldata() %>% filter(ID %in% hl$ID)
        
        # Display GG Plot using the fitted values
        livPlotPeptideFunction(hl=hl,
                               hldata=hldata,
                               labelMethod=input$livLabel,
                               plotModel="TwoCompartment",
                               plotType=input$livPlotType,
                               toggleResidual=input$livToggleResidual,
                               pss=input$livFitting_pss,
                               kp=input$livFitting_kp,
                               fit=livFitProteinPeptide_CC1(),
                               isProtein=T,
                               R2filter=input$livFilter_pro_pep_R2)
        
        }, bg="transparent")


output$livDisplayProteinPeptideSS1Text <- renderText({
        if(!is.null(livFitProteinPeptide_SS1())){
                
                # Melting the nested list from the fitted value. 1 is k, 2 is dk, 3 is R2
                df <- melt(livFitProteinPeptide_SS1())
                
                # Filter the datafarme with R2
                results_to_let_in <-  df %>% filter(L2 == 3) %>% filter(value > input$livFilter_pro_pep_R2[1]) %>% filter(value < input$livFilter_pro_pep_R2[2]) %>% dplyr::select(L1) %>% unlist()
                
                # Pick out all the remaining k values
                df <- df %>% filter(L2 == 1) %>% filter(L1 %in% results_to_let_in)
                
                # Calculate geometric mean/log-average (exponent of mean of logs)
                geom.mean = 2^(mean(log2(unlist(df$value)))) 
                
                # Calculate geometric CV if more than 1 peptide
                if (nrow(df) >1) {
                        geom.cv = sqrt(exp(sd(log(df$value))^2) - 1) * 100 
                } else {geom.cv = Inf}
                
                paste0('After filtering, ', nrow(df), ' peptides remain. Geometric mean of k is ', geom.mean %>% round(3), ' and geometric CV is ', geom.cv %>% round(2)," %.")
        }
})

output$livDisplayProteinPeptideCC1Text <- renderText({
        
        if(!is.null(livFitProteinPeptide_CC1())){
                
                # Melting the nested list from the fitted value. 1 is k, 2 is dk, 3 is R2
                df <- melt(livFitProteinPeptide_CC1())
                
                # Filter the datafarme with R2
                results_to_let_in <-  df %>% filter(L2 == 3) %>% filter(value > input$livFilter_pro_pep_R2[1]) %>% filter(value < input$livFilter_pro_pep_R2[2]) %>% dplyr::select(L1) %>% unlist()
                
                # Pick out all the remaining k values
                df <- df %>% filter(L2 == 1) %>% filter(L1 %in% results_to_let_in)
  
                # Calculate geometric mean/log-average (exponent of mean of logs)
                geom.mean = 2^(mean(log2(unlist(df$value)))) 
                
                # Calculate geometric CV if more than 1 peptide
                if (nrow(df) >1) {
                        geom.cv = sqrt(exp(sd(log(df$value))^2) - 1) * 100 
                        } else {geom.cv = Inf}
                
                paste0('After filtering, ', nrow(df), ' peptides remain. Geometric mean of k is ', geom.mean %>% round(3), ' and geometric CV is ', geom.cv %>% round(2)," %.")
                
        }
        
})
