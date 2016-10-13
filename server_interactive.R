######
######
###### Interactive Refitting Functions
######
######


livRefitSinglePeptideFunction <- function(hl, 
                                          hldata, 
                                          row_id, 
                                          fitMethod="Nelder-Mead", 
                                          fitModel="SteadyStateOneParameter", 
                                          labelMethod="aa",
                                          pss,
                                          kp){
        
        
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
        
        # Get the correct kp, Amax, and a value
        if(labelMethod == "aa"){
                Amax <- pss
                a <- 0
                kp <- kp
        } else {
                Amax <- pss # Placeholder - calcualte Amax using heavy water formula later
                a <- 0 # Placeholder - calcualte a using heavy water formula later
                kp <- kp
        }
        
        
        # Slicing hl and hldata.out by the ID
        hl <- hl %>% filter(ID == row_id) %>% slice(1)
        ds <- hldata %>% filter(ID == row_id)
        
        # LOOP VERSION
        withProgress(message = 'Refitting Data...', value = 0.5, {
               
                if(fitModel == "SteadyStateOneParameter"){

                        if(fitMethod == "Brent"){
                                Optimize <- optim(0.249, SS1_Refit, method = "Brent", lower = 0, upper = 4, control = list(maxit=15))
                        } else {
                                Optimize <- optim(0.249, SS1_Refit, method = fitMethod, control = list(maxit=100))
                        }
                        
                        k <- Optimize$par
                        SS <- Optimize$value
                        SE <- (Optimize$value/(hl$DP-1))^0.5
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
                        SE <- (Optimize$value/(hl$DP-1))^0.5
                        dk <- sapply(ds$t, function(x) CCModel_dk(x, Optimize$par, kp, a, Amax, SE)) %>% abs() %>% min()
                        
                }
                
                R2 <- 1- (SS/(sum((ds$A0 - mean(ds$A0))^2)))

        })
        
        list(round(k, 3),
             round(dk, 3),
             round(R2, 3))
}



######
######
###### Interactive Refit SERVER reactive functions
######
######

livRead_hl <- reactive ({
       
        organ <- input$livOrgan 
        label_type <- input$livLabel 
        # 
        file_path <- paste0("data/", label_type, "_", organ, ".txt")
        
        withProgress(message = 'Reading Data...', value = 0.2, {
                
        hl <- read.table(file_path, header=T, fill=T, as.is=T, sep="\t") 
        
        incProgress(0.5)
        
        hl <- hl %>% group_by(ID) %>% do({
                
                concat = paste0(.$Peptide, "[", .$z, "+]") 
                
                #N = calcLabelingSite(.$Peptide)[1]
                #a = calc_thr_a(.$Peptide)
                #Amax = a * (1 - .$pss)^N
                
                ess = calc_essentiality(.$Peptide)
                
                len = nchar(gsub("\\(.*\\)","", .$Peptide))
                num_k = stringr::str_count(.$Peptide, "K")
                
                data.frame(., concat, len, num_k, ess)
                
        }) %>% ungroup()
        
        })
        hl 
        
})  
        

livRead_hldata <- reactive({

        organ <- input$livOrgan
        label_type <- input$livLabel
        
        file_path <- paste0("data/", label_type, "_", organ, "_data.txt")
        
        hldata <- read.table(file_path, header=T, fill=T, as.is=T, sep="\t")

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
        
        if(!is.null(livFilterPeptides())) livFilterPeptides() %>% dplyr::select(ID, Uniprot, Peptide, z, DP)
        
        }, selection="single",  options = list(pageLength=10,
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

# Refit the selected peptide using the SteadyStateOneParameter model
livRefitSelectedPeptide_SS <- reactive({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl) == 0) return(NULL)
        
        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        row_id <- hl %>% slice(s) %>% dplyr::select(ID) %>% unlist %>% as.numeric()
        
        
       livRefitSinglePeptideFunction(hl=isolate(livFilterPeptides()), 
                               hldata=isolate(livRead_hldata()), 
                               row_id=row_id, 
                               fitMethod=isolate(input$livFitMethod),
                               fitModel="SteadyStateOneParameter",
                               labelMethod=isolate(input$livLabel),
                               pss=input$livFitting_pss,
                               kp=input$livFitting_kp)

        # Return a list containing (k, dk, R2) 
 
})  

livRefitSelectedPeptide_CC <- reactive({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl) == 0) return(NULL)
        
        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        row_id <- hl %>% slice(s) %>% dplyr::select(ID) %>% unlist %>% as.numeric()
        
        livRefitSinglePeptideFunction(hl=isolate(livFilterPeptides()), 
                                   hldata=isolate(livRead_hldata()), 
                                   row_id=row_id, 
                                   fitMethod=isolate(input$livFitMethod),
                                   fitModel="TwoCompartmentOneParameter",
                                   labelMethod=isolate(input$livLabel),
                                   pss=input$livFitting_pss,
                                   kp=input$livFitting_kp)
        
        # Return a list containing (k, dk, R2) 
})  


# Display individual peptides from the refitting results
output$livDisplayPeptide_SS1 <- renderPlot({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl) == 0) return(NULL)
        
        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the correct kp, Amax, and a value
        if(isolate(input$livLabel) == "aa"){
                Amax <- input$livFitting_pss
                a <- 0
        } else {
                Amax <- input$livFitting_pss # Placeholder - calcualte Amax using heavy water formula later
                a <- 0 # Placeholder - calcualte a using heavy water formula later
        }
        
        # Get the hl-data RIA values
        hldata <- livRead_hldata()
        hldata <- hldata %>% filter(ID %in% hl$ID)
        
        # Get the fitted values
        fit <- livRefitSelectedPeptide_SS()
        k <- fit[1] %>% unlist()
        dk <- fit[2] %>% unlist()
        R2 <- fit[3] %>% unlist()
        kp <- input$livFitting_kp
        annot <- paste0("First-Order Kinetics Model\n", "Peptide: ", hl$concat, "\n k: ", k, " R2: ", R2)
        
        hldata <- hldata %>% rowwise() %>% do({
                predicted <- SSModel(.$t, k, a, Amax)
                residual <- .$A0 - predicted
                data.frame(., predicted, residual)                  
        })
        
        # Plot out the curves for the SteadyState model
        g <- ggplot(data = hldata, aes(x = t, y = A0))
        g <- g + geom_point() + ggtitle(annot)
        g <- g + stat_function(fun = function(x) SSModel(x, k, a, Amax))
        g <- g + stat_function(fun = function(x) SSModel(x, k+dk, a, Amax), color = "red")
        g <- g + stat_function(fun = function(x) SSModel(x, k^2/(k+dk), a, Amax), color="red")
        
        # Dropping a line from the data point to the curve
        #g <- g + geom_segment(aes(x=t, xend=t, y=predicted, yend=A0))
        
        #Find x limit
        gx <- ggplot_build(g)$panel$ranges[[1]]$x.range[2]
        gy <- ggplot_build(g)$panel$ranges[[1]]$y.range[2]
        g <- g + xlim(0, gx)
        
        
        # Plot residual plot instead of best-fit curve
        if (input$livPlotType == 2){

                
                g <- ggplot(data = hldata, aes(x = t, y = residual))
                g <- g + geom_point() + ggtitle("First-Order Kientics Model - Residual Plot")
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
}, bg="transparent")

# Display individual peptides from the refitting results
output$livDisplayPeptide_CC1 <- renderPlot({
        
        # Get the dataset
        hl <- livFilterPeptides()
        if(is.null(hl)) return(NULL)
        if(nrow(hl) == 0) return(NULL)
        
        # Get the selected row
        s <- input$livDisplayPeptides_rows_selected %>% as.integer()
        if(!length(s)) return(NULL)
        hl <- hl %>% slice(s)
        
        # Get the correct kp, Amax, and a value
        if(isolate(input$livLabel) == "aa"){
                Amax <- input$livFitting_pss
                a <- 0
        } else {
                Amax <- input$livFitting_pss # Placeholder - calcualte Amax using heavy water formula later
                a <- 0 # Placeholder - calcualte a using heavy water formula later
        }
        
        # Get the hl-data RIA values
        hldata <- livRead_hldata()
        hldata <- hldata %>% filter(ID %in% hl$ID)
        
        # Get the fitted values
        fit <- livRefitSelectedPeptide_CC()
        k <- fit[1] %>% unlist()
        dk <- fit[2] %>% unlist()
        R2 <- fit[3] %>% unlist()
        kp <- input$livFitting_kp
        annot <- paste0("Two-Compartment Model\n", "Peptide: ", hl$concat, "\n k: ", k, " R2: ", R2)
        
        # Plot out the curves for the TwoComartment model
        g <- ggplot(data = hldata, aes(x = t, y = A0))
        g <- g + geom_point() + ggtitle(annot)
        g <- g + stat_function(fun = function(x) CCModel(x, k, kp, a, Amax))
        g <- g + stat_function(fun = function(x) CCModel(x, k+dk, kp, a, Amax), color = "red")
        g <- g + stat_function(fun = function(x) CCModel(x, k^2/(k+dk), kp, a, Amax), color="red")

        #Find x limit
        gx <- ggplot_build(g)$panel$ranges[[1]]$x.range[2]
        gy <- ggplot_build(g)$panel$ranges[[1]]$y.range[2]
        g <- g + xlim(0, gx)

        # Plot residual plot instead of best-fit curve
        if(input$livPlotType == 2)
        {
                hldata <- hldata %>% rowwise() %>% do({
                        residual <- .$A0 - CCModel(.$t, k, kp, a, Amax)
                        data.frame(., residual)
                })
                
                g <- ggplot(data = hldata, aes(x = t, y = residual))
                g <- g + geom_point() + ggtitle("Two-Compartment Model - Residual Plot") 
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
}, bg="transparent")

