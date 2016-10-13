######
######
###### Refitting Functions
######
######


refitFunction <- function(hl, hldata, fitMethod="Nelder-Mead", fitModel="NS", max_do=500, refit_choice="Graph"){
        
        
                        NS_Refit <- function(var){
                                
                                Refitting_Predicted <- sapply(ds$t, function(x) NSModel(x, var, hl$N[c], hl$kp[c], hl$pss[c], hl$a[c]))
                                Refitting_SS <- sum((ds$A0-Refitting_Predicted)^2)
                                
                                #Use FS rather than A0 for least squares if optimizing for the rescaling function
                                #Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
                                return(Refitting_SS)
                        }
                        
                        # var[1] is k, var[2] is Amax
                        SS_Refit <- function(var){
                                Refitting_Predicted <- sapply(ds$t, function(x) SSModel(x, var[1], hl$a[c], var[2]))
                                Refitting_SS <- sum((ds$A0-Refitting_Predicted)^2)
                                return(Refitting_SS)
                        }
                        
                        # var[1] is ksyn, var[2] is kdeg
                        CC_Refit <- function(var){
                                Refitting_Predicted <- sapply(ds$t, function(x) CCModel(x, var[1], a, var[2], Amax))
                                Refitting_SS <- sum((ds$A0-Refitting_Predicted)^2)
                                return(Refitting_SS)
                        }
                        
                        
                        # To redo - calcaulate FS for refitting, later. 
                        Calculate_FS <- function(x){
                                A0. <- hl$a[c]
                                Ainf. <- hl$a[c]*(1-hl$pss[c])^hl$N[c]
                                FS. <- (x-A0.)/(Ainf.-A0.)
                                return(FS.)
                        }
        
        if (refit_choice == "Graph") return(hl)

        
        # Only take the first 500
        hl <- hl %>% slice (1:max_do)
        
        print(fitModel)
        
        ###
        ### dplyr rowwise version of the for loop (May be marginally faster, but cannot report shiny-style progress bar)
        ###
        # hl <- hl %>% rowwise() %>% do ({    
        #         # Subsetting the hl_data to get only the particular ID being considered
        #         ds <- filter(hldata, ID == .$ID)
        #         # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
        #         Optimize <- optim(.$k, function(var) sum((ds$A0 -sapply(ds$t, function(y) NSModel(y, var, .$N, .$kp, .$pss, .$a)))^2), method = fitMethod) #optim(start value, fxn) # Use optim() for Nelder-Mead
        #         Refitted_SS <- Optimize$value
        #         Refitted_R2 <- 1- (Refitted_SS/(sum((ds$A0 - mean(ds$A0))^2)))
        #         dk <- sapply(ds$t, function(x) NSModel_dk(x, Optimize$par, .$N, .$kp, .$pss, .$a, .$SS)) %>% abs() %>% min()
        #         data.frame(., Optimize$par, dk, Refitted_R2, Refitted_SS)
        # })   
        
        # LOOP VERSION
        withProgress(message = 'Refitting Data...', value = 0, {
                for (c in 1:nrow(hl)) {
                        incProgress(1/nrow(hl))
                        
                        ds <- filter(hldata, ID == hl$ID[c])
                        # ## Calculate Fractional Synthesis
                        # FS <- sapply(ds$A0, Calculate_FS) # Calculate FS only if you are interested in the rescaling (combining peptides for fitting)
                        # ds <- cbind(ds,FS)
                        
                        # This is the optimization function block that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K.
                        
                        # The SS model
                        if(fitModel == "SS"){
                                Optimize <- optim(c(0.25,0.05), SS_Refit, method = fitMethod, control = list(maxit=100))
                                
                                k <- Optimize$par[1]
                                SS <- Optimize$value
                                SE <- (Optimize$value/(hl$DP[c]-1))^0.5
                                dk <- sapply(ds$t, function(x) SSModel_dk(x, Optimize$par[1], hl$a[c], Optimize$par[2], SE)) %>% abs() %>% min()
                                Amax <- Optimize$par[2] # This is actually the optimized Amax. Just putting it here for now.
                                
                        }
                        
                        # THe NS model
                        if(fitModel == "NS"){
                                if(fitMethod == "Brent"){
                                        Optimize <- optim(0.25, NS_Refit, method = "Brent", lower = 0, upper = 4, control = list(maxit=15))
                                } else {
                                        Optimize <- optim(0.25, NS_Refit, method = fitMethod, control = list(maxit=20)) 
                                }
                                
                                k <- Optimize$par
                                SS <- Optimize$value
                                SE <- (Optimize$value/(hl$DP[c]-1))^0.5
                                
                                # Calculate dk (fitting error) using analytical solution. Note: there is a dk for every t - we calculate it from the time point where
                                # k would be most sensitive to A.

                                dk <- sapply(ds$t, function(x) NSModel_dk(x, Optimize$par, hl$N[c], hl$kp[c], hl$pss[c], hl$a[c], SE)) %>% abs() %>% min()
                                Amax <- hl$a[c] * (1 - hl$pss[c])^hl$N[c]
                        }
                        
                        if(fitModel == "CC"){
                                Amax <- hl$a[c] * (1 - hl$pss[c])^hl$N[c]
                                
                                Optimize <- optim(c(0.25,1), CC_Refit, method = fitMethod, control = list(maxit=100))
                                
                                k <- Optimize$par[1]
                                hl$k2[c] <- Optimize$par[2]
                                SS <- Optimize$value
                                SE <- (Optimize$value/(hl$DP[c]-1))^0.5
                                dk <- sapply(ds$t, function(x) CCModel_dk(x, Optimize$par[1], a, Optimize$par[2], Amax, SE)) %>% abs() %>% min()
                                ksyn <- Optimize$par[2] # This is actually the optimized Amax. Just putting it here for now.
                                
                        }
                        
                        R2 <- 1- (SS/(sum((ds$A0 - mean(ds$A0))^2)))
                        
                        hl$Amax[c] <- Amax %>% round(3)
                        hl$k[c] <- k %>% round(4)
                        hl$dk[c] <- dk %>% round(4)
                        hl$R2[c] <- R2 %>% round(3)
                        hl$SS[c] <- SS %>% round(4)
                        hl$SE[c] <- SE %>% round(4)
                        
                }
        })
        
        glimpse(hl)        
        
        hl
}



######
######
###### ProTurn Refit SERVER reactive functions
######
######

readhl <- reactive({
        hlout_File <- input$file1
        
        if (is.null(input$file1)){return(NULL)}
        
        hl <- read.table(hlout_File$datapath, header=T, fill=T, as.is=T, sep="\t") %>% arrange(Uniprot)
        hl <- hl %>% group_by() %>% do({
                SE = (.$SS/(.$DP-1))^0.5
                concat = paste0(.$Peptide, "[", .$z, "+]") 
                oldk = .$k
                oldSS = .$SS
                olddk = .$dk
                Amax = 0
                k2 = .$k
                len = nchar(gsub("\\(.*\\)","", .$Peptide))
                data.frame(., SE, concat, oldk, Amax, oldSS, olddk, k2, len)
        })
        glimpse(hl)          
        hl %>% ungroup()
})   

readhldata <- reactive({
        hldataout_File <- input$file2
        if (is.null(input$file2)){return(NULL)}
        
        hldata <- read.table(hldataout_File$datapath, header=T, fill=T, as.is=T, sep="\t")
        hldata[,1:3]
})   


refitProturn <- reactive({
        
        # Dependency on the Refit button
        input$startRefitButton
        
        if (is.null(isolate(readhl()))){return(NULL)}
        if (is.null(isolate(readhldata()))){return(NULL)}
        
        hl <-    isolate(readhl())
        hldata <- isolate(readhldata())
        
        if(!nrow(hl) > 0) {return(NULL)}                
        refitFunction(hl, hldata, isolate(input$fitMethod),isolate(input$fitModel), 2000, isolate(input$refit_choice))
        
})  

output$refitStatus <- renderText({
        refitted_data <- refitProturn()
        if(is.null(refitted_data)) return(NULL)
        paste(nrow(refitted_data), " peptides found. Proceed to the next step.")
        
})

# Filter the refitting output by R2, or DP, for display     
filterRefit <- reactive({
        
        refitted_data <- refitProturn()
        if(is.null(refitted_data))  {return(NULL)}
        
        refitted_data %>% dplyr::filter(R2 >= input$result_R2filter[1], 
                                        R2 <= input$result_R2filter[2],
                                        len >= input$result_lenfilter[1],
                                        len <= input$result_lenfilter[2],
                                        DP >= input$result_DPfilter[1],
                                        DP <= input$result_DPfilter[2])
        
})


# Summarize the output of refitting into a data table (Protein level selection table)
getRefitProteins <- reactive({
        if(!is.null(filterRefit()))  {
                refitted_data <- filterRefit() 
                refitted_data <- refitted_data %>% ungroup() %>% group_by(Uniprot) %>% 
                        summarize(k.median = median(k) %>% round(4), 
                                  k.log.mean = 2^(mean(log2(k))) %>% round(4), 
                                  k.mad = mad(k) %>% round(4), 
                                  geom.cv = sqrt(exp((sd(log2(k))*log(2))^2)-1) %>% round(4) , 
                                  n = n()
                                  ) 
                
        }
})

# Display the protein table
output$refitProteinSummary <- DT::renderDataTable({
        if(!is.null(getRefitProteins())){
                getRefitProteins()
        }
}, selection="single",  options = list(lengthMenu = c(10,20), 
                                       pageLength = 10,
                                       filter="top",
                                       columnDefs = list(list(width = '200px', targets= "_all"))
))     

# Display the protein table's summary statistics
output$refitProteinText <- renderText({
        if(!is.null(getRefitProteins())){
                refitted_data <- getRefitProteins()
                paste0("After filtering, there remains ", nrow(refitted_data), " proteins and ", sum(refitted_data$n, na.rm=T), " peptides. The median CV is ", median(refitted_data$geom.cv, na.rm=T)*100,"%.")
        }
        
})

# Subset the filtered result with the selected Protein
getRefitPeptides <- reactive({
        if(!is.null(filterRefit()))  {
                i <- input$refitProteinSummary_rows_selected %>% as.integer()
                if(!length(i)) return(NULL)
                selected_uniprot <- getRefitProteins() %>% slice(i) %>% dplyr::select(Uniprot) %>% unlist()
                filterRefit() %>% filter(Uniprot == selected_uniprot)
        }
})

# Summarize the output of refitting into a data table (Peptide level selection table)
output$refitPeptideSummary <- DT::renderDataTable({
        if(is.null(getRefitPeptides())) return(NULL)
        
        getRefitPeptides() %>% dplyr::select(ID, Uniprot, Seq = concat, DP, k, oldk, dk, olddk, oldSS, SS, SE, R2) #%>% 
        
}, selection="single",  options = list(lengthMenu = c(10,20), 
                                       pageLength = 10,
                                       columnDefs = list(list(width = '200px', targets= "_all"))
))    # filter="top"

# Display individual peptides from the refitting results
output$displayPeptide <- renderPlot({
        
        # Dependency on the Refit button
        input$startRefitButton
        
        refitted <- getRefitPeptides()
        
        if(is.null(refitted)) return(NULL)
        if(nrow(refitted) == 0) return(NULL)
        
        s <- input$refitPeptideSummary_rows_selected %>% as.integer()
        
        if(!length(s)) return(NULL)
        print(s)
        
        refitted_data <- refitted %>% slice(s)
        
        hldata <- isolate(readhldata())
        hldata <- hldata %>% filter(ID %in% refitted_data$ID)
        
        
        k <- refitted_data$k
        k2 <- refitted_data$k2
        dk <- refitted_data$dk
        N <- refitted_data$N
        kp <- refitted_data$kp
        pss <- refitted_data$pss
        a <- refitted_data$a
        Amax <- refitted_data$Amax
        
        annot <- paste(refitted_data$concat, as.character(refitted_data$k), as.character(refitted_data$R2))
        
        
        g <- ggplot(data = hldata, aes(x = t, y = A0))
        g <- g + geom_point()
        
        # Plot out the curves for the NS model
        if(input$fitModel == "NS"){
                g <- g + stat_function(fun = function(x) NSModel(x, k, N, kp, pss, a))
                g <- g + stat_function(fun = function(x) NSModel(x, k+dk, N, kp, pss, a), color = "red")
                g <- g + stat_function(fun = function(x) NSModel(x, k^2/(k+dk), N, kp, pss, a), color="red")
        }
        
        # Plot out the curves for the SS model
        if(input$fitModel == "SS"){
                g <- g + stat_function(fun = function(x) SSModel(x, k, a, Amax))
                g <- g + stat_function(fun = function(x) SSModel(x, k+dk, a, Amax), color = "red")
                g <- g + stat_function(fun = function(x) SSModel(x, k^2/(k+dk), a, Amax), color="red")
        }
        
        # Plot out the curves for the SS model
        if(input$fitModel == "CC"){
                g <- g + stat_function(fun = function(x) CCModel(x, k, k2, a, Amax))
                g <- g + stat_function(fun = function(x) CCModel(x, k+dk, k2, a, Amax), color = "red")
                g <- g + stat_function(fun = function(x) CCModel(x, k^2/(k+dk), k2, a, Amax), color="red")
        }
        
        
        #Find x limit
        gx <- ggplot_build(g)$panel$ranges[[1]]$x.range[2]
        gy <- ggplot_build(g)$panel$ranges[[1]]$y.range[2]
        g <- g + xlim(0, gx)
        g <- g + annotate("text", label = annot , x = gx/2, y = gy, size = 6)
        g <- g + ggtitle(annot)
        
        #g <- ggplotly(g)
        g
})

# Display individual peptides from the refitting results
output$displayResidual <- renderPlot({
        
        # Dependency on the Refit button
        input$startRefitButton
        
        refitted <- getRefitPeptides()
        
        if(is.null(refitted)) return(NULL)
        if(nrow(refitted) == 0) return(NULL)
        
        s <- input$refitPeptideSummary_rows_selected %>% as.integer()
        
        if(!length(s)) return(NULL)
        
        refitted_data <- refitted %>% slice(s)
        
        #filter(refitted_data, Uniprot == as.character(input$pickRefitProtein)) %>% filter(concat == as.character(input$pickRefitPeptide))
        #if(is.null(refitted_data)) return(NULL)
        
        hldata <- isolate(readhldata())
        hldata <- hldata %>% filter(ID %in% refitted_data$ID)
        
        
        k <- refitted_data$k
        k2 <- refitted_data$k2
        dk <- refitted_data$dk
        N <- refitted_data$N
        kp <- refitted_data$kp
        pss <- refitted_data$pss
        a <- refitted_data$a
        Amax <- refitted_data$Amax
        
        if(input$fitModel == "NS"){
                hldata <- hldata %>% rowwise() %>% do({
                        residual <- .$A0 - NSModel(.$t, k, N, kp, pss, a)
                        data.frame(., residual)
                })
        }
        
        if(input$fitModel == "SS"){
                hldata <- hldata %>% rowwise() %>% do({
                        residual <- .$A0 - SSModel(.$t, k, a, Amax)
                        data.frame(., residual)                  
                })
        }
        
        if(input$fitModel == "CC"){
                hldata <- hldata %>% rowwise() %>% do({
                        residual <- .$A0 - CCModel(.$t, k, k2, a, Amax)
                        data.frame(., residual)                  
                })
        }
        
        
        g <- ggplot(data = hldata, aes(x = t, y = residual))
        g <- g + geom_point()
        g <- g + geom_hline(yintercept=0, col="black")
        g <- g + geom_hline(yintercept= c(abs(max(hldata$A0)-min(hldata$A0))/2,-abs(max(hldata$A0)-min(hldata$A0))/2), col="red")
        g <- g + ggtitle("Residual Plot")
        g
})


# Display the protein table on the refit Protein page.
# Shiny doesn't seem to allow the same output DT table to be displayed twice, so we have to make a duplicate

# Filter the refitting output by R2, or DP, for display     
filterRefit2 <- reactive({
        
        refitted_data <- refitProturn()
        if(is.null(refitted_data))  {return(NULL)}
        
        refitted_data %>% dplyr::filter(R2 >= input$result_R2filter2[1], 
                                        R2 <= input$result_R2filter2[2],
                                        len >= input$result_lenfilter[1],
                                        len <= input$result_lenfilter[2],
                                        DP >= input$result_DPfilter2[1],
                                        DP <= input$result_DPfilter2[2])
        
})


# Summarize the output of refitting into a data table (Protein level selection table)
getRefitProteins2 <- reactive({
        if(is.null(filterRefit2()))  return(NULL)
                refitted_data <- filterRefit2() 
                refitted_data <- refitted_data %>% ungroup() %>% group_by(Uniprot) %>% 
                        summarize(k.median = median(k) %>% round(4), 
                                  k.log.mean = 2^(mean(log2(k))) %>% round(4), 
                                  k.mad = mad(k) %>% round(4), 
                                  geom.cv = sqrt(exp((sd(log2(k))*log(2))^2)-1) %>% round(4), 
                                  n = n()
                        ) 
                
})

output$refitProteinSummary2 <- DT::renderDataTable({
        if(!is.null(getRefitProteins2())){
                getRefitProteins2()
        }
        
}, selection="single",  options = list(lengthMenu = c(10,20), 
                                       pageLength = 10,
                                       filter="top",
                                       columnDefs = list(list(width = '200px', targets= "_all"))
))    

# Display the protein table's summary statistics
output$refitProteinText2 <- renderText({
        if(!is.null(getRefitProteins2())){
                refitted_data <- getRefitProteins2()
                paste0("After filtering, there remains ", nrow(refitted_data), " proteins and ", sum(refitted_data$n, na.rm=T), " peptides. The median CV is ", median(refitted_data$geom.cv, na.rm=T)*100,"%.")
        }
        
})

# Subset the filtered result with the selected Protein
getRefitPeptides2 <- reactive({
        if(!is.null(filterRefit2()))  {
                i <- input$refitProteinSummary2_rows_selected %>% as.integer()
                if(!length(i)) return(NULL)
                selected_uniprot <- getRefitProteins2() %>% slice(i) %>% dplyr::select(Uniprot) %>% unlist()
                filterRefit2() %>% filter(Uniprot == selected_uniprot)
        }
})

# Summarize the output of refitting into a data table (Peptide level selection table)
output$refitPeptideSummary2 <- DT::renderDataTable({
        if(is.null(getRefitPeptides2())) return(NULL)
        
        getRefitPeptides2() %>% dplyr::select(ID, Uniprot, Seq = concat, DP, k, oldk, dk, olddk, oldSS, SS, SE, R2) #%>% 
        
}, selection="multiple",  options = list(lengthMenu = c(10,20), 
                                       pageLength = 10,
                                       columnDefs = list(list(width = '200px', targets= "_all"))
))    # filter="top"
