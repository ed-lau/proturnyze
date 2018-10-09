######
######
###### Helper Functions that are Shared Across Multiple Software Functions
######
######

# Take in a gene list and a background list (mouse uniprot IDs), then perform GO enrichment analysis
enrichFunction <- function(genelist_go, background_go, pfilter = 0.05){

    withProgress(message = 'Calculating Significance...', value = 0, {
    # This is weird. Should this count be the genelist/background length, or only those that are associated with at least one GO?
    #genelist_count <- nrow(genelist)
    #background_count <- nrow(background)
        genelist_count <- genelist_go %>% distinct(uniprot) %>% nrow()
        background_count <- background_go %>% distinct(uniprot) %>% nrow()
    
    
        genelist_go_summary <- genelist_go %>% group_by(go_id) %>% summarize(proteins = length(unique(uniprot)), go_name = first(go_name), aspect=first(aspect)) %>% arrange(-proteins)
    
    
        for (i in 1:nrow(genelist_go_summary)){
          ls_go <- genelist_go_summary$proteins[i]
          ls_not_go <- genelist_count - genelist_go_summary$proteins[i]
          bg_go <- sum(background_go$go_id == genelist_go_summary$go_id[i])
          bg_not_go <- background_count - bg_go
      
          test <- fisher.test(matrix(c(ls_go,ls_not_go,bg_go,bg_not_go), nrow=2, ncol=2))
          genelist_go_summary$p[i] <- test$p.value
          genelist_go_summary$est[i] <- test$estimate
          genelist_go_summary$conf1[i] <- test$conf.int[1]
          genelist_go_summary$conf2[i] <- test$conf.int[2]
          incProgress(0.5/nrow(genelist_go_summary)) 
          incProgress(1/i)
          }
    
    genelist_go_summary$BH <- p.adjust(genelist_go_summary$p, method="BH")
    output <- filter(genelist_go_summary, BH <= pfilter, proteins >= 2) %>% dplyr::arrange(BH)
    
    # Get the names and aspect of the enriched GO terms through web service (alt: cache them too?)
    if (!nrow(output) > 0){return(NULL)}
      
    # Retrieve the GO names things automatically
#      for (i in 1:nrow(output)){
#          go_details <- read.table(paste("http://www.ebi.ac.uk/QuickGO/GAnnotation?format=tsv&limit=1&goid=",output$go_id[i],"&col=goID,goName,aspect",sep=""), fill=T, header=T, as.is=T, sep="\t")
#          output$go_name[i] <- go_details$GO.Name[1]
#          output$aspect[i] <- go_details$Aspect[1]
#          incProgress(1/i)
#        }
    
  })   
  return(output)
}

######
######
###### Enrich helper
######
######

getTax <- function(gene){
  tax <- read.table(paste0("http://www.uniprot.org/uniprot/?query=uniprotkb:", gene, "&format=tab&limit=1&columns=organism-id"), fill=T, header=T, as.is=T)
  tax <- tax[1,1]
  return(tax)
}

getTaxUniprot <-function(tax){
  taxUniprot <- read.table(paste0("http://www.uniprot.org/uniprot/?query=taxonomy:", tax ,"+reviewed:yes&format=tab&columns=id,database(OMA)"), header=T, fill=T, as.is=T) %>% dplyr::select(tax_uniprot = Entry, OMA = Cross.reference)
  taxUniprot <- taxUniprot %>% filter(OMA != "") %>% distinct() #%>% filter(tax_uniprot %in% genelist$uniprot) %>% left_join(allUniprot) %>% filter(uniprot != "NA")
  return(taxUniprot)
}


getUniprotGN <- function(Uniprot){ 
        taxUniprot <- read.table(paste0("http://www.uniprot.org/uniprot/?query=accession:",Uniprot,"&format=tab&columns=id,genes"), header=T, fill=T, as.is=T, skip=1, sep="\t")
        colnames(taxUniprot) <- c("Uniprot", "GN")
        taxUniprot$GN <- gsub(" .*","", taxUniprot$GN)
        return(taxUniprot$GN[1])
}