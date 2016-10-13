###
require("dplyr"); require("ggplot2"); require("tidyr") 
require("reshape2"); require("heatmap3"); require("gplots"); require(caret)
###


annot <- read.table("annot.tab", fill=T, as.is=T, header=T, sep="\t", quote="")

# Remove all the gene names except teh first one
annot$GN <- gsub(" .*", "", annot$GN)

all_protein_k <- read.table("all_protein_k.txt", fill=T, as.is=T, header=T, sep = "\t", quote="")
names(all_protein_k) <- c("Uniprot", "Strain", "Factor", "median.k", "mad.k")

all_protein_k <- left_join(all_protein_k, annot)
all_protein_k <- all_protein_k %>% dplyr::select(Uniprot,GN,PN,Strain,Factor,median.k, mad.k)


write.table(all_protein_k, "all_protein_k.txt", quote=F, row.names=F, sep="\t")

trinity_flat$Strain <- gsub("cej", "CE/J", trinity_flat$Strain)
trinity_flat$Strain <- gsub("fvb", "FVB/NJ", trinity_flat$Strain)
trinity_flat$Strain <- gsub("dba", "DBA/2J", trinity_flat$Strain)
trinity_flat$Strain <- gsub("c57", "C57BL/6J", trinity_flat$Strain)
trinity_flat$Strain <- gsub("balbc", "BALB/cJ", trinity_flat$Strain)
trinity_flat$Strain <- gsub("aj", "A/J", trinity_flat$Strain)
trinity_flat$Factor <- gsub("ctrl", "Normal", trinity_flat$Factor)
trinity_flat$Factor <- gsub("iso", "Hypertrophy", trinity_flat$Factor)

# Add the cellular compartments to trinity
trinity_flat <- read.table("trinity_flat.txt", fill=T, as.is=T, header=T, sep = "\t", quote="")

all_mouse_go <- read.table("data/all_mouse_go.txt", fill=T, as.is=T, header=T, sep = "\t", quote="")
# Get only terms from 10 major organelles (Sci Data manuscript)
all_mouse_cc <- all_mouse_go %>% filter(aspect == "Component")
all_mouse_cc <- all_mouse_cc %>% filter(go_id %in% c("GO:0005739", "GO:0005829", "GO:0005634", "GO:0005886",
                                                     "GO:0005615", "GO:0005783", "GO:0005794", "GO:0005768",
                                                     "GO:0005777", "GO:0005764"))
colnames(all_mouse_cc)[1] <- "Uniprot"
trinity_flat_cc <- trinity_flat %>% left_join(all_mouse_cc)
# Raplace no compartments (NA) proteins with "others" label
trinity_flat_cc$go_name <- replace(trinity_flat_cc$go_name, is.na(trinity_flat_cc$go_name), "others")
trinity_flat_cc <- trinity_flat_cc %>% dplyr::select(-aspect, -go_id)

write.table(trinity_flat_cc, "trinity_flat_cc.txt", quote=F, row.names=F, sep="\t")
