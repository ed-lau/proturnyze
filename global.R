
######
######
###### Globals
######
######

###
### Big datasets that are shared across sessions
###

proteinValue <- read.table("data/trinity_flat_cc.txt", header=T, sep="\t", quote="", fill=T)
demoData <- read.table("data/trinity.txt", header=T, sep="\t", quote="", fill=T)


###
### Chat box
###


vars <- reactiveValues(chat=NULL, users=NULL) # Globally define a place where all users can share some reactive data.

if (file.exists("chat.Rds")){
        vars$chat <- readRDS("chat.Rds")
} else {
        vars$chat <- "Welcome to the live help chat!"
}
linePrefix <- function(){
        if (is.null(isolate(vars$chat))){
                return("")
        }
        return("<br />")
}

