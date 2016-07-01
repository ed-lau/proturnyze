###
### These are wrapper functions that separate parts of the UIs into separate pages for tidiness.
###


chat_page <- function(){
        
        
        tabPanel("Chat",
                 tagline(),
                 sidebarLayout(
                         sidebarPanel(
                                 h3("Chatroom"),
                                 tags$hr(),
                                 p("Ask questions and share your findings here."),
                                 br()
                                 ), 
                         mainPanel(
                                 includeCSS("shinychat.css"),
                                 
                                 # And custom JavaScript -- just to send a message when a user hits "enter"
                                 # and automatically scroll the chat window for us. Totally optional.
                                 includeScript("sendOnEnter.js"),
                                 
                                 div(
                                         # Setup custom Bootstrap elements here to define a new layout
                                         class = "container-fluid", 
                                         div(class = "row-fluid",
                                             # Set the page title
                                             tags$head(tags$title("Chatroom"))
                                            
                                         ),
                                         # The main panel
                                         div(
                                                 class = "row-fluid", 
                                                 mainPanel(
                                                         # Create a spot for a dynamic UI containing the chat contents.
                                                         uiOutput("chat"),
                                                         
                                                         # Create the bottom bar to allow users to chat.
                                                         fluidRow(
                                                                 div(class="span10",
                                                                     textInput("entry", "")
                                                                 ),
                                                                 div(class="span2 center",
                                                                     actionButton("send", "Send")
                                                                 )
                                                         )
                                                 ),
                                                 # The right sidebar
                                                 sidebarPanel(
                                                         # Let the user define his/her own Name
                                                         textInput("user", "Your Name:", value=""),
                                                         tags$hr(),
                                                         h5("Connected Users"),
                                                         # Create a spot for a dynamic UI containing the list of users.
                                                         uiOutput("userList"),
                                                         tags$hr(),
                                                         helpText("Available Helpers: ","Maggie Lam", "Brian Bleakley")
                                                 )
                                         )
                                 )
                                 )
                                 )
                 )
        
        
}
