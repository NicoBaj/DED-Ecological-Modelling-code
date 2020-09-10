# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# A DED simulator
# 
# Based on the paper "Assessing the risk of COVID-19 importation and the effect of quarantine"
# by Julien Arino, Nicolas Bajeux, Stephanie Portet and James Watmough
# https://www.medrxiv.org/content/10.1101/2020.08.12.20173658v1

library(shiny)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Spread of Dutch Elm Disease in Winnipeg"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("Neighbourhood", 
                  label = "Choose a neighbourhood",
                  choices = c("North River Heights (NRH)", 
                              "Mixed Pulberry-Victoria Crescent-Crescent Park (PVC)"),
                  selected = "North River Heights (NRH)"),
      
      selectInput("IC", 
                  label = "Choose an initial condition",
                  choices = c("Cluster", 
                              "Two clusters",
                              "Random"),
                  selected = "Cluster"),
      
      sliderInput("R_B",
                  "Maximum beetle dispersal distance (m):",
                  min = 20,
                  max = 380,
                  step = 40,
                  value = 100),
      sliderInput("s_dt",
                  "Proportion of beetles surviving one week:",
                  min = 0.97,
                  max = 0.99,
                  step = 0.01,
                  value = 0.99),
      sliderInput("p_i",
                  "Proba. for a beetle to succesfully infect a tree:",
                  min = 0.01,
                  max = 0.05,
                  value = 0.03,
                  step = 0.01),
      sliderInput("p_r",
                  "Proba. for an inf. tree to succesfully infect another one:",
                  min = 0,
                  max = 1,
                  value = 0.1),
      actionButton("Generate", "Run one simulation")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Neighbourhood",
                           textOutput("info_ngh"),
                           imageOutput("plot_ngh",height=10)),
                  tabPanel("Maps", 
                           plotOutput("plot_maps", width = "600px", height = "450px")),
                  # tabPanel("Proportions",
                  #          plotOutput("plot_proportions", width = "800px", height = "600px"))
                  #tabPanel("Effect of quarantine", verbatimTextOutput("summary"))
                  #tabPanel("Table", tableOutput("table"))
      )
    )
  )
)

# Define server logic required to draw the outputs of one simulation
server <- function(input, output) {
  
  # Reactive expression to generate the requested simuation ----
  # This is called whenever the inputs change. The output functions
  # defined below then use the value computed from this expression
  RESULTS_NEW <- reactive({
    # Get parameters from console
    params <- list()
    # Set the neighbourhood
    params$NGH = input$Neighbourhood
    # Set the IC
    params$IC = input$IC
    # Set main parameters
    params$p_i = input$p_i
    params$R_B = input$R_B
    params$p_r = input$p_r
    params$s_dt = input$s_dt
  
    # # Set the gates (false here)
    # params$PLOT_SIM = FALSE
    # params$PLOT_PROP = FALSE
    return(params)
  })
  
  output$info_ngh <- renderText({
    exp_text = "Two areas can be used."
    exp_text = paste(exp_text, "North River Heights (NRH) is an actual City of Winnipeg neighbourhood.")
    exp_text = paste(exp_text, "Pulberry, Victoria Crescent and Crescent Park (PVC) is the agglomeration of three City of Winnipeg neighbourhoods.")
    exp_text = paste(exp_text, "PVC contains 1513 trees, while NRH contains 2004.")
    exp_text = paste(exp_text, "The dots represent the elm trees present in the area.")
    exp_text = paste(exp_text, "Their size is increasing as a function of their root neighbours, i.e., the trees whose roots are connected with, see the paper for more information.")
    exp_text = paste(exp_text, "")
    exp_text = paste(exp_text, "")
  })
  
  # Display the currently chosen neighbourhood
  output$plot_ngh <- renderImage({
    OUT = RESULTS_NEW()
    if(OUT$NGH == "Mixed Pulberry-Victoria Crescent-Crescent Park (PVC)"){
      filename <- "PCP.png"
    } else if (OUT$NGH == "North River Heights (NRH)"){
      filename <- "NRH.png"
    }

    # Return a list containing the filename and alt text
    list(src = filename,
         width = 500,
         height = 500,
         alt = paste("Neighbourhood : ", input$NGH))
}, deleteFile = FALSE)
  
  
  observeEvent(input$Generate, {
    
  output$plot_maps <- renderPlot({
    OUT = RESULTS_NEW()
    PLOT_SIM  = TRUE
    PLOT_PROP = FALSE
    source("run_sim.R",local = TRUE)
  })
  
  # output$plot_proportions <- renderPlot({
    # OUT = RESULTS_NEW()
    # par(mfrow = c(1,1))
    # prop_Di_Wi = mat.or.vec(nr=1,nc=length(tree_states))
    # for (i in 1:length(prop_Di_Wi)){
    #   prop_Di_Wi[i] = length(which(tree_states[,i]=="I_W" & tree_states[,i]=="I_D"))/length(tree_states[,i])
    # }
    # plot(0:11,prop_Di_Wi) 
    # PLOT_SIM = FALSE
    # PLOT_PROP = TRUE
    # source("C://Users/nicol/Documents/Work/DED-Ecological-Modelling-code/CODE/run_sim.R",local = TRUE)
  # })

  })
}

# Run the application 
shinyApp(ui = ui, server = server)
