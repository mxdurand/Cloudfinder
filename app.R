library(shiny)
library(bslib)

# Define UI ----
ui <- fluidPage(
  titlePanel("Cloudfinder"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "inputFile",
        label = "File upload:",
        multiple = FALSE,
        accept = ".csv",
        buttonLabel = "Browse...",
        placeholder = "No file selected",
      ),
      sliderInput(inputId = "CSthres",
                  label = "Clear Sky Threshold (%):",
                  min = 0,
                  max = 100,
                  value = 70
      ),
      sliderInput(inputId = "CStimelimit",
                  label = "Clear Sky Time limit (min):",
                  min = 0,
                  max = 100,
                  value = 45
      ),
      sliderInput(inputId = "OVthres",
                  label = "Overcast Threshold (%):",
                  min = 0,
                  max = 100,
                  value = 60
      ),
      sliderInput(inputId = "OVtimelimit",
                  label = "Overcast Time limit (min):",
                  min = 0,
                  max = 100,
                  value = 30
      ),
      selectInput(inputId = "leftovers",
                  label = "Second pass for remaining data?",
                  choices = c("TRUE", "FALSE"),
                  selected = "TRUE",
                  multiple = FALSE
      ),
      actionButton("go", "Plot")
    ),
    mainPanel(
      plotOutput(outputId = "plot")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  v <- reactiveValues(doPlot = FALSE)
  
  observeEvent(input$go, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$go
  })
  
  observeEvent(input$tabset, {
    v$doPlot <- FALSE
  })  
  
  output$plot <- renderPlot({
    if (v$doPlot == FALSE) return()
    
    isolate({
      source("Cloudfinder.R")
      l <- procClassif(file = input$inputFile$datapath, input$CSthres/100, input$OVthres/100, input$CStimelimit, input$OVtimelimit, leftovers = ifelse(input$leftovers == "TRUE", TRUE, FALSE))
      
      dfi <- l[[1]]
      dfj <- l[[2]]
      
      par(mar = c(4,5,2,0.5), bty = "L")
      plot(-500, ylim = c(-50,2500), xlim = c(0, 86400), xaxt = "n", yaxt = "n", xlab  = "", ylab = "", main = paste(dfi$date[1]), font.main = 4)
      
      for(j in 1:nrow(dfj))
      {
        if(is.na(dfj[j,"sky"]) | dfj[j,"sky"] == ""){
          icol = "white"
        } else if(dfj[j,"sky"] == "AMB"){
          icol = "lightgreen"
        } else if(dfj[j,"sky"] == "PC"){
          icol = "#BFEFFF"
        } else if(dfj[j,"sky"] == "CS"){
          icol = "#FFEC8B"
        } else if(dfj[j,"sky"] == "OV"){
          icol = "gray90"
        }
        
        itsec <- sapply(strsplit(dfj[j,"hour"], split = ":"), function(x) {x <- as.numeric(x) ; x[1]*60*60 + x[2]*60 + x[3]})
        itsec2 <- itsec + 3600
        polygon(x = c(itsec, itsec, itsec2, itsec2), y = c(-100, 3000, 3000, -100), border = NA, density = NULL, col = icol)
      }
      abline(v = c(7*60*60, 19*60*60), col = "gray30")
      abline(v = c(13*60*60), col = "gray30", lty = 5)
      # PARlimit <- dfi$PARscl * growinglimit(dfi$PARscl)
      # points(PARlimit ~ tsec, dfi, type = "l", col = "green4", lwd = 2, lty = 3)
      points(PARscl ~ tsec, dfi, type = "l", col = "red", lwd = 2, lty = 3)
      points(PARmes ~ tsec, dfi, type = "l", col = "black")
      axis(side = 1, font = 2, at = seq(3600,82800, 2*3600), labels = seq(1,23,2))
      axis(side = 2, font = 2, las = 2)
      mtext (side = 1, line = 2.5, text = "Time (h)", font = 4, cex = 1.2)
      mtext (side = 2, line = 3.5, text = "PPFD", font = 4, cex = 1.2)
      
    })
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)