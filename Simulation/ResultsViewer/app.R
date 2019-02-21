#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Beta Binomial Coverage Experiment"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
      sidebarPanel(
         sliderInput("Alpha",
                     "Alpha Prior Value",
                     min = 0,
                     max = 5,
                     step=.1, value=0),
         sliderInput("Beta",
                     "Beta Prior Value",
                     min = 0,
                     max = 5,
                     step=.1, value=0),
         sliderInput("Prob", "Probability of success",
                     .01, .99, step=.01, value=.5),
         numericInput("Nmin", "N Minimum", value=10, min=1, max=100),
         numericInput("Nmax", "N Maximum", value=15, min=1, max=100),
         actionButton("Update", "Update")
      ),

      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("Plot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #print(getwd())
  require(ggplot2)
   #load("../ShinyData.RData")

  randomVals <- eventReactive(input$Update, {
    runif(1)
  })

   output$Plot <- renderPlot({
     randomVals()
     isolate({
     range =input$Nmin:input$Nmax
     df=data.frame(N=range, p=input$Prob, Alpha=input$Alpha, Coverage=0)

          for (i in range){
     dataSets=rbinom(10000, i, input$Prob)
     lowers=qbeta(c(.025), input$Alpha+dataSets, input$Beta+i-dataSets)
     uppers=qbeta(c(.975), input$Alpha+dataSets, input$Beta+i-dataSets)
     df$Coverage[df$N==i]<-
       mean(lowers<input$Prob & input$Prob<uppers)
          }
     ggplot(df, aes(x=N, y=Coverage))+geom_point()+geom_line()+ylab("Coverage")+
       geom_hline(yintercept = .95)+geom_hline(yintercept = .945, linetype=2)+
       geom_hline(yintercept = .955, linetype=2)+ggtitle(paste0("Beta(", input$Alpha,", ",input$Beta,
                                                               ") prior with p=", input$Prob))+
       theme(text=element_text(size=15))
     })
   })
}

# data<-meltedCov[meltedCov$Alphas==input$Alpha & meltedCov$ProbSuccess==input$Prob,]
# #print(head(data))
# ggplot(data, aes(x=Ns, y=value))+geom_point()+geom_line()+ylab("Coverage")+
#   geom_hline(yintercept = .95)+geom_hline(yintercept = .94, linetype=2)+
#   geom_hline(yintercept = .96, linetype=2)
# Run the application
shinyApp(ui = ui, server = server)

