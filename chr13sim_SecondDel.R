#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Define UI for application that draws a MAD-like plot of the UPD13
ui <- shinyUI(ui = fluidPage(
   
   # Application title
   titlePanel("UPD13 Cellularity Simulator"),
   
   # Sidebar with a slider input for the different cellularities and parameters
   sidebarLayout(
      sidebarPanel(
         textInput(inputId = "sample",
                   label = "Sample ID:",
                   value = ""),
         selectInput(inputId = "type",
                     label = "Type of plot:", 
                     choices=c("LRR from -1 to 1", "LRR from -2 to 2")),
         numericInput(inputId = "nprobes",
                     label = "number of probes in the array:",
                     value = 36590),
         numericInput(inputId = "LRRdev",
                      label = "Log R Ratio deviation:",
                      value = 0.2, 
                      step = 0.01),
         numericInput(inputId = "BAFdev",
                      label = "B Allele Frequency deviation:",
                      value = 0.03, 
                      step = 0.01),
         numericInput(inputId = "CelDelHet",
                      label = "Cellularity of the heterozygous deletion:",
                      min = 0,
                      max = 100,
                      value = 0,
                      step=1),
         numericInput(inputId = "DelLBP",
                      label = "Deletion left Breakpoint:",
                      value = 49500000),
         numericInput(inputId = "DelRBP",
                      label = "Deletion rigth Breakpoint:",
                      value = 50500000),
         checkboxInput(inputId = "DelHet2",
                       label = strong("Second heterozygous deletion"),
                       value = FALSE),
         conditionalPanel(condition = "input.DelHet2 == true",
                          numericInput(inputId = "CelDelHet2",
                                       label = "Cellularity of the second heterozygous deletion:",
                                       min = 0, max = 100, value = 0, step = 1),
                          numericInput(inputId = "DelLBP2",
                                       label = "2nd hit deletion left Breakpoint:",
                                       value = 49500000),
                          numericInput(inputId = "DelRBP2",
                                       label = "2nd hit deletion rigth Breakpoint:",
                                       value = 50500000)),
         checkboxInput(inputId = "UPDdel",
                       label = strong("Homozygous deletion by UPD"),
                       value = FALSE),
         conditionalPanel(condition = "input.UPDdel == true",
                          numericInput(inputId = "CelUPDDel",
                                      label = "Cellularity of the Homozygous deletion by UPD:",
                                      min = 0, max = 100, value = 0, step = 1),
                          numericInput(inputId = "UPDHomoDelbreak",
                                       label = "Breakpoint of the Homozygous deletion by UPD:",
                                       min = 17920000, max = 114100000, value = 22000000, step = 1000000),
                          numericInput(inputId = "CelUPDDelRec",
                                       label = "Cellularity of the reciprocal homozygous deletion by UPD:",
                                       min = 0, max = 100, value = 0, step = 1)),
         checkboxInput(inputId = "UPDres",
                       label = strong("Rescue UPD"),
                       value = FALSE),
         conditionalPanel(condition = "input.UPDres == true",
                          numericInput(inputId = "CelUPDRes",
                                       label = "Cellularity of the rescue UPD:",
                                       min = 0, max = 100, value = 0, step = 1),
                          numericInput(inputId = "UPDRescuebreak",
                                       label = "Breakpoint of the rescue UPD:",
                                       min = 17920000, max = 114100000, value = 17920000, step = 1000000),
                          numericInput(inputId = "CelUPDResRec",
                                       label = "Cellularity of the reciprocal rescue UPD:",
                                       min = 0, max = 100, value = 0, step = 1),
                          checkboxInput(inputId = "UPDres2",
                                        label = strong("Second Rescue UPD"),
                                        value = FALSE),
                          conditionalPanel(condition = "input.UPDres2 == true",
                                           numericInput(inputId = "CelUPDRes2",
                                                        label = "Cellularity of the second rescue UPD:",
                                                        min = 0, max = 100, value = 0, step = 1),
                                           numericInput(inputId = "UPDRescuebreak2",
                                                        label = "Breakpoint of the second rescue UPD:",
                                                        min = 17920000, max = 114100000, value = 19920000, step = 1000000),
                                           numericInput(inputId = "CelUPDRes2Rec",
                                                        label = "Cellularity of the reciprocal second rescue UPD:",
                                                        min = 0, max = 100, value = 0, step = 1)))
      ),
      
      # Show a plot of the simulated mosaic, a zoom on the deletion region and a table with the cellularity values
      mainPanel(
         plotOutput("MosaicPlot"),
         plotOutput("MosaicDeletionPlot"),
         dataTableOutput("table")
      )
   )
))

# Define server logic required to draw the plots
server <- shinyServer(function(input, output, session) {
  BAFcalc <- function(n, Bdev, BAFdev){
    BAFAB <- c(rnorm(100000, 0.5 - Bdev, BAFdev), rnorm(100000, 0.5+Bdev, BAFdev))
    res <- sample(BAFAB, n, replace = TRUE)
    return(res)
  }
  UPDCelCalc <- function(x) 100*2*x
  UPDBdevCalc <- function(x) x/(2*100)
  LRRCalc <- function(x) log(x/2)/1.5
 
   observe({
    if(!input$UPDdel){
      updateNumericInput(session,
                         inputId = "CelUPDDel",
                         value = 0)
      updateNumericInput(session,
                         inputId = "CelUPDDelRec",
                         value = 0)
    }
    if(!input$DelHet2){
      updateNumericInput(session,
                         inputId = "CelDelHet2",
                         value = 0)
      # updateNumericInput(session,
      #                    inputID = "DelLBP2",
      #                    value = 49500000)
      # updateNumericInput(session,
      #                    inputID = "DelRBP2",
      #                    value = 50500000)
    }
    if(!input$UPDres) {
      updateNumericInput(session,
                         inputId = "CelUPDRes",
                         value = 0)
      updateNumericInput(session,
                         inputId = "CelUPDResRec",
                         value = 0)
    } 
    if(!input$UPDres2) {
      updateNumericInput(session,
                         inputId = "CelUPDRes2",
                         value = 0)
      updateNumericInput(session,
                         inputId = "CelUPDRes2Rec",
                         value = 0)
    }
  })
  
   output$MosaicPlot <- renderPlot({
      # generate interactive plot based on inputs from ui.R
      set.seed(1113)
      pos <- sort(runif(input$nprobes, min=17920000, max=114100000))
      ChrUPDDel <- input$CelUPDDel*2
      ChrDelHet2 <- input$CelDelHet2
      ChrDelHet <- input$CelDelHet
      ChrDel <- (ChrDelHet + ChrUPDDel + ChrDelHet2)
      LRRDel <- LRRCalc(2 - ChrDel/100)
      LRRDel1 <- LRRCalc(2 - (ChrDelHet + ChrUPDDel)/100)
      LRRDel2 <- LRRCalc(2 - (ChrDelHet2)/100)
      CelWT <- 100 - input$CelDelHet - input$CelUPDDel - input$CelUPDDelRec - input$CelUPDRes - input$CelUPDResRec - input$CelDelHet2- input$CelUPDRes2 - input$CelUPDRes2Rec

      BAFAA <- rnorm(100000, 0, input$BAFdev)
      BAFAA[BAFAA < 0] <- abs(BAFAA[BAFAA < 0])
      BAFAB <- rnorm(100000, 0.5, input$BAFdev)
      BAFBB <- rnorm(100000, 1, input$BAFdev)
      BAFBB[BAFBB > 1] <- 2 - BAFBB[BAFBB > 1]
      
      BdevDel <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet + 2*input$CelUPDDelRec + 2*input$CelUPDResRec + 2*input$CelUPDRes2Rec)) # Depèn de molts factors el BAF de la delecció
      BdevDel1 <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet2) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDResRec + 2*input$CelUPDRes2 + 2*input$CelUPDRes2Rec + input$CelDelHet + 2*input$CelUPDDelRec + 2*input$CelDelHet2))
      BdevDel2 <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet + 2*input$CelUPDDel) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDResRec + 2*input$CelUPDRes2 + 2*input$CelUPDRes2Rec + 2*input$CelDelHet + 2*input$CelUPDDel + 2*input$CelUPDDelRec + input$CelDelHet2))
      BdevUPDDel <- UPDBdevCalc(input$CelUPDDel - input$CelUPDDelRec)
      BdevUPDRes <- UPDBdevCalc(input$CelUPDRes - input$CelUPDResRec)
      BdevUPDRes2 <- UPDBdevCalc(input$CelUPDRes2 - input$CelUPDRes2Rec)
      
      BAF <- sample(c(BAFAA, BAFAB, BAFBB), size = length(pos), replace = T)
      if (input$UPDHomoDelbreak >= input$UPDRescuebreak){
        if (input$UPDHomoDelbreak >= input$UPDRescuebreak){
          BAF[pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDRes, input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2), replace = T)
          BAF[pos >= input$UPDRescuebreak2 & pos < input$UPDHomoDelbreak] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDRes2), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak2 & pos < input$UPDHomoDelbreak), replace = T)
        } else {
          BAF[pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDRes, input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak), replace = T)
          BAF[pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDDel), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak), replace = T)
        }
      } else {
        BAF[pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDDel, input$BAFdev), BAFBB), size = sum(pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak), replace = T)
        BAF[pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDDel+BdevUPDRes), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2), replace = T)
      }
      BAF[pos >= input$UPDRescuebreak & pos >= input$UPDHomoDelbreak & pos >= input$UPDRescuebreak2 ] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDDel+BdevUPDRes2), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos >= input$UPDHomoDelbreak & pos >= input$UPDRescuebreak2), replace = T)
      
      del1 <- pos >= input$DelLBP & pos <= input$DelRBP
      del2 <- pos >= input$DelLBP2 & pos <= input$DelRBP2
      BAF[del1 & del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel, input$BAFdev), BAFBB), size = sum(del1 & del2), replace = T)
      BAF[del1 & !del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel1, input$BAFdev), BAFBB), size = sum(del1 & !del2), replace = T)
      BAF[!del1 & del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel2, input$BAFdev), BAFBB), size = sum(!del1 & del2), replace = T)
      #if (input$Homdel) BAF[pos > input$DelLBP2 & pos <= input$DelRBP2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel2, input$BAFdev), BAFBB), size = sum(pos > input$DelLBP2 & pos <= input$DelRBP2), replace = T)
      
      BAF[BAF < 0] <- abs(BAF[BAF < 0])
      BAF[BAF > 1] <- (2 - BAF[BAF > 1])
      
      LRR <- rnorm(length(pos), 0, input$LRRdev)

      LRR[del1 & del2] <- rnorm(sum(del1 & del2), LRRDel, input$LRRdev)
      LRR[del1 & !del2] <- rnorm(sum(del1 & !del2), LRRDel1, input$LRRdev)
      LRR[!del1 & del2] <- rnorm(sum(!del1 & del2), LRRDel2, input$LRRdev)
      
      par(mar = c(5, 4, 4, 4) + 0.1)
      if (input$type == "LRR from -1 to 1") {
        plot(pos, LRR, ylim = c(-1, 1), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = input$sample, xaxt="n", yaxt="n")
        axis( 2, at = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", -0.5, "0.0", 0.5, "1.0"), las = 1, col = "black", col.axis = "black")
      } else {
        plot(pos, LRR, ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = input$sample, xaxt="n", yaxt="n")
        axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c("-2.0", "-1.5", "-1.0", -0.5, "0.0", 0.5, "1.0", "1.5", "2.0"), las = 1, col = "black", col.axis = "black")
      }
      par(new = TRUE)
      plot(pos, BAF, col = 2, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F)
      abline(h=0.5, col=8)
      abline(h=c(0.33, 0.66), col=8)
      
      mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
      mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
      axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), las = 1, col = "black", col.axis = "black")
      xaxis <- seq(min(pos), max(pos), length.out=10)
      axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
      mtext("position (Mb)", side = 1, col = "black", line = 2.5)
   })
   output$MosaicDeletionPlot <- renderPlot({
     # generate interactive plot based on inputs from ui.R
     set.seed(1113)
     pos <- sort(runif(input$nprobes, min=17920000, max=114100000))
     if(input$DelHet2){
       xlimits <- c(min(c(input$DelLBP, input$ DelLBP2))-500000, max(c(input$DelRBP, input$DelRBP2))+500000)
     } else {
       xlimits <- c(min(input$DelLBP)-500000, max(input$DelRBP)+500000)
     }
     
     pos <- sort(runif(input$nprobes, min=17920000, max=114100000))
     ChrUPDDel <- input$CelUPDDel*2
     ChrDelHet2 <- input$CelDelHet2
     ChrDelHet <- input$CelDelHet
     ChrDel <- (ChrDelHet + ChrUPDDel + ChrDelHet2)
     LRRDel <- LRRCalc(2 - ChrDel/100)
     LRRDel1 <- LRRCalc(2 - (ChrDelHet + ChrUPDDel)/100)
     LRRDel2 <- LRRCalc(2 - (ChrDelHet2)/100)
     CelWT <- 100 - input$CelDelHet - input$CelUPDDel - input$CelUPDDelRec - input$CelUPDRes - input$CelUPDResRec - input$CelDelHet2- input$CelUPDRes2 - input$CelUPDRes2Rec
     
     BAFAA <- rnorm(100000, 0, input$BAFdev)
     BAFAA[BAFAA < 0] <- abs(BAFAA[BAFAA < 0])
     BAFAB <- rnorm(100000, 0.5, input$BAFdev)
     BAFBB <- rnorm(100000, 1, input$BAFdev)
     BAFBB[BAFBB > 1] <- 2 - BAFBB[BAFBB > 1]
     
     BdevDel <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet + 2*input$CelUPDDelRec + 2*input$CelUPDResRec + 2*input$CelUPDRes2Rec)) # Depèn de molts factors el BAF de la delecció
     BdevDel1 <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet2) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDResRec + 2*input$CelUPDRes2 + 2*input$CelUPDRes2Rec + input$CelDelHet + 2*input$CelUPDDelRec + 2*input$CelDelHet2))
     BdevDel2 <- abs(0.5 - (CelWT + 2*input$CelUPDRes + 2*input$CelUPDRes2 + input$CelDelHet + 2*input$CelUPDDel) / (2*CelWT + 2*input$CelUPDRes + 2*input$CelUPDResRec + 2*input$CelUPDRes2 + 2*input$CelUPDRes2Rec + 2*input$CelDelHet + 2*input$CelUPDDel + 2*input$CelUPDDelRec + input$CelDelHet2))
     BdevUPDDel <- UPDBdevCalc(input$CelUPDDel - input$CelUPDDelRec)
     BdevUPDRes <- UPDBdevCalc(input$CelUPDRes - input$CelUPDResRec)
     BdevUPDRes2 <- UPDBdevCalc(input$CelUPDRes2 - input$CelUPDRes2Rec)
     
     BAF <- sample(c(BAFAA, BAFAB, BAFBB), size = length(pos), replace = T)
     if (input$UPDHomoDelbreak >= input$UPDRescuebreak){
       if (input$UPDHomoDelbreak >= input$UPDRescuebreak){
         BAF[pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDRes, input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2), replace = T)
         BAF[pos >= input$UPDRescuebreak2 & pos < input$UPDHomoDelbreak] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDRes2), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak2 & pos < input$UPDHomoDelbreak), replace = T)
       } else {
         BAF[pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDRes, input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak), replace = T)
         BAF[pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDDel), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDHomoDelbreak), replace = T)
       }
     } else {
       BAF[pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak] <- sample(c(BAFAA, BAFcalc(100000, BdevUPDDel, input$BAFdev), BAFBB), size = sum(pos >= input$UPDHomoDelbreak & pos < input$UPDRescuebreak), replace = T)
       BAF[pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDDel+BdevUPDRes), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos < input$UPDRescuebreak2), replace = T)
     }
     BAF[pos >= input$UPDRescuebreak & pos >= input$UPDHomoDelbreak & pos >= input$UPDRescuebreak2 ] <- sample(c(BAFAA, BAFcalc(100000, abs(BdevUPDRes+BdevUPDDel+BdevUPDRes2), input$BAFdev), BAFBB), size = sum(pos >= input$UPDRescuebreak & pos >= input$UPDHomoDelbreak & pos >= input$UPDRescuebreak2), replace = T)
     
     del1 <- pos >= input$DelLBP & pos <= input$DelRBP
     del2 <- pos >= input$DelLBP2 & pos <= input$DelRBP2
     BAF[del1 & del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel, input$BAFdev), BAFBB), size = sum(del1 & del2), replace = T)
     BAF[del1 & !del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel1, input$BAFdev), BAFBB), size = sum(del1 & !del2), replace = T)
     BAF[!del1 & del2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel2, input$BAFdev), BAFBB), size = sum(!del1 & del2), replace = T)
     #if (input$Homdel) BAF[pos > input$DelLBP2 & pos <= input$DelRBP2] <- sample(c(BAFAA, BAFcalc(100000, BdevDel2, input$BAFdev), BAFBB), size = sum(pos > input$DelLBP2 & pos <= input$DelRBP2), replace = T)
     
     BAF[BAF < 0] <- abs(BAF[BAF < 0])
     BAF[BAF > 1] <- (2 - BAF[BAF > 1])
     
     LRR <- rnorm(length(pos), 0, input$LRRdev)
     
     LRR[del1 & del2] <- rnorm(sum(del1 & del2), LRRDel, input$LRRdev)
     LRR[del1 & !del2] <- rnorm(sum(del1 & !del2), LRRDel1, input$LRRdev)
     LRR[!del1 & del2] <- rnorm(sum(!del1 & del2), LRRDel2, input$LRRdev)
     par(mar = c(5, 4, 4, 4) + 0.1)
     if (input$type == "LRR from -1 to 1") {
       plot(pos, LRR, ylim = c(-1, 1), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = input$sample, xaxt="n", yaxt="n", xlim=xlimits)
       axis( 2, at = c(-1, -0.5, 0, 0.5, 1), labels = c("-1.0", -0.5, "0.0", 0.5, "1.0"), las = 1, col = "black", col.axis = "black")
     } else {
       plot(pos, LRR, ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = input$sample, xaxt="n", yaxt="n", xlim=xlimits)
       axis( 2, at = c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c("-2.0", "-1.5", "-1.0", -0.5, "0.0", 0.5, "1.0", "1.5", "2.0"), las = 1, col = "black", col.axis = "black")
     }
     par(new = TRUE)
     plot(pos, BAF, col = 2, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, xlim=xlimits)
     abline(h=0.5, col=8)
     abline(h=c(0.33, 0.66), col=8)
     
     mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
     mtext("BAF", side = 4, col = "black", line = 2.5, adj = 0.5)
     axis( 4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", 0.2, 0.4, 0.6, 0.8, "1.0"), las = 1, col = "black", col.axis = "black")
     xaxis <- seq(xlimits[1], xlimits[2], length.out=10)
     axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
     mtext("position (Mb)", side = 1, col = "black", line = 2.5)
   })
   
   output$table <- renderDataTable({
     pos <- sort(runif(input$nprobes, min=17920000, max=114100000))
     ChrUPDDel <- input$CelUPDDel*2
     ChrDelHet2 <- input$CelDelHet2
     ChrDelHet <- input$CelDelHet
     ChrDel <- (ChrDelHet + ChrUPDDel + ChrDelHet2)
     LRRDel <- LRRCalc(2 -ChrDel/100)
     LRRDel1 <- LRRCalc(2 - (ChrDelHet + ChrUPDDel)/100)
     LRRDel2 <- LRRCalc(2 - (ChrDelHet2)/100)
     CelWT <- 100 - input$CelDelHet - input$CelUPDDel - input$CelUPDDelRec - input$CelUPDRes - input$CelUPDResRec - input$CelDelHet2 - input$CelUPDRes2 - input$CelUPDRes2Rec
     #BdevDel <-  abs(0.5 - (CelWT + 2*input$CelUPDRes) / (2*CelWT + 2*input$CelUPDRes + input$CelDelHet)) # Depends on a lot of elements
     #BdevUPDDel <- UPDBdevCalc(abs(input$CelUPDDel - input$CelUPDDelRec))
     #BdevUPDRes <- UPDBdevCalc(abs(input$CelUPDRes - input$CelUPDResRec))
     data <- data.frame(WT=round(CelWT, 2), 
                        HetDel=round(input$CelDelHet, 2), 
                        HetDel2=round(input$CelDelHet2, 2), 
                        UPDDel=round(input$CelUPDDel, 2),
                        UPDDelRec=round(input$CelUPDDelRec, 2), 
                        UPDRes=round(input$CelUPDRes, 2),
                        UPDResRec=round(input$CelUPDResRec, 2), 
                        UPDRes2=round(input$CelUPDRes2, 2),
                        UPDRes2Rec=round(input$CelUPDRes2Rec, 2),
                        LRRDel1=round(LRRDel1, 2),
                        LRRDel2=round(LRRDel2, 2),
                        LRRDel=round(LRRDel, 2))
                        #UPDDelRec=round(input$CelUPDDelRec, 2),
                        #UPDResRec=round(input$CelUPDResRec, 2),
                        #UPDRes2Rec=round(input$CelUPDRes2Rec, 2))
     data
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

