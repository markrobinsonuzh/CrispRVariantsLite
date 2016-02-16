################################################################################
# UI
################################################################################


output$plots <- renderUI({
    plotOutput("crispplots")
})


################################################################################
# BEHAVIOUR
################################################################################

# Disable running of the plots until Reference and Bams defined
observe({
   if(is.null(d$ref)){
      updateButton(session, "run_plot", style ="default", icon = icon("ban"), disable = TRUE )
    }else{
      updateButton(session,"run_plot", 'Plot', icon =  icon("area-chart"), style = "success", block = TRUE, disable = FALSE ) 
    }
})


################################################################################
# FUNCTIONS
################################################################################

createCripSet <- reactive({
  if(!is.null(t$DF)){
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  CrispR set ", value = 0)
    
    n <- 20
    
    for (i in 1:10){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = "A")
      Sys.sleep(0.5)
    }
    
    md <- t$DF
    d$ref <- setRef()
    guide <- setGuides()
    
    for (i in 1:10){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = "B")
      Sys.sleep(0.5)
    }

    cset <- readsToTarget(v$bm_fnames, target = guide, reference = d$ref, names = md$label, target.loc = input$target_loc)
    
    return(cset)
  }
})


################################################################################
# PLOT FUNCTION
################################################################################

createCrispPlot <- reactive({
  pcrisp = NULL
  output$crispplots <- renderPlot({
    
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  the plot ", value = 0)
    
    # Number of times we'll go through the loop
    n <- 20
    
    for (i in 1:20){
      #Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = "plotting")
      Sys.sleep(0.5)
    }
     
     try({
       plotVariants(d$cset, txdb = d$txdb, 
      col.ht.ratio = c(1,6),
      left.plot.margin = grid::unit(c(0.1,0,6,0.2), "lines"),
      
      plotFreqHeatmap.args = list(
        top.n = input$top.n,
        min.freq = input$min.freq,
        min.count = input$min.count,
        x.size = input$x.size, 
        plot.text.size = input$plot.text.size, 
        legend.text.size = input$legend.text.size, 
        x.angle = input$x.angle
        ),
      
      plotAlignments.args = list(
        top.n = input$top.n,
        min.freq = input$min.freq,
        min.count = input$min.count,
        axis.text.size = input$axis.text.size, 
        ins.size = input$ins.size, 
        legend.symbol.size = input$legend.symbol.size, 
        legend.text.size = input$legend.text.size
        ), 
        row.ht.ratio = c(1,6)
      )  
     }, silent = TRUE)
     

  })
  return(pcrisp)
})

#create the plot
 observeEvent(input$run_plot,{    
    d$cset <- createCripSet()
    d$txdb <- setTxdb()
    createCrispPlot()
    toggleModal(session, "modal_2", toggle = "close")
  })

