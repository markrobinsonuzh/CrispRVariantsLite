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
   if(is.null(d$guide)){
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
      Sys.sleep(0.005)
    }
    
    md <- t$DF
    
    for (i in 1:10){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = "B")
      Sys.sleep(0.005)
    }
    print(d$guide)
    
    d$cset <- readsToTarget(v$bm_fnames, target = d$guide,
                 reference = setRef(), names = md$label, 
                 target.loc = input$target_loc, verbose = FALSE)
    
    return(d$cset)
  }
})


################################################################################
# PLOT FUNCTION
################################################################################

createCrispPlot <- reactive({
  pcrisp = NULL
  
  output$crispplots <- renderPlot({
    
    
    ## Warn if the CrisprSet is NULL
    validate(
      need(!is.null(d$cset), 
           paste(c("CrisprSet could not be created.", 
                   "No on-target reads?"), sep = "\n"))
    )
    
    input$replot
    
    isolate({
    
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Creating  the plot ", value = 0)
    
      # Number of times we'll go through the loop
      n <- 20
    
      for (i in 1:20){
        #Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = "plotting")
        Sys.sleep(0.005)
      }
    
      r_ht <- as.numeric(strsplit(input$row.ht.ratio, ":")[[1]])
      c_wd <- as.numeric(strsplit(input$col.wdth.ratio, ":")[[1]])
     
      try({
        group <- as.factor(t$DF$group)

        plotVariants(d$cset, txdb = d$txdb, 
          col.wdth.ratio = c_wd, row.ht.ratio = r_ht,
          gene.text.size = 8,
          left.plot.margin = ggplot2::unit(c(1,0,5,2), "lines"),
      
          plotAlignments.args = list(
            top.n = input$top.n,
            min.freq = input$min.freq,
            min.count = input$min.count,
            target.loc = d$t.loc,
            guide.loc = IRanges(
            start = d$seq.width + 1,
            end = end(d$guide) - start(d$guide) - d$seq.width + 1),
            axis.text.size = input$axis.text.size, 
            ins.size = input$ins.size,
            plot.text.size = input$plot.text.size, 
            legend.symbol.size = input$legend.symbol.size,
            legend.text.size = input$legend.text.size
          ), 
        
          plotFreqHeatmap.args = list(
            top.n = input$top.n,
            min.freq = input$min.freq,
            min.count = input$min.count,
            type =  input$plot.type,
            x.size = input$x.size, 
            plot.text.size = input$plot.text.size, 
            legend.text.size = input$legend.text.size,
            x.angle = input$x.angle,
            group = group
          )
        ) + theme(legend.position="none")  
      }, silent = TRUE) 
    })
  }, height = 600)
  return(pcrisp)
})


# create the plot
 observeEvent(input$run_plot,{  
    d$cset <- createCripSet()
    print(d$cset)
    d$txdb <- setTxdb()
    createCrispPlot()
    toggleModal(session, "modal_2", toggle = "close")  
  })
