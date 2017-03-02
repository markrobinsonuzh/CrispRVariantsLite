#####################################################################################
# Render the final allele summary plot
#####################################################################################


# Disable running of the plots until Reference and Bams defined
observe({
   if(is.null(d$guide)){
      updateButton(session, "run_plot", style ="default", icon = icon("ban"),
                   disable = TRUE )
      updateButton(session, "run_plot_guide", style ="default", icon = icon("ban"),
                   disable = TRUE )
      updateButton(session, "run_guide", 'Create guides', style = "primary",
                   block = TRUE, disable = FALSE)
    }else{
      updateButton(session,"run_plot", 'Plot', icon =  icon("area-chart"),
                   style = "success", block = TRUE, disable = FALSE)
      updateButton(session,"run_plot_guide", 'Plot', icon =  icon("area-chart"),
                  style = "success", block = TRUE, disable = FALSE)
      updateButton(session, "run_guide", "Update guide", style ="info",
                   block = TRUE, disable = FALSE)
    }
})


#________________________________________________________________
# Initialise CrisprSet object

createCrisprSet <- reactive({
  if(!is.null(t$DF)){
    # t$DF is the metadata table
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  CrisprSet object ", value = 0)
    md <- t$DF
    increment_prog(progress, 10, "")
    
    d$cset <- CrispRVariants::readsToTarget(v$bm_fnames, target = d$guide,
                 reference = Biostrings::DNAString(d$ref), names = md$label,
                 target.loc = d$t.loc, verbose = FALSE)
    
    return(d$cset)
  }
})


#________________________________________________________________
# Creation and arrangement of main plot

# Create heatmap
frequency_heatmap <- reactive({
  
  group <- as.factor(t$DF$group)
  
    CrispRVariants::plotFreqHeatmap(d$cset,
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
})


# Create allele plot
allele_plot <- reactive({
    CrispRVariants::plotAlignments(d$cset,
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
    )
})


# Create plots, arrange
createCrispPlot <- reactive({
  
  output$crispplots <- renderPlot({

    ## Warn if the CrisprSet is NULL
    validate(
      need(!is.null(d$cset),
       
           paste(c("CrisprSet could not be created.", 
                   "No on-target reads?"), sep = "\n"))
    )

    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  the plot ", value = 0)
  
    # Number of times we'll go through the loop
    increment_prog(progress, 20, "plotting")
  
    input$replot
  
    isolate({
  
      r_ht <- as.numeric(strsplit(input$row.ht.ratio, ":")[[1]])
      c_wd <- as.numeric(strsplit(input$col.wdth.ratio, ":")[[1]])
  
      plot_margin <- switch(input$plot.margins,
                       "Less" = ggplot2::unit(c(1,0,2,2), "lines"),
                       "Normal" = ggplot2::unit(c(1,0,5,2), "lines"),
                       "More" = ggplot2::unit(c(1,0,10,2), "lines"))
  
      try({
      CrispRVariants:::arrangePlots(
        annotation_plot(), allele_plot(), frequency_heatmap(),
        col.wdth.ratio = c_wd, row.ht.ratio = r_ht,
        left.plot.margin = plot_margin)
      }, silent = TRUE) 
    
    })
  
  })
})

#______________________________________________________________
# Observe button click "run plot", generate and display main plots

# Create the box in which the plot is displayed
output$plots <- renderUI({
    if(is.null(d$cset)) return()
    plotOutput("crispplots", width="auto", height ="600px")
})

# Render plot
 observeEvent(input$run_plot,{
    d$cset <- createCrisprSet()
    print(d$cset)
    d$txdb <- setTxdb()
    createCrispPlot()
    toggleModal(session, "modal_2", toggle = "close")
  })

