

setGuides <- reactive({
  seq.start <- as.numeric(input$g.start)
  seq.end <- seq.start + 22 
  seq.width <- as.numeric(input$g.length)
  
  guide <- GRanges(
    seqnames = input$g.chr, 
    ranges = IRanges(
      start = seq.start, 
      end = seq.end
      ), 
    strand = input$g.strand
    )
  
  guide <- guide + seq.width
  tloc <- 17 + seq.width
  updateTextInput(session, "target_loc", value = paste0(tloc))
  
  return(guide)
})

setTxdb <- reactive({
  f <- paste0("data/txdb/", input$txDb)
  txdb <- AnnotationDbi::loadDb(f)
  return(txdb)
})

setRef <- reactive({
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating  Reference ", value = 0)
  
  n <- 15
  
  for (i in 1:8){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.5)
  }
  
  d$guide <- setGuides()
  gd <- d$guide
  genome_index <- paste0("/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7/", input$select_Refgenome)
  
  #genome_index <- paste0("data/genome/", input$select_Refgenome)
  cmd <- paste0("samtools faidx ", genome_index)
  cmd <- paste0(cmd, " %s:%s-%s")
  ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]
  
  switch (input$g.strand,
    "-" = ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
    "+" = ref <- Biostrings::DNAString(ref)
  )
  
  d$ref <- ref
  
  for (i in 1:7){
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.5)
  }

  return(ref)
})

createCripSet <- reactive({
  if(!is.null(t$DF) && !is.null(d$ref) && !is.null(d$guide)){
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  CrispR set ", value = 0)
    
    n <- 20
    
    for (i in 1:10){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n)
      Sys.sleep(0.1)
    }
    
    md <- t$DF
    
    for (i in 1:10){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n)
      Sys.sleep(0.1)
    }

    cset <- readsToTarget(v$bm_fnames, target = d$guide, reference = d$ref, names = md$label, target.loc = input$target_loc)
    
    return(cset)
  }
})

observe({
  if(is.null(d$ref)){
    updateButton(session, "run_plot", disabled = TRUE, style = "default" )
  }else{
    updateButton(session, "run_plot", disabled = FALSE, style = "success" ) 
  }
})

createCrispPlot <- reactive({
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
        row.ht.ratio = c(1,4)
      )
  })
})

output$plots <- renderUI({
    plotOutput("crispplots")
})


creatPlotRef <- reactive({
  output$ref_plot <- renderPlot({
    plotAlignments(setRef(), alns = NULL, target.loc = input$target_loc, ins.sites = data.frame())
    })
})

output$error1 <- renderUI({
  p(paste0("Reference : " , d$ref, " strand : ", input$g.strand))
})
