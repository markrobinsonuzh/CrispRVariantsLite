################################################################################
# UI
################################################################################

output$error1 <- renderUI({
  p(paste0("Reference : " , d$ref, " strand : ", input$g.strand))
})

################################################################################
# FUNCTIONS
################################################################################

setGuides <- reactive({
  
  seq.start <- as.numeric(input$g.start)
  seq.end <- input$target_loc+6
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
  
  gd <- setGuides()
  
  genome_index <- paste0(genome,"/",input$select_Refgenome)
  
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

################################################################################
# PLOTS
################################################################################

creatPlotRef <- reactive({
  output$ref_plot <- renderPlot({
    start <- as.numeric(input$g.start)
    plotAlignments(
        setRef(),
        alns = NULL,
        guide.loc = IRanges(start, input$target_loc+6),  
        target.loc = input$target_loc, 
        ins.sites = data.frame())
    })
})

observeEvent(input$run_guide,{
    creatPlotRef()
    toggleModal(session, "modal_ref", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
  })

