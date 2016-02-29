################################################################################
# UI
################################################################################

output$error1 <- renderUI({
  p(paste0("Reference : " , d$ref, " strand : ", input$g.strand))
})

################################################################################
# BEHAVIOUR
################################################################################

# Disable creating the guides until Reference and Bams defined
observe({
   if(is.null(v$bm_fnames)){
      updateButton(session,"create_guides", style ="default", icon = icon("ban"), disable = TRUE )
    }else{
      updateButton(session,"create_guides", style = "success",  icon = icon("area-chart"), block = TRUE, disable = FALSE ) 
    }
})



################################################################################
# FUNCTIONS
################################################################################

setGuides <- reactive({
  print("setting guide")
  
  seq.start <- as.numeric(input$g.start)
  seq.target.loc <- as.numeric(input$target_loc)
  seq.end <-  seq.start + seq.target.loc + 5
  
  guide <- GRanges(
    seqnames = input$g.chr, 
    ranges = IRanges(
      start =  seq.start,
      end = seq.end
      ), 
    strand = input$g.strand
    )
    
  d$seq.width <- as.numeric(input$g.length)
  d$guide <- guide + d$seq.width
  d$t.loc <- seq.target.loc + d$seq.width
  
  return(d$guide)
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
    Sys.sleep(0.05)
  }
  
  gd <- setGuides()
  
  genome_index <- paste0("./data/genome/", input$select_Refgenome)
  cmd <- paste0("samtools faidx ", genome_index)
  cmd <- paste0(cmd, " %s:%s-%s")
 
  
  ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]
  print(sprintf("reference is %s",ref))
  print(sprintf("length of ref is %s", nchar(ref)))
  
  switch (input$g.strand,
    "-" =  d$ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
    "+" =  d$ref <- Biostrings::DNAString(ref)
  )
  

  
  for (i in 1:7){
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.5)
  }
  return( d$ref)
})

################################################################################
# PLOTS
################################################################################

creatPlotRef <- reactive({  
  output$ref_plot <- renderPlot({ 
       
    ref <- setRef() 
    box_end <- end(d$guide) - start(d$guide) - d$seq.width + 1
  
    plotAlignments(
        ref,
        alns = NULL,
        target.loc = d$t.loc,
        guide.loc = IRanges( 
          start = max(d$seq.width + 1), 
          end = box_end),
        ins.sites = data.frame(),
        axis.text.size = 14
        )
    }, height = 200)
    
})

observeEvent(input$run_guide,{
    creatPlotRef()
    toggleModal(session, "modal_ref", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
    d$cset <- createCripSet()
    print(d$cset)
    d$txdb <- setTxdb()
    createHTable()
 })

