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
  seq.target.loc <- as.numeric(input$g.start) + as.numeric(input$target_loc)
  seq.start <-  seq.target.loc - 16
  seq.end <-   seq.target.loc + 6
  seq.width <- as.numeric(input$g.length)
  
  
  guide <- GRanges(
    seqnames = input$g.chr, 
    ranges = IRanges(
      start =  seq.start,
      end = seq.end
      ), 
    strand = input$g.strand
    )
    
  d$guide <- guide + seq.width
 # tloc <- as.numeric(input$target_loc) + seq.width
  
 # updateTextInput(session, "target_loc", value = paste0(tloc))
  
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
    Sys.sleep(0.5)
  }
  
  gd <- setGuides()
  
  genome_index <- paste0("./data/genome/", input$select_Refgenome)
  cmd <- paste0("samtools faidx ", genome_index)
  cmd <- paste0(cmd, " %s:%s-%s")
 
  
  ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]
  
  
  switch (input$g.strand,
    "-" =  d$ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
    "+" =  d$ref <- Biostrings::DNAString(ref)
  )
  
  updateTextInput(session, "ref_seqs", value = d$ref)
  
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
      t.loc <- as.numeric(input$target_loc)
    plotAlignments(
        setRef(),
        alns = NULL,
        target.loc = input$target_loc,
        #guide.loc = IRanges(
        #    start = 
        #    end = t.loc +6,
        #    width = t.loc + 6,
         #   ),  
        ins.sites = data.frame()
        )
    })
    
})

observeEvent(input$run_guide,{
    creatPlotRef()
    toggleModal(session, "modal_ref", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
    createHTable()
 })

