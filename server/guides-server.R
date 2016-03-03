################################################################################
# UI
################################################################################

output$error1 <- renderUI({
  p(paste0("Reference : " , d$ref, " strand : ", input$g.strand))
})

output$guide <- renderUI({
  plotOutput("guide_plot")
})

output$next_step <- renderUI({
    #bsButton("ok ", 'ok', style = "info", block = TRUE)
    actionButton("next_step", "Next step" , width="100%", style="primary")
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
  #f <- paste0("data/txdb/", input$txDb)
  # Note: requires that the file names of the txdb match those of ref genome and
  # this can be achieved by symbolic links
  f <- paste0("./data/txdb/", gsub(".fa",".sqlite",input$select_Refgenome))
  txdb <- AnnotationDbi::loadDb(f)
  return(txdb)
})

setRef <- reactive({
    print("setRef")
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
  switch (input$g.strand,
    "-" =  ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
    "+" =  ref <- Biostrings::DNAString(ref)
  )
    
  print(ref)

  
  updateTextInput(session, "ref_seqs", value = ref)
  
  for (i in 1:7){
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.05)
  }
  return(ref)
})

################################################################################
# PLOTS
################################################################################

creatPlotRef <- reactive({
 output$ref_plot <- renderPlot({
    plot_reference()
 }, height=200)
 
 output$guide_plot <- renderPlot({
   plot_reference()
  }, height = 200)
})


observeEvent(input$run_guide,{
    isolate({
       creatPlotRef()  
    })
 })
 
 observeEvent(input$next_step,{
        toggleModal(session, "modal_ref", toggle = "close")
        toggleModal(session, "modal_2", toggle = "open")
 })
 
plot_reference <- reactive({

        
   ref <- setRef() 
   d$ref <- ref

   box_end <- end(d$guide) - start(d$guide) - d$seq.width + 1

   plotAlignments(
       ref,
       alns = NULL,
       target.loc = d$t.loc,
       guide.loc = IRanges(
         start = d$seq.width + 1,
         end = box_end),
       ins.sites = data.frame(),
       axis.text.size = 14
       )
       
})
