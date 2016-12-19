################################################################################
# UI
################################################################################

output$error1 <- renderUI({
    tags$div(
        p(paste0("coordinates : " , input$g.start, " strand : ", input$g.strand, " chr : ", input$g.chr))
    )
})

output$guide <- renderUI({
    tags$div(
        plotOutput("guide_plot", width="100%", height="200px"),
        uiOutput("error1")
    )
    
})

output$next_step <- renderUI({
    bsButton("next_step", "Back" , type="action", style = "default", block = TRUE)
})


################################################################################
# BEHAVIOUR
################################################################################

# Disable creating the guides until Reference and Bams defined
observe({
   if(is.null(v$bm_fnames)){
      updateButton(session,"create_guides", style ="default", icon = icon("ban"), disable = TRUE )
    }else{
      updateButton(session,"create_guides", style = "success",  icon = icon("area-chart"),
                   block = TRUE, disable = FALSE ) 
    }
})

# open reference modal
  observeEvent(input$create_guides, {
   #toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_ref", toggle = "open")
  })

 observeEvent(input$next_step,{
      toggleModal(session, "modal_ref", toggle = "close")
     #toggleModal(session, "modal_2", toggle = "open")
 })
 
################################################################################
# FUNCTIONS
################################################################################

setGuideCoords <- function(guide = NULL, window = NULL){
  # Set guide coordinates including window around guide
  
  if(input$run_guide == 0) return()
 
  seq.start <- as.numeric(input$g.start)
  target.loc <- as.numeric(input$target_loc)
  seq.end <-  seq.start + target.loc + 5
  
  if (is.null(guide)){
  guide <- GRanges(
      seqnames = input$g.chr, 
      ranges = IRanges(start =  seq.start, end = seq.end),
      strand = input$g.strand
    )
  }   

  # Add window around guide, adjust target location 
  d$seq.width <- ifelse(is.null(window), as.numeric(input$g.length), window)
  d$guide <- guide + d$seq.width
  d$t.loc <- target.loc + d$seq.width
      
  return(d$guide)
}


setTxdb <- reactive({
  # Note: requires that the file names of the txdb match those of ref genome and
  # this can be achieved by symbolic links
  f <- paste0("./data/txdb/", gsub(".fa",".sqlite",input$select_Refgenome))
  d$txdb <- AnnotationDbi::loadDb(f)
  return(d$txdb)
})


setGuidesFromCoords <- function(genome_index, progress){
    # Set the guide and fetch the guide sequence directly from the reference    
 
    gd <- setGuideCoords()
    cmd <- paste0("samtools faidx ", genome_index)
    cmd <- paste0(cmd, " %s:%s-%s")
    
    # Try to extract the sequence matching this location from the reference    
    err <- try({
        ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), 
                      intern = T)
        # ref is returned in fasta format, may be multi-line
        ref <- paste0(ref[2:length(ref)], collapse = "")
    }, silent = TRUE)
 
     # Error if sequence not found
    isolate({
        if(is.error(err)){
          createAlert(session, "alertRef", alertId = "alertRef3", title = "INFO",
            content = paste0("the sequence starting at the coordinate : ", input$g.start,
                             " in ", input$select_Refgenome," at ", input$g.chr, " not found " ),
          style = "info", dismiss = TRUE, append = FALSE)
          return(NULL)
        }
    })
            
    increment_prog(progress, 5, "finding strand")                  
                
    # Reverse complement the sequence if a negative strand was specified
    switch (input$g.strand,
         "-" =  ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
         "+" =  ref <- Biostrings::DNAString(ref))
  
   #updateTextInput(session, "ref_seqs", value = toString(ref)) 
   return(ref)
}


writeGuideFasta <- function(ref){
  fa <- paste0(MHmakeRandomString(),gsub("[- :]", "", Sys.time()),"_reference.fasta")
  fa <- file.path(v$fasta_temp, fa)
  write.fasta(sequences = ref, names = "reference", file.out = fa )
  return(fa)
}


mapGuide <- function(ref, idx){
   fa <- writeGuideFasta(ref) 

   cmd <- paste("bwa aln %s %s | bwa samse %s - %s | grep -v '^@' | awk -F \"\t\"",
              "'{if ($2 == 0)print $3, $4, length($10), \"+\";",
              "else if ($2 == 16) print $3, $4, length($10), \"-\"}' && rm %s")

   result <- NULL
        
   err <- try({
     result <- strsplit(system(sprintf(cmd, idx, fa, idx, fa, fa), intern = TRUE), " ")[[1]]
     print(result)
   })

   isolate({
     if(is.error(err)){
       createAlert(session, "alertRef", alertId = "alertRef2", title = "INFO",
         content = paste0("specified sequence : ", ref , " does not exist in ",
                        input$select_Refgenome),
         style = "warning",  dismiss = TRUE, append = FALSE)
       return(NULL)
     }  
   })
  return(result)
}




observeEvent(input$run_guide,{
  
    # If bwa cannot be run, stop the app
    if (system("bwa", ignore.stderr = TRUE) == 127){
      createAlert(session, "alertRef","alertRef5", title = "WARNING",
         content = "Please make sure BWA is available", 
         style = "warning",
         append = FALSE, dismiss = FALSE)
      stopApp()
    }
  
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  Reference ", value = 0)
  
    n <- 20
 
    increment_prog(progress, 3, "Prepare for mapping sequence")    
 
    ref <- NULL
    genome_index <- paste0("./data/genome/", input$select_Refgenome)
    
    # If text was entered into the coordinates list, prioritise over
    # a sequence entered   
    if(nchar(input$g.start) > 0 ){
    
      ref <- setGuidesFromCoords(genome_index, progress)
      if (is.null(ref)) return()

    } else if (nchar(input$ref_seqs) >=  10 ) {
      
      # This update to slider doesn't seem to be propagated until
      # return from function, possibly because of tabs
      updateSliderInput(session, "g.length", value = 0)
      increment_prog(progress, 2, "getting coordinates")

       # Map input sequence and parse the results
      idx <- genome_index
      ref <- input$ref_seqs
      
      result <- mapGuide(ref, idx)
      if (is.null(result)) return()

      chr <- result[1]
      sq.start <- as.numeric(result[2])
      sq.length <- as.numeric(result[3])
      strd <- result[4]
      
      # Update the text input
      increment_prog(progress, 3, "updating text input")

      updateTextInput(session, "g.chr", value = paste(chr))
      updateTextInput(session, "g.start", value = paste(sq.start))
      updateTextInput(session, "g.strand", value = paste(strd))
      
      increment_prog(progress, 2, "creating the guide") 
      
      # Store the guide coordinates
      guide <- GenomicRanges::GRanges(chr,
             IRanges(sq.start , width = sq.length), strand = strd)

      setGuideCoords(guide, 0)
      ref <- Biostrings::DNAString(ref)
    
    } else{
        createAlert(session, "alertRef", alertId = "alertRef1", title = "WARNING",
        content = "Please enter the sequence or set the coordinates", style = "info",
        dismiss = TRUE, append = FALSE)
        
        return()
    }
        
  increment_prog(progress, 7, "")  
  
  d$ref <- Biostrings::DNAString(ref)
  
})


################################################################################
# PLOTS
################################################################################

observeEvent(input$run_plot_guide,{    

    d$cset <- createCrisprSet()
    d$txdb <- setTxdb()
    createCrispPlot()
    toggleModal(session, "modal_ref", toggle = "close")
    Sys.sleep(0.05)
    toggleModal(session, "modal_2", toggle = "close")
  })


  reference_plot <- reactive({
   
    if(is.null(d$ref)) return(NULL)
     
    box_end <- end(d$guide) - start(d$guide) - d$seq.width + 1
   
    plotAlignments(
       d$ref,
       alns = NULL,
       target.loc = d$t.loc,
       guide.loc = IRanges(
         start = d$seq.width + 1,
         end = box_end),
       ins.sites = data.frame(),
       axis.text.size = 14,
       plot.text.size = 3,
    )
  })
 
annotation_plot <- reactive({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Annotating transcript location ", value = 0)
 
    increment_prog(progress, 5, "")
  
    if(is.null(d$guide)){    
        return()
    } 
    
    d$txdb <- setTxdb()
     
    CrispRVariants:::annotateGenePlot(txdb = d$txdb, target = d$guide, 
                          gene.text.size = 8)      
  })
 
  output$guide_plot <- renderPlot({
    reference_plot()
  })
  
  output$ref_plot <- renderPlot({
    reference_plot()
  })
 
 output$plot_anot <- renderPlot({
   annotation_plot()
 })
