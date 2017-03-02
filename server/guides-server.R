################################################################################
# UI
################################################################################

# Box between plots showing coordinates
output$error1 <- renderUI({
    tags$div(
        p(sprintf("Coordinates: %s  start: %s  strand: %s",
                  input$g.chr, input$g.start, input$g.strand))
    )
})

# Guide plot collects the plot and the coordinates (stored as error1)
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

# Disable creating the guides until bams defined
# This button appears on the figure setting screen
# (either remove or add a "Settings" button from guide screen)
observe({
   if(is.null(v$bm_fnames)){
      updateButton(session,"create_guides", style ="default", icon = icon("ban"),
                   disable = TRUE)
    }else{
      updateButton(session,"create_guides", style = "success",
                   icon = icon("area-chart"), block = TRUE, disable = FALSE)
    }
})

# Open reference modal
observeEvent(input$create_guides, {
    toggleModal(session, "modal_ref", toggle = "open")
})

# Close reference modal when finished creating guide
observeEvent(input$next_step,{
      toggleModal(session, "modal_ref", toggle = "close")
})

# Switch indicating whether sequence or coordinates were last edited
observeEvent({
    input$g.start
    input$g.chr
    input$g.strand
    input$g.length
   }, {d$use.coords = TRUE})


observeEvent({
   input$ref.seqs
   }, {d$use.coords = FALSE})

################################################################################
# FUNCTIONS
################################################################################


# Set genome and transcription data base
setTxdb <- reactive({
  cat(sprintf("Setting txdb to %s \n", input$select_Refgenome),
      file = "crv_scratch.txt", append = TRUE)

  # Note: requires that the file names of the txdb match those of ref genome and
  # this can be achieved by symbolic links
  f <- paste0("./data/txdb/", gsub(".fa",".sqlite", input$select_Refgenome))
  AnnotationDbi::loadDb(f)
})


genome_index <- reactive({
    paste0("./data/genome/", input$select_Refgenome)
})


observeEvent(input$select_Refgenome, {
    chrs <- read.table(paste0(genome_index(), ".fai"),
                   stringsAsFactors = FALSE)[,1]
    
    updateSelectInput(session, "g.chr", choices = chrs, selected = chrs[1])
    
    cat("New chromosome options:\n")
    cat(chrs)
    cat("\n")

})


setGuideCoords <- function(window = NULL, seq.start = NULL, target.loc = NULL,
                           strand = NULL, chromosome = NULL){
  # Store guide coordinates including window around guide

  cat("Setting guide coordinates \n", file = "crv_scratch.txt", append = TRUE)
  
  seq.start <- ifelse(is.null(seq.start), as.numeric(input$g.start), seq.start)
  target.loc <- ifelse(is.null(target.loc), as.numeric(input$target_loc), target.loc)
  chromosome <- ifelse(is.null(chromosome), input$g.chr, chromosome)
  strand <- ifelse(is.null(strand), input$g.strand, strand)
  
  seq.end <-  seq.start + target.loc + 5
  
  ncat <- function(x){
    cat(sprintf("%s\n", x))
  }
  
  ncat(seq.start)
  ncat(target.loc)
  ncat(seq.end)
  ncat(input$g.chr)
  ncat(input$g.strand)
  ncat(input$g.length)
  
  guide <- GenomicRanges::GRanges(
      seqnames = chromosome,
      ranges = IRanges::IRanges(start = seq.start, end = seq.end),
      strand = strand
  )

  # Add window around guide, adjust target location 
  d$seq.width <- ifelse(is.null(window), as.numeric(input$g.length), window)
  d$guide <- guide + d$seq.width
  d$t.loc <- target.loc + d$seq.width
  
  ncat(d$seq.width)
  ncat(d$guide)
  ncat(d$t.loc)
  
  return(d$guide)
}


setGuidesFromCoords <- function(genome_idx, progress){
    # This function gets the sequence matching the specified coordinates
    
    cat("Setting guide from coords \n", file = "crv_scratch.txt", append = TRUE)
 
    # Set the guide and fetch the guide sequence directly from the reference
 
    increment_prog(progress, 5, "Fetching sequence")
    
    gd <- setGuideCoords()
    cmd <- "samtools faidx %s %s:%s-%s"
    
    # Try to extract the sequence matching this location from the reference
    ref <- system(sprintf(cmd, genome_idx, seqnames(gd)[1], start(gd)[1], end(gd)[1]),
                  intern = TRUE)
 
    if(length(ref) < 2){
        cat("Faidx error \n", file = "crv_scratch.txt", append = TRUE)

        createAlert(session, "alertRef", alertId = "alertRef3", title = "INFO",
            content = sprintf("No sequence found at %s:%s in genome %s",
                              input$g.chr, input$g.start, input$select_Refgenome),
             style = "info", dismiss = TRUE, append = FALSE)
             return(NULL)
    }
    
    ref <- paste0(ref[2:length(ref)], collapse = "")
    
    # Reverse complement the sequence if a negative strand was specified
    switch (input$g.strand,
         "-" =  ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
         "+" =  ref <- Biostrings::DNAString(ref))
  
   updateTextInput(session, "ref_seqs", value = toString(ref)) 
   return(ref)
}


writeGuideFasta <- function(ref){
    cat("Writing guide fasta \n", file = "crv_scratch.txt", append = TRUE)

    fa <- paste0(MHmakeRandomString(), 
                 gsub("[- :]", "", Sys.time()), "_reference.fasta")
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
    
    cat("Run guide clicked \n", file = "crv_scratch.txt", append = TRUE)

    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Creating  Reference ", value = 0) 


    # If coordinates exist and were edited most recently than sequence
    if (isTRUE(d$use.coords) & nchar(input$g.start) > 0){
        ref <- setGuidesFromCoords(genome_index(), progress)
    
    } else if (nchar(input$ref_seqs) >=  10) {
    
        # Check bwa is available
        if (system("bwa", ignore.stderr = TRUE) == 127){
            createAlert(session, "alertRef","alertRef5", title = "WARNING",
            content = "Please make sure BWA is available", 
            style = "warning",
            append = FALSE, dismiss = FALSE)
            stopApp()
        }
     
        increment_prog(progress, 3, "Prepare for mapping sequence")

        # This update to slider doesn't seem to be propagated until
        # return from function, possibly because of tabs
        updateSliderInput(session, "g.length", value = 0)
        increment_prog(progress, 2, "getting coordinates")

        # Map input sequence and parse the results
        idx <- genome_index()
        ref <- input$ref_seqs
        print(ref)

        result <- mapGuide(ref, idx)
        if (is.null(result)) return()

        chr <- result[1]
        sq.start <- as.numeric(result[2])
        sq.length <- as.numeric(result[3])
        strd <- result[4]
      
        # Update the text input
        increment_prog(progress, 3, "updating text input")

        updateSelectInput(session, "g.chr", selected = paste(chr))
        updateTextInput(session, "g.start", value = paste(sq.start))
        updateTextInput(session, "g.strand", value = paste(strd))
          
        increment_prog(progress, 2, "creating the guide") 
        setGuideCoords(seq.start = sq.start, chromosome = chr, strand = strd,
                       target.loc = (nchar(input$ref_seqs) - 6), window = 0)
      
    } else{
        createAlert(session, "alertRef", alertId = "alertRef1", title = "WARNING",
        content = "Please enter the sequence or set the coordinates", style = "info",
        dismiss = TRUE, append = FALSE)
        
        return()
    }
    
    d$use.coords <- NULL
    #increment_prog(progress, 7, "")
  
    cat(sprintf("Setting reference from %s \n", isolate(as.character(d$ref))),
        file = "crv_scratch.txt", append = TRUE)
  
    d$ref <- as.character(ref)

    cat(sprintf("Reference set to %s\n", ref ) , file = "crv_scratch.txt", append = TRUE)
  
    reference_plot()
  
})


################################################################################
# PLOTS
################################################################################


observeEvent(input$run_plot_guide,{
    #if (input$run_guide == 0) return(NULL)

    cat("Making allele plot \n", file = "crv_scratch.txt", append = TRUE)

    d$cset <- createCrisprSet()
    setTxdb()
    createCrispPlot()
    toggleModal(session, "modal_ref", toggle = "close")
    Sys.sleep(0.05)
    toggleModal(session, "modal_2", toggle = "close")
  })




#_________________________________________________________________________________ 
# Plot the transcripts
output$plot_anot <- renderPlot({

    cat("rendering annotation plot \n", file = "crv_scratch.txt", append = TRUE)

    annotation_plot()

})


annotation_plot <- reactive({
    req(input$run_guide, d$guide)
    
    cat("making annotation plot \n", file = "crv_scratch.txt", append = TRUE)

    cat(sprintf("Txdb guide loc: %s - %s, width %s\n",
        start(d$guide), end(d$guide), d$seq.width),
        file = "crv_scratch.txt", append = TRUE)

    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Annotating transcript location ", value = 0)
 
    increment_prog(progress, 5, "")
     
    CrispRVariants:::annotateGenePlot(txdb = setTxdb(), target = d$guide,
                          gene.text.size = 8)
})

#_________________________________________________________________________________ 
# Plot the reference sequence (appears on two modals)
reference_plot <- reactive({
    req(input$run_guide, d$ref, d$guide, d$seq.width, d$t.loc)

    cat("making reference plot \n", file = "crv_scratch.txt", append = TRUE)

    cat(sprintf("Guide loc: %s - %s, width %s\n",
        start(d$guide), end(d$guide), d$seq.width),
        file = "crv_scratch.txt", append = TRUE)
    box_end <- end(d$guide) - start(d$guide) - d$seq.width + 1
    
    cat(sprintf("Ref: %s\n", as.character(d$ref)),
        file = "crv_scratch.txt", append = TRUE)

    p <- CrispRVariants::plotAlignments(
       Biostrings::DNAString(d$ref),
       alns = NULL,
       target.loc = d$t.loc,
       guide.loc = IRanges(
         start = d$seq.width + 1,
         end = box_end),
       ins.sites = data.frame(),
       axis.text.size = 14,
       plot.text.size = 3,
    )
    return(p)
})


output$guide_plot <- renderPlot({
      cat("rendering guide plot \n", file = "crv_scratch.txt", append = TRUE)
      reference_plot()
})

output$ref_plot <- renderPlot({
      cat("rendering reference plot \n", file = "crv_scratch.txt", append = TRUE)
      reference_plot()
})



