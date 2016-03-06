################################################################################
# UI
################################################################################

output$error1 <- renderUI({
    tags$div(
        p(paste0("Reference : " , input$g.start, " strand : ", input$g.strand))
    )
})

output$guide <- renderUI({
    tags$div(
        plotOutput("guide_plot", width="100%", height="200px"),
        uiOutput("error1")
    )
    
})

output$next_step <- renderUI({
    if(!is.null(d$ref)){  
       bsButton("next_step", "Next" , type="action", style = "success", block = TRUE)    
    }else{
       bsButton("next_step", "Back" , type="action", style = "default", block = TRUE)
    }
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

# open reference modal
  observeEvent(input$create_guides, {
    # toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_ref", toggle = "open")
  })

 observeEvent(input$next_step,{
      toggleModal(session, "modal_ref", toggle = "close")
     #toggleModal(session, "modal_2", toggle = "open")
 })
 



################################################################################
# FUNCTIONS
################################################################################

setGuidesFromCoords <- function(){
  if(input$run_guide == 0) return()
 
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
}


setTxdb <- reactive({
  #f <- paste0("data/txdb/", input$txDb)
  # Note: requires that the file names of the txdb match those of ref genome and
  # this can be achieved by symbolic links
  f <- paste0("./data/txdb/", gsub(".fa",".sqlite",input$select_Refgenome))
  d$txdb <- AnnotationDbi::loadDb(f)
  return(d$txdb)
})

observeEvent(input$run_guide,{

  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating  Reference ", value = 0)
  
  n <- 20
  
  for (i in 1:3){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail="Prepare for mapping sequence")
    Sys.sleep(0.05)
  }
  
    ref <- NULL
    genome_index <- paste0("./data/genome/", input$select_Refgenome)

    if(nchar(input$g.start) > 0 ){
            
            for (i in 1:5){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="mapping sequence")
             Sys.sleep(0.05)
            }
            
            print(input$ref_seqs)
            gd <- setGuidesFromCoords()
            
            cmd <- paste0("samtools faidx ", genome_index)
            cmd <- paste0(cmd, " %s:%s-%s")
            ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]   
            
             for (i in 1:5){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="finding strand")
             Sys.sleep(0.05)
            }
                
            switch (input$g.strand,
                "-" =  ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
                "+" =  ref <- Biostrings::DNAString(ref)
                )
            
            updateTextInput(session, "ref_seqs", value = toString(ref))
    
    } else if (nchar(input$ref_seqs) >=  10 ) {
        
        updateSliderInput(session, "g.length", value = 0)
        
            for (i in 1:2){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="getting coordinates")
             Sys.sleep(0.05)
            }
                
        idx <- genome_index
        ref <- input$ref_seqs
                     
        fa <- paste0(MHmakeRandomString(),gsub("[- :]", "", Sys.time()),"_reference.fasta")
        fa <- file.path(v$fasta_temp, fa)

        write.fasta(sequences = ref, names = "reference", file.out = fa )  
        
        #Helen implementaion BWA MEM read to shorts
        cmd <- paste("bwa aln %s %s | bwa samse %s - %s | grep -v '^@' | awk -F \"\t\"",
              "'{if ($2 == 0)print $3, $4, length($10), \"+\";",
              "else if ($2 == 16) print $3, $4, length($10), \"-\"}' && rm %s")
        
        # Run the mapping
        print(sprintf(cmd, idx, fa, idx, fa, fa))
        
        
         #bwa fastmap ./data/genome/hg19.fa test.fa | grep EM | awk '{print $5}'
        #cmd <- sprintf("bwa fastmap %s %s | grep EM | awk '{print $5}'", idx, fa)
  
        # Run the mapping
        result <- strsplit(system(sprintf(cmd, idx, fa, idx, fa, fa), intern = TRUE), " ")[[1]]
        
        chr <- result[1]
        sq.start <- as.numeric(result[2])
        sq.length <- as.numeric(result[3])
        strd <- result[4]
       
         for (i in 1:2){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="creating the guide")
             Sys.sleep(0.05)
            }
        
        guide <- GenomicRanges::GRanges(chr,
                   IRanges( sq.start , end = (sq.start + sq.length - 1)), strand = strd)
        
        d$seq.width <- 0
        d$guide <- guide 
        d$t.loc <- as.numeric(input$target_loc) 
                
        for (i in 1:3){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="updating text input")
             Sys.sleep(0.05)
            }
        
         # update the text input
          updateTextInput(session, "g.chr", value = paste(chr))
          updateTextInput(session, "g.start", value = paste(sq.start))
          updateTextInput(session, "g.strand", value = paste(strd))
                
        for (i in 1:3){
             #Increment the progress bar, and update the detail text.
             progress$inc(1/n, detail="finding strand")
             Sys.sleep(0.05)
            }
        
        ref <- Biostrings::DNAString(ref)
                
    } else{
        createAlert(session, "alertRef", alertId = "alertRef1", title = "WARNING",
        content = "Please enter the sequence or set the coordinates", style = "danger", 
        dismiss = TRUE, append = FALSE)
    }
        
 
  
      
  for (i in 1:7){
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.05)
  }
  
  d$ref <- Biostrings::DNAString(ref)
  
})

################################################################################
#  ###############################################################################
 
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
    for (i in 1:7){
      # Increment the progress bar, and update the detail text.
      progress$inc(1/i)
      Sys.sleep(0.05)
    }
  
    if(is.null(d$guide)){
        return(NULL)
    } 
    d$txdb <- setTxdb()
    print("output$plot_anot")
    print(d$guide)
    print(d$txdb)
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
 