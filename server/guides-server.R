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
    if(!is.null(d$ref)){
       bsButton("next_step", "Next step" , type="action", style = "success", block = TRUE)    
    }
    #actionButton("next_step", "Next step" , width="100%", style="primary")
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
   print("setRef")
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating  Reference ", value = 0)
  
  n <- 15
  
  for (i in 1:8){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail="Prepare for mapping sequence")
    Sys.sleep(0.05)
  }
  
  
    ref <- NULL
    genome_index <- paste0("./data/genome/", input$select_Refgenome)
    
    if(nchar(input$g.start) > 2 && nchar(input$g.chr) > 2){
            gd <- setGuides()
            
            cmd <- paste0("samtools faidx ", genome_index)
            cmd <- paste0(cmd, " %s:%s-%s")
            ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]   
    } 
    

    
     if(nchar(input$ref_seqs) >= 23){
        idx <- genome_index
        ref <- input$ref_seqs
        
        print(sprintf("this is the index genome : %s ",idx))
        print(sprintf("this is the reference sequence from the user : %s ", ref))
        
        #fa <- paste0(MHmakeRandomString(),gsub("[- :]", "", Sys.time()),"_reference.fasta")
        #fa <- file.path(v$fasta_temp, fa)
        fa <- sprintf("%s.fasta", gsub(" ", "", date()))
        # Write the fasta
        cat(sprintf(">ref\n%s\n", ref), file=fa)
        #write.fasta(sequences = ref, names = "reference", file.out = fa ) 
        
        #Helen implementaion BWA MEM read to shorts
        cmd <- paste("bwa aln %s %s | bwa samse %s - %s |  grep -v '^@' | awk -F \"\t\"",
              "'{if ($2 == 0)print $3, $4, length($10), \"+\";",
              "else if ($2 == 16) print $3, $4, length($10), \"-\"}' && rm %s")
        
        # Run the mapping
        result <- strsplit(system(sprintf(cmd, idx, fa, idx, fa, fa), intern = TRUE), " ")[[1]]
        
        print(sprintf("result %s", result))
        
        # update the text input
        updateTextInput(session, "g.chr", value = result[[1]])
        updateTextInput(session, "g.start", value = result[[2]])
        updateTextInput(session, "g.strand", value = result[[4]])
        updateSliderInput(session, "g.length", value = 17)
        
        
        
  
        gd <- GenomicRanges::GRanges(result[1], IRanges(as.numeric(result[2]), 
            width = as.numeric(result[3])), strand = result[4])
        
    }
        
  
  
  switch (input$g.strand,
    "-" =  ref <- Biostrings::reverseComplement(Biostrings::DNAString(ref)),
    "+" =  ref <- Biostrings::DNAString(ref)
  )
      
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
 output$guide_plot <- renderPlot({
   plot_reference()
    
  }, height = 200)
  
  output$ref_plot <- renderPlot({
    plot_reference()
 }, height = 200)

})


observeEvent(input$run_guide,{
        print("create guide")
       creatPlotRef()  
 })
 
 observeEvent(input$next_step,{
        toggleModal(session, "modal_ref", toggle = "close")
        toggleModal(session, "modal_2", toggle = "open")
 })
 
plot_reference <- reactive({

   print("plot_reference")     
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

