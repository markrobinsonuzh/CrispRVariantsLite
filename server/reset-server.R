   # open metadata pannel
  observeEvent(input$confirm, {
    
    # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Reseting The Session", value = 0)
  
  # Number of times we'll go through the loop
  n <- 15
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "deleting data")
    Sys.sleep(0.005)
  }
  
    toggleModal(session,"modal_reset", toggle = "close")
  
    v$sq_nms = NULL
    v$ab1_fnames = NULL
    v$fq_fnames = NULL
    
    unlink(file.path(v$fq_dir, list.files(v$fq_dir)), recursive= TRUE) 
    unlink(file.path(v$bam_dir, list.files( v$bam_dir)), recursive= TRUE) 
    unlink(file.path(v$ab1_dir, list.files(v$ab1_dir)), recursive= TRUE) 
    
    unlink(v$fq_dir)
    unlink(v$bam_dir)
    unlink(v$ab1_dir)
    
    v$seqs = NULL
    v$bm_fnames = NULL
    v$srt_bm_names = NULL
    d$cset = NULL
    d$ref = NULL
    d$guide = NULL
    state$bam = F
    d$mds = NULL
    d$txdb = NULL
    d$bm_fnames = NULL
    t$DF = NULL
    v$inFile = NULL
    
    state$reset <- TRUE
    
    
    updateTextInput(session, "g.chr", value = NULL)
    updateTextInput(session, "g.start", value = NULL)
    updateTextInput(session, "g.strand", value = NULL)
    
    
    
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "Compiling data")
    Sys.sleep(0.005)
  }
  
  setDir()
  toggleModal(session,"modal_2", toggle = "open")  
 
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.005)
  }
    
      

  })
   
  observeEvent(input$cancel, {
          toggleModal(session,"modal_reset", toggle = "close")
  })