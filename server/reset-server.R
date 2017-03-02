# open metadata pannel
observeEvent(input$confirm, {
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Reseting The Session", value = 0)
  
    increment_prog(progress, 15, n.inc = 5, detail = "deleting data")
  
    toggleModal(session,"modal_reset", toggle = "close")
  
    v$sq_nms = NULL
    v$ab1_fnames = NULL
    v$fq_fnames = NULL
    
    unlink(v$fq_dir, recursive= TRUE)
    unlink(v$bam_dir, recursive= TRUE)
    unlink(v$ab1_dir, recursive= TRUE)
    
    v$bm_fnames = NULL
    v$srt_bm_names = NULL
    d$cset = NULL
    d$ref = NULL
    d$guide = NULL
    d$mds = NULL
    d$txdb = NULL
    d$bm_fnames = NULL
    t$DF = NULL
    v$inFile = NULL
        
    # Note that nulls in updateTextInput are ignored, so values must be set
    updateTextInput(session, "g.chr", value = "")
    updateTextInput(session, "g.start", value = "")
    updateTextInput(session, "ref_seqs", value = "")

    increment_prog(progress, 15, n.inc = 5, detail = "Compiling data")
  
    setDir()
    toggleModal(session,"modal_2", toggle = "open")
 
    increment_prog(progress, 15, n.inc = 5)

})


observeEvent(input$cancel, {
          toggleModal(session,"modal_reset", toggle = "close")
})