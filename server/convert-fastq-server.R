observeEvent(input$run_fastq,{
    # File type zip is enforced by input button
    req(input$fq_input)

    closeAlert(session, "alertFASTQ")

    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Preprocessing  FASTQ ", value = 0)
    increment_prog(progress, 15, detail = "unzip FASTQ files", n.inc = 5)

    # Unzip uploaded file, check that fastq files are found
    temp <- unzip(input$fq_input$datapath, exdir = v$fq_dir)
    v$fq_fnames <- dir(v$fq_dir, "fastq$|fq$|fq.gz$|fastq.gz$", recursive = TRUE,
                      full.names = TRUE)
    
    # Warn if no files found
    if (length(v$fq_fnames) == 0){
        createAlert(session, "alertFASTQ", "prepAlertFASTQ",
        title = "WARNING", style = "warning", append = FALSE,
        content = "No files ending in '.fastq' or '.fq' found")
        return(NULL)
    }

    # Map the fastq files to the selected genome
    increment_prog(progress, 15, detail = "Map FastQ reads", n.inc = 6)
    mapFastQ()

    toggleModal(session, "modal_FASTQ", toggle = "close")
    

    createHTable("metadata")
    
    #closeAlert(session, "prepAlertFASTQ")
     
   # }else{
   #   createAlert(session, "alertFASTQ", "prepAlertFASTQ", title = "WARNING",
   #     content = "FASTQ files (.zip) not loaded", style = "warning", append = FALSE)
   # }
    
    
  })
