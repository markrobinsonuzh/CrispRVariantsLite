observeEvent(input$run_fastq,{
    if(!is.null(input$fastq_files)){
      File <- input$fastq_files
      rfiles <- (tools::file_ext(File$name) == c("zip"))
      
      if(!rfiles){
        createAlert(session, "alertFASTQ", "prepAlertZIP", title = "WARNING",
          content = "Wrong  Format / upload a zip file", style = "warning", append = FALSE)
      }else{
        
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Preprocessing  AB1 ", value = 0)
        
        # Number of times we'll go through the loop
        n <- 15
        
        for (i in 1:5){
          #Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = "unzip Fast FASTQs files")
          Sys.sleep(0.5)
        }
        
        uploadFastq()
        
        for (i in 1:3){
          progress$inc(1/n, detail = "Map FastQs reads")
          Sys.sleep(0.5)
        }
        
        mapFastQ()
        
        for (i in 1:3){
          progress$inc(1/n, detail = "Map FastQs reads")
          Sys.sleep(0.5)
        }
        
        state$ini = TRUE
        
        toggleModal(session, "modal_FASTQ", toggle = "close")
        toggleModal(session, "modal_2", toggle = "open")
        
        for (i in 1:4){
          progress$inc(1/n, detail = "Create Metadata")
          Sys.sleep(0.5)
        }
        
        createHTable()
        
      }
      
      closeAlert(session, "prepAlertFASTQ")
     
    }else{
      createAlert(session, "alertFASTQ", "prepAlertFASTQ", title = "WARNING",
        content = "FASTQ files (.zip) not loaded", style = "warning", append = FALSE)
    }
    
    
  })