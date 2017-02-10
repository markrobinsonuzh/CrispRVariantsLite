################################################################################
# FUNCTIONS
################################################################################


output$ab1 <- renderUI({
    fileInput(v$ab1_input, 'Upload ZIP file with .AB1 files in directories',
              multiple = F, width = "100%")
})

 #---------------------
 # Converting AB1 file to FastQ
 #---------------------
   
 observeEvent(input$run_prep,{
   #check if file is compressed
   if(!is.null(input[[v$ab1_input]])){
     File <- input[[v$ab1_input]]
     rfiles <- (tools::file_ext(File$name) == c("zip"))
     
     if(!rfiles){
       createAlert(session, "alertAB1", "prepAlertZIP", title = "WARNING",
         content = "Wrong  Format / upload a zip file", style = "warning", append = FALSE)
     }else{
      # closeAlert(session, "prepAlertAB1")
      
       # Create a Progress object
       progress <- shiny::Progress$new()
       # Make sure it closes when we exit this reactive, even if there's an error
       on.exit(progress$close())
       progress$set(message = "Preprocessing  AB1 ", value = 0)
       
       # Number of times we'll go through the loop
       increment_prog(progress, 15, "Convert .ab1 to FastQ", n.inc = 5)
       #n <- 15
       # 
       #for (i in 1:5){
       #  #Increment the progress bar, and update the detail text.
       #  progress$inc(1/n, detail = "Convert .ab1 to FastQ")
       #  Sys.sleep(0.05)
       #}
       
       convertAb1toFasq()
       
       increment_prog(progress, 15, "Map FASTQ reads", n.inc = 5)       

       #for (i in 1:5){
       #  progress$inc(1/n, detail = "Map FastQs reads")
       #  Sys.sleep(0.05)
       #}
       
       mapFastQ()
       
       increment_prog(progress, 15, "Compiling data", n.inc = 5) 
       #for (i in 1:5){
       #  progress$inc(1/n, detail = "compiling data")
       #  Sys.sleep(0.05)
       #}
       
       # Close the upload panel and create the metadata table
       toggleModal(session, "modal_AB1", toggle = "close")
       state$reset <- FALSE
       print(sprintf("d$id %s #1", d$id))
       
       if(!is.null(d$id)){
        createHTable(d$id)    
       } 
       
     }
   }else{
      createAlert(session, "alertAB1", "prepAlertAB1", title = "WARNING",
        content = "AB1 files (.zip) not loaded", style = "warning", append = FALSE)
   }
 })
