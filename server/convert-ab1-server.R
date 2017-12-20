################################################################################
# Upload AB1, convert and map
# Functions for conversion are defined in preprocessing-server.R
################################################################################

#___________________________________________________________________
# Converting AB1 file to FastQ then bam
observeEvent(input$run_prep,{
    # To do - better error handling if no ab1 files extracted
    
    # Require a zip file to be loaded before running
    req(input$ab1_input)
   
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Preprocessing  AB1 ", value = 0) 
    increment_prog(progress, 15, "Convert .ab1 to FastQ", n.inc = 5)
    
    found_ab1 <- convertAb1toFastq()
    if (! isTRUE(found_ab1)) return()
    
    increment_prog(progress, 15, "Map FASTQ reads", n.inc = 5)
    
    mapFastQ()
    increment_prog(progress, 15, "Compiling data", n.inc = 5) 

    # Close the upload panel and create the metadata table
    closeAlert(session, "alertFASTQ")
    toggleModal(session, "modal_AB1", toggle = "close")
    createHTable("metadata")

 })
