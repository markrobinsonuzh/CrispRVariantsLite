################################################################################
# UI
################################################################################

createHTable <- reactive({
  
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating Table ", value = 0)
  
  # Number of times we'll go through the loop
  n <- 15
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.1)
  }
  
  
  output$htable <- renderRHandsontable({
    if (!is.null(input$htable)) {
      t$DF <- hot_to_r(input$htable)
    }else{
      t$DF <- getMetadata()
    }
    setHot(t$DF)
    rhandsontable(t$DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "Compiling data")
    Sys.sleep(0.1)
  }
  
  toggleModal(session, "modal_table", toggle = "open")
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.1)
  }
  
  output$table <- renderUI({
    rHandsontableOutput("htable")
  })
  
  output$metadata <- renderUI({
    bsButton("edit_xls", "metadata", icon =  icon("table"), style = "info", block = TRUE) 
  })
  
  
})

################################################################################
# BEHAVIOUR
################################################################################

observe({
  # list BAM files
  data_dir <- input$upload_bams
  
  # if the file doesn't exist
  if(!is.null(data_dir) && is.null(d$bm_fnames))
  {
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Unzipping the BAM files", value = 0)
    
    # Number of times we'll go through the loop
    n <- 20
    
    for (i in 1:10){
      #Increment the progress bar, and update the detail text.
      progress$inc(1/n)
      Sys.sleep(0.1)
    }
    
    downloadbm <- file.path(v$bam_dir)
    #ifelse(!dir.exists(downloadbm), dir.create(downloadbm,  showWarnings = FALSE), FALSE)
    temp <- unzip(data_dir$datapath, exdir = downloadbm)
    v$bm_fnames <- dir(downloadbm, ".bam$", full.names = TRUE, recursive = T)
    for (i in 1:10){
      #Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail("storing files on the server"))
      Sys.sleep(0.1)
    }
    createHTable()
  }
  
})


################################################################################
# FUNCTIONS
################################################################################

# ----------------------
#  make table of all sequences
# ----------------------
getMetadata <- reactive({
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Reading the metadata", value = 0)
  
  # Number of times we'll go through the loop
  n <- 15
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.01)
  }
  bm_fnames <- dir(v$bam_dir,".bam$", full.names = TRUE, recursive=T)  # re-define
  v$bm_fnames <- bm_fnames
  #mapped sequence
  mp_reads <- lapply(bm_fnames, function(file){
      cmd = paste0( "samtools view -F 0x04 -c ", file)
      system(cmd, intern = TRUE)
  })
  
  
  #UNmapped sequence
  ump_reads <- lapply(bm_fnames, function(file){
      cmd = paste0( "samtools view -f 0x04 -c ", file)
      system(cmd, intern = TRUE)
  })
  
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.01)
  }
  
  l = length(bm_fnames)
  x <- basename(bm_fnames)
  d_fnames <- unlist(lapply(x,function(x) paste0("/bam/", x )))
  t.DF = data.frame(file.name = d_fnames, 
                    label = gsub(".bam$","",basename(bm_fnames)), 
                    group = rep(1,l),
                    number.of.sqs = unlist(mp_reads),
                    remaining.number.of.sqs = unlist(ump_reads)
                    )
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.01)
  }
  return(t.DF)
})


# downloadHandler() takes two arguments, both functions.
# The content function is passed a filename as an argument, and
#   it should write out data to that filename.
output$downloadTable <- downloadHandler(
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    paste("metadata", input$filetype, sep = ".")
  },
  
  # This function should write data to a file given to it by
  # the argument 'file'.
  content = function(file) {
    sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")
    # Write to a file specified by the 'file' argument
    write.table(t$DF, file, sep = sep,
      row.names = FALSE)
  }
)
