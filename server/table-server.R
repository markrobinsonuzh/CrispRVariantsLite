################################################################################
# UI
################################################################################
# uncomment lines below if action button is used to commit changes
# values = list()
# setHot = function(x) values[["hot"]] <<- x

# comment lines below if action button is used to commit changes
values = reactiveValues()
setHot = function(x) values[["htable"]] = x


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
    Sys.sleep(0.05)
  }
  
  
  
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "Compiling data")
    Sys.sleep(0.05)
  }
  
  #toggleModal(session, "modal_table", toggle = "open")
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.05)
  }
  
  output$table <- renderUI({
    rHandsontableOutput("htable")
  })
  
  output$metadata <- renderUI({
    bsButton("edit_xls", "metadata", icon =  icon("table"), style = "success", block = TRUE) 
  })
  
  toggleModal(session, "modal_table", toggle = "open")
  
})

################################################################################
# BEHAVIOUR
################################################################################

#downland the bams files on the server
observe({
  # list BAM files
  data_dir <- input$upload_bams
  
  # if the file doesn't exist
  if(!is.null(data_dir) && is.null(v$bm_fnames))
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
      Sys.sleep(0.05)
    }
    
    downloadbm <- file.path(v$bam_dir)
    temp <- unzip(data_dir$datapath, exdir = downloadbm)
    v$bm_fnames <- dir(downloadbm, ".bam$", full.names = TRUE, recursive = T)
    for (i in 1:10){
      #Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail("storing files on the server"))
      Sys.sleep(0.05)
    }
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
    Sys.sleep(0.05)
  }
  bm_fnames <- dir(v$bam_dir,".bam$", full.names = TRUE, recursive=T)  # re-define
  v$bm_fnames <- bm_fnames
  #mapped sequence
  mp_reads <- lapply(bm_fnames, function(file){
      cmd = paste0( "samtools view -F 2052 -c ", file)
      system(cmd, intern = TRUE)
  })
  
  
  #UNmapped sequence
  ump_reads <- lapply(bm_fnames, function(file){
      cmd = paste0( "samtools view -f 2052 -c ", file)
      system(cmd, intern = TRUE)
  })
  
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.05)
  }
  
  l = length(bm_fnames)
  x <- basename(bm_fnames)
  d_fnames <- unlist(lapply(x,function(x) paste0("/bam/", x )))
  lbl <- gsub("_s.bam$","",basename(bm_fnames))
  stub <- gsub("[\ |\\/]", "_",v$bam_dir)
  stub <- gsub("_bam","",stub)
  lbl <- gsub(stub,"",lbl)
  t.DF = data.frame(file.name = d_fnames, 
                    label = lbl, 
                    group = rep(1,l),
                    number.of.sqs = unlist(mp_reads),
                    stringsAsFactors = FALSE )
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n)
    Sys.sleep(0.05)
  }
  return(t.DF)
})



