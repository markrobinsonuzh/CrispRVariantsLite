library(shiny)
library(shinydashboard)


################################################################################
# Options, default settings, and load packages
################################################################################
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 30MB.
options(shiny.maxRequestSize = 10000*1024^2)
# Set Shiny Reaction Log to TRUE
options(shiny.reactlog=TRUE)

# Run the auto-installer/updater code:
# source("install.R", local = TRUE)

# Graphic-saving utilities
source("core/ggsave.R", local = TRUE)


shinyServer(function(input, output, session) {
  
  
  # This code will be run once per user
  users_data <- data.frame(START = Sys.time())
  
  #Reactive values to monitor the state of the app
  state <- reactiveValues(
    ini = FALSE,
    bam = F,
    err = F,
    info = 1
  )
  
  #Reactive values for preprocessing
  v <- reactiveValues(
    sq_nms = NULL, # name of the sequence
    ab1_fnames = NULL, # ab1 files name
    fq_fnames = NULL, #fastq files name
    fq_dir = NULL, # fastq directories
    bam_dir = NULL, # bam directoiry
    ab1_dir = NULL,
    seqs = NULL, # data frame of the data
    bm_fnames = NULL, #bam file
    srt_bm_names = NULL,  #bam directories
    bam_temp = NULL
  )
  
  #Reactive values for producing the plots
  d <- reactiveValues(
    cset = NULL,
    ref = NULL,
    mds = NULL,
    txdb = NULL,
    bm_fnames = NULL,
    guide = NULL
  )

  

  source("server/preprocessing-server.R", local = T)
  source("server/table-server.R", local = T)
  source("server/figures-server.R", local = T)
  source("server/help-tooltip-server.R", local = T)
  source("server/save-data-server.R", local = T)

  # open the modal options
  
  # create the temp dir for the files
  setDir <- reactive({
    temp.dir <- file.path(tempdir(), paste0(MHmakeRandomString(),gsub("[- :]", "", Sys.time())))
    ifelse(!dir.exists(temp.dir), dir.create(temp.dir, showWarnings = FALSE), FALSE)
    
    bam_dir <- file.path(temp.dir, "bam")
    ifelse(!dir.exists(bam_dir), dir.create(bam_dir,  showWarnings = FALSE), FALSE)
    v$bam_dir <- bam_dir
    
    bam_temp <- file.path(temp.dir, "bam_temp")
    ifelse(!dir.exists(bam_temp), dir.create(bam_temp,  showWarnings = FALSE), FALSE)
    v$bam_temp <- bam_temp
    
    
    fq_dir <- file.path(temp.dir, "fastq")
    ifelse(!dir.exists(fq_dir), dir.create(fq_dir,  showWarnings = FALSE), FALSE)
    v$fq_dir <- fq_dir
    
    ab1_dir <- file.path(temp.dir, "ab1")
    ifelse(!dir.exists(ab1_dir), dir.create(ab1_dir,  showWarnings = FALSE), FALSE)
    v$ab1_dir <- ab1_dir
    
  })

  reset <- reactive({
    bam_dir <- file.path(temp.dir, "bam")
    ifelse(dir.exists(bam_dir),unlink(v$bam_dir, recursive = T), FALSE)
    fq_dir <- file.path(temp.dir, "fastq")
    ifelse(dir.exists(fq_dir), unlink(v$fq_dir, recursive = T), FALSE)
    ab1_dir <- file.path(temp.dir, "ab1")
    ifelse(dir.exists(ab1_dir), unlink(v$ab1_dir, recursive = T), FALSE)
    v$sq_nms = NULL
    v$ab1_fnames = NULL
    v$fq_fnames = NULL
    v$fq_dir = NULL
    v$bam_dir = NULL
    v$ab1_dir = NULL
    v$seqs = NULL
    v$bm_fnames = NULL
    v$srt_bm_names = NULL
    d$cset = NULL
    d$ref = NULL
    d$mds = NULL
    d$txdb = NULL
    d$bm_fnames = NULL
    t$DF = NULL
    t$inFile = NULL
  })
  

  
  observe(
    if(!state$ini){
      setDir()
      #file.remove(dir(tempdir(), full.names=TRUE))
      toggleModal(session,"modal_1", toggle = "open")
    }
  )
  
  observeEvent(input$select_AB1,{
      toggleModal(session, "modal_2", toggle = "close")
      toggleModal(session, "modal_AB1", toggle = "open")
      state$ini = TRUE
      state$bam = T
  })
  
  observeEvent(input$select_FastQ,{
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_FASTQ", toggle = "open")
  })
  
  
  observeEvent(input$start_btn,{
    toggleModal(session, "modal_1", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
    state$ini = TRUE
    state$bam = F
  })
  
  observeEvent(input$info_btn,{
    toggleModal(session, "modal_2", toggle = "open")
  })
  
  observeEvent(input$update_xls,{
    toggleModal(session, "modal_table", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
    state$ini = TRUE
    state$bam = F
  })
  
  #---------------------
  # Converting AB1 file to FastQ
  #---------------------
  
 observeEvent(input$run_prep,{
   if(!is.null(input$ab1_files)){
     File <- input$ab1_files
     rfiles <- (tools::file_ext(File$name) == c("zip"))
     
     if(!rfiles){
       createAlert(session, "alertAB1", "prepAlertZIP", title = "WARNING",
         content = "Wrong  Format / upload a zip file", style = "warning", append = FALSE)
     }else{
       closeAlert(session, "prepAlertAB1")
       # Create a Progress object
       progress <- shiny::Progress$new()
       # Make sure it closes when we exit this reactive, even if there's an error
       on.exit(progress$close())
       progress$set(message = "Preprocessing  AB1 ", value = 0)
       
       # Number of times we'll go through the loop
       n <- 15
       
       for (i in 1:5){
         #Increment the progress bar, and update the detail text.
         progress$inc(1/n, detail = "Convert .ab1 to FastQ")
         Sys.sleep(0.5)
       }
       
       convertAb1toFasq()
       
       for (i in 1:5){
         progress$inc(1/n, detail = "Map FastQs reads")
         Sys.sleep(0.5)
       }
       
       mapFastQ()
       
       for (i in 1:5){
         progress$inc(1/n, detail = "compiling data")
         Sys.sleep(0.5)
       }
       
       state$ini = TRUE
       toggleModal(session, "modal_AB1", toggle = "close")
       toggleModal(session, "modal_2", toggle = "open")
       
       createHTable()
       
     }
          }else{
      createAlert(session, "alertAB1", "prepAlertAB1", title = "WARNING",
        content = "AB1 files (.zip) not loaded", style = "warning", append = FALSE)
    }
  })
  
  
  
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
  
  
  
  
  observeEvent(input$select_data, {
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_1", toggle = "open")
  })
  
  observeEvent(input$create_guides, {
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_ref", toggle = "open")
  })
  
  observeEvent(input$edit_xls, {
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_table", toggle = "open")
  })
  
  observeEvent(input$run_guide,{
    creatPlotRef()
    toggleModal(session, "modal_ref", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
  })
  
  observeEvent(input$info_btn, {
    toggleModal(session, "modal_help", toggle = "open")
  })
  
  
 
    
    output$bams <- renderUI({
        if(!is.null(v$bm_fnames))
        {
            tags$div(
                p("BAM files stored"),
                downloadButton('downloadBAM', 'Download')
            )
        }
        else
        {
            fileInput('upload_bams', 'Upload Bams', multiple = F, width = "100%")
        }
    })
 
  
  
  observeEvent(input$reset,{
    reset()
    setDir()
  })
  
  observeEvent(input$run_plot,{
    
    d$cset <- createCripSet()
    d$txdb <- setTxdb()
    createCrispPlot()
    
    toggleModal(session, "modal_2", toggle = "close")
  })
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadBAM <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("BAM","zip", sep = ".")
    },
    
    # the argument 'file'.
    content = function(file){
      # Write to a file specified by the 'file' argument
        fs <- dir(v$bam_dir, pattern = ".bam$", full.names = TRUE, recursive = TRUE)
        zip(zipfile = file, files = fs)
    },
     contentType = "application/zip"
  )

  
 # This code will be run after the client has disconnected
  session$onSessionEnded(function(){
    users_data$END <- Sys.time()
  })
  
})