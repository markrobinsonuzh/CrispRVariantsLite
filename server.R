library(shiny)
library(shinydashboard)

################################################################################
# Options, default settings, and load packages
################################################################################
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 30MB.
options(shiny.maxRequestSize = 10000*1024^2)
# Set Shiny Reaction Log to TRUE
#options(shiny.reactlog=TRUE)
#options(shiny.error = browser)
# Run the auto-installer/updater code:
# source("install.R", local = TRUE)

shinyServer(function(input, output, session) {

  t <- reactiveValues(
      DF = NULL # The metadata table
  )
  
  #Reactive values for preprocessing
  v <- reactiveValues(
      ab1_fnames = NULL, # ab1 files name
      ab1_dir = NULL,
      fq_fnames = NULL, #fastq file names
      fq_dir = NULL, # fastq directories
      bm_fnames = NULL, #bam file
      bam_dir = NULL, # bam directoiry
      srt_bm_names = NULL,  # sorted bam files
      sq_nms = NULL, # name of the sequence
      fasta_temp = NULL,
      inFile = NULL,
      temp.dir = NULL
  )
  
  #Reactive values for producing the plots
  d <- reactiveValues(
     cset = NULL,
     ref = NULL,
     mds = NULL,
     txdb = NULL,
     bm_fnames = NULL,
     guide = NULL,
     seq.width = NULL,
     t.loc = NULL,
     id = NULL,
     use.coords = NULL
  )

  # Create the temp directories for the files
  setDir <- function(){
      v$temp.dir <- file.path(tempdir(), paste0(MHmakeRandomString(),
                              gsub("[- :]", "", Sys.time())))
      
      v$bam_dir <- file.path(v$temp.dir, "bam")
      v$fasta_temp <- file.path(v$temp.dir, "fasta_temp")
      v$fq_dir <- file.path(v$temp.dir, "fastq")
      v$ab1_dir <- file.path(v$temp.dir, "ab1")
      
      dummy <- lapply(c(v$temp.dir, v$bam_dir, v$fasta_temp,
                         v$fq_dir, v$ab1_dir), create_dir)

  }

  cat(sprintf("start sourcing %s \n", simpletime()), file = "crv_scratch.txt", append = TRUE)

  #  source("server/warning-server.R", local = T)
  source("server/preprocessing-server.R", local = T)
  source("server/convert-ab1-server.R", local = T)
  source("server/convert-fastq-server.R", local = T)
  source("server/table-server.R", local = T)
  source("server/guides-server.R", local = T)
  source("server/figures-server.R", local = T)
  source("server/help-tooltip-server.R", local = T)
  source("server/save-data-server.R", local = T)
  source("server/reset-server.R", local =T)
  cat(sprintf("finished sourcing %s \n", simpletime()), file = "crv_scratch.txt", append = TRUE)

  
  # Open the modal options.  modal_1 is the data selection modal
  observe({
      setDir()
      toggleModal(session,"modal_1", toggle = "open")
  })
  
  # Open/close the AB1 modal
  observeEvent(input$select_AB1,{
      toggleModal(session, "modal_AB1", toggle = "open")
  })
  
  observeEvent(input$back_ab1,{
    toggleModal(session, "modal_AB1", toggle = "close")
  })
  
  # Open/close the FASTQ modal
  observeEvent(input$select_FastQ,{
    toggleModal(session, "modal_FASTQ", toggle = "open")
  })
  
  observeEvent(input$back_fastq,{
    toggleModal(session, "modal_FASTQ", toggle = "close")
  })
  
  # Open data input screen, close welcome screen
  observeEvent(input$start_btn,{
    toggleModal(session, "modal_1", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
  })
  
  # Open instructions modal
  observeEvent(input$select_data, {
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_1", toggle = "open")
  })

  # Open metadata pannel
  observeEvent(input$edit_xls, {
    toggleModal(session, "modal_table", toggle = "open")
  })

})

################################################################################
# UTILITY FUNCTIONS
################################################################################


create_dir <- function(dir.name){
    if (!dir.exists(dir.name)) { dir.create(dir.name, showWarnings = FALSE) }
}

increment_prog <- function(prog, n, message = NULL, time.int = 0.05,
                           n.inc = n, detail = NULL){
   for (i in 1:n.inc){
    #Increment the progress bar, and update the detail text.
    prog$inc(1/n, message = message, detail = detail)
    Sys.sleep(time.int)
  } 
}


