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
#options(shiny.error = browser)
# Run the auto-installer/updater code:
# source("install.R", local = TRUE)

shinyServer(function(input, output, session) {

# uncomment lines below if action button is used to commit changes
# values = list()
# setHot = function(x) values[["hot"]] <<- x

# comment lines below if action button is used to commit changes
#values = reactiveValues()
#setHot = function(x) values[["htable"]] = x

t <- reactiveValues(
    DF = NULL,
    inFile = NULL
)
<<<<<<< HEAD

ID_ = NULL
=======
>>>>>>> 629f251f2facd2f7fa86cb82d084cf72571c6ca5
  
  # This code will be run once per user
  users_data <- data.frame(START = Sys.time())
 
  #Reactive values to monitor the state of the app
  state <- reactiveValues(
    ini = FALSE,
    bam = F,
    err = F,
    info = 1,
    reset = FALSE
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
    fasta_temp = NULL,
<<<<<<< HEAD
    inFile = NULL,
    ab1_input = NULL
=======
    inFile = NULL
>>>>>>> 629f251f2facd2f7fa86cb82d084cf72571c6ca5
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
<<<<<<< HEAD
   plot = NULL,
   id = NULL
=======
   plot = NULL
>>>>>>> 629f251f2facd2f7fa86cb82d084cf72571c6ca5
 )

   # create the temp dir for the files
  setDir <- function(){
    d$id <- MHmakeRandomString()
    v$ab1_input <- paste0("ab1_files_",d$id)    

    print(sprintf("d$id %s #0", d$id))

    
    temp.dir <- file.path(tempdir(), paste0(MHmakeRandomString(),gsub("[- :]", "", Sys.time())))
    ifelse(!dir.exists(temp.dir), dir.create(temp.dir, showWarnings = FALSE), FALSE)
    
    bam_dir <- file.path(temp.dir, "bam")
    ifelse(!dir.exists(bam_dir), dir.create(bam_dir,  showWarnings = FALSE), FALSE)
    v$bam_dir <- bam_dir
    
    fasta_temp <- file.path(temp.dir, "fasta_temp")
    ifelse(!dir.exists(fasta_temp), dir.create(fasta_temp,  showWarnings = FALSE), FALSE)
    v$fasta_temp <- fasta_temp
    
    
    fq_dir <- file.path(temp.dir, "fastq")
    ifelse(!dir.exists(fq_dir), dir.create(fq_dir,  showWarnings = FALSE), FALSE)
    v$fq_dir <- fq_dir
    
    ab1_dir <- file.path(temp.dir, "ab1")
    ifelse(!dir.exists(ab1_dir), dir.create(ab1_dir,  showWarnings = FALSE), FALSE)
    v$ab1_dir <- ab1_dir
    
  }


<<<<<<< HEAD
=======

>>>>>>> 629f251f2facd2f7fa86cb82d084cf72571c6ca5

  source("server/preprocessing-server.R", local = T)
  source("server/convert-ab1-server.R", local = T)
  source("server/convert-fastq-server.R", local = T)
  source("server/table-server.R", local = T)
  source("server/guides-server.R", local = T)
  source("server/figures-server.R", local = T)
  source("server/help-tooltip-server.R", local = T)
  source("server/save-data-server.R", local = T)
  source("server/reset-server.R", local =T)
  
  # open the modal options
  
  observe(
    if(!state$ini){
      setDir()
      #file.remove(dir(tempdir(), full.names=TRUE))
      toggleModal(session,"modal_1", toggle = "open")
    }
  )
  
  # open the AB1 modal
  observeEvent(input$select_AB1,{
      #toggleModal(session, "modal_2", toggle = "close")
      toggleModal(session, "modal_AB1", toggle = "open")
      state$ini = TRUE
  })
  
  observeEvent(input$back_ab1,{
    #toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_AB1", toggle = "close")
  })
  
  # open the FASTQ modal
  observeEvent(input$select_FastQ,{
    #toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_FASTQ", toggle = "open")
  })
  
  observeEvent(input$back_fastq,{
    #toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_FASTQ", toggle = "close")
  })
  
  # open start modal 
  observeEvent(input$start_btn,{
    toggleModal(session, "modal_1", toggle = "close")
    toggleModal(session, "modal_2", toggle = "open")
    state$ini = TRUE
  })
  
  # open instructions modal
  observeEvent(input$select_data, {
    toggleModal(session, "modal_2", toggle = "close")
    toggleModal(session, "modal_1", toggle = "open")
  })

  # open metadata pannel
  observeEvent(input$edit_xls, {
    toggleModal(session, "modal_table", toggle = "open")
<<<<<<< HEAD
    })


=======
  })
  
  
   
  
  
  
  
  
   # open metadata pannel
  observeEvent(input$reset, {
      
    v$sq_nms = NULL
    v$ab1_fnames = NULL
    v$fq_fnames = NULL
    #v$fq_dir = NULL
    #v$bam_dir = NULL
    #v$ab1_dir = NULL
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
      
      
      #empty_dir(v$bam_dir)
      #empty_dir(v$fq_dir)
      #empty_dir(v$ab1_dir)
      #setDir()
      updateTextInput(session, "g.chr", value = NULL)
        updateTextInput(session, "g.start", value = NULL)
        updateTextInput(session, "g.strand", value = NULL)
      toggleModal(session,"modal_1", toggle = "open")
  })
    
>>>>>>> 629f251f2facd2f7fa86cb82d084cf72571c6ca5
})
