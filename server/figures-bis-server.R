d <- reactiveValues(
  cset = NULL,
  refs = NULL,
  mds = NULL,
  txdb = NULL
)

setGuides <- reactive({
  inFile <- input$upload_guide
  if(!is.null(inFile)){
    gd <- rtracklayer::import.bed(inFile$datapath) + 5
  }
  return(gd)
})


setTxdb <- reactive({
  f <- paste0("data/txdb/", input$txDb)
  txdb <- AnnotationDbi::loadDb(f)
  return(txdb)
})

setMds <- reactive({
  inFile <- input$upload_metadata
  mds <- read.xls(inFile$datapath)
  # add mds
  d$mds <- split(mds, mds$sgRNA1)
})

setRef <- function(gd){
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating  Reference ", value = 0)
  
  n <- 15
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "A")
    Sys.sleep(0.5)
  }
  
  genome_index <- paste0("/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7/", input$select_Refgenome)
  
  cmd <- paste0("samtools faidx ", genome_index)
  system(cmd)
  cmd <- paste0(cmd, " %s:%s-%s") 
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "B")
    Sys.sleep(0.5)
  }
  
  ref <- system(sprintf(cmd, seqnames(gd)[1], start(gd)[1], end(gd)[1]), intern = T )[[2]]
  
  for (i in 1:5){
    #Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = "C")
    Sys.sleep(0.5)
  }
    return(ref)
}


createPlots <- reactive ({
  if(!is.null(t$DF)){
    mds <- t$DF
    mds <- split(mds, mds$sgRNA1)
    guides <- setGuides()
    
    # list BAM files
    data_dir <- input$upload_bams
    if(is.null(data_dir) && state$bam){
      temp <- unzip(data_dir$datapath, exdir = v$bam_dir)
      v$bm_fnames <-  dir(v$bam_dir, ".bam$", recursive = TRUE, full.names = TRUE)
    }
    for (i in seq_along(mds)){
      md <- mds[[i]]
      bams <- paste0(gsub("[\ |\\/]", "_", md$directory),"_s.bam")
      bam_fnames <- file.path(v$bam_dir, bams)
      refs <- setRef(guides)
      guide <- guides[guides$name == unique(md$sgRNA1)]
      cset <- readsToTarget(bam_fnames, guide, reference = refs[guides$name == unique(md$sgRNA1)], target.loc = input$target_loc, names = as.character(md$Short.name))
      local({
        my_i <- i
        plotname <- paste("plot", my_i, sep="")
        isolate({
          # create the scatter plot
          output[[plotname]] <- renderPlot({
            isolate({
              plotVariants(cset, txdb = setTxdb(), left.plot.margin = grid::unit(c(0.1,0,6,0.2), "lines"))
            })
          })
        })
      })
    }
  }
  
})

# create dynamically the plots
renderUIPlots <- reactive({
  plots <- renderUI({
    norm$plots <- createPlots()
    numberOfPlots <- length(mds)
    plot_output_list <- lapply(1:numberOfPlots, function(i){
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname)
    })
    do.call(tagList, plot_output_list)
  })
  return(plots)
})
