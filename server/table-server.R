################################################################################
# Create the metadata table object
################################################################################


values = reactiveValues()
val_ = NULL

table <- function(name){
   
    y <- paste0(name) 
    output[[y]] <- renderRHandsontable({
    #solution from : https://github.com/jrowen/rhandsontable/issues/27
       if (is.null(input[[y]]) || val_ != y ) 
       {   
           val_ <<- y
           values[[y]] <- getMetadata()
        }
        else if(!is.null(input[[y]]) && val_ == y)
        {
            values[[y]] <- hot_to_r(input[[y]])
        } 
    
    t$DF <- values[[y]]
    rhandsontable(t$DF) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
})
    
}

output$table <- renderUI({
     rHandsontableOutput("metadata")
})


createHTable <- function(id){  
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Creating Table ", value = 0)
  
  increment_prog(progress, 15, "", n.inc = 5)
  table(id)

  increment_prog(progress, 15, "Compiling data", n.inc = 10)

  toggleModal(session, "modal_table", toggle = "open")  
}

################################################################################
# BEHAVIOUR
################################################################################

# Disable creating the guides until Reference and Bams defined
observe({
   if(is.null(v$bm_fnames)){
      updateButton(session,"edit_xls", style = "default", icon = icon("ban"),
                   disable = TRUE)
    } else {
      updateButton(session,"edit_xls", style = "success",  icon = icon("table"),
                  block = TRUE, disable = FALSE ) 
    }
})

# close metadata panel
  observeEvent(input$guide_from_table, {
    toggleModal(session, "modal_table", toggle = "close")
    toggleModal(session, "modal_ref", toggle = "open")
  })

#downland the bams files on the server
observeEvent(input$upload_bams, {
   v$inFile <- input$upload_bams
    
  # if the file doesn't exist
  if(!is.null(v$inFile))
  {
    # list BAM files
    data_dir <- v$inFile
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Unzipping the BAM files", value = 0)
   
    increment_prog(progress, 20, "", n.inc = 10) 
    
    downloadbm <- file.path(v$bam_dir)
    temp <- unzip(data_dir$datapath, exdir = downloadbm)
    temp <- dir(downloadbm, ".bam$", full.names = TRUE, recursive = T)
    
    base_nms <- gsub(sprintf("%s[\\/]*", downloadbm), "", temp)
    safe_bm_fnames <- gsub("[^[:alnum:]\\/\\.]", "_",
                           iconv(base_nms, to = "ASCII//TRANSLIT"))
    safe_bm_fnames <- gsub("_{2,}", "_", safe_bm_fnames)
    safe_bm_fnames <- file.path(downloadbm, safe_bm_fnames)
    file.rename(temp, safe_bm_fnames)
    v$bm_fnames <- safe_bm_fnames
    
    increment_prog(progress, 20, "Storing files on the server", n.inc = 10)
    createHTable("metadata")
  
  }
})




################################################################################
# FUNCTIONS
################################################################################

# ----------------------
#  make table of all sequences
# ----------------------
getMetadata <- function(){

  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Reading the metadata", value = 0)
  
  increment_prog(progress, 15, "", n.inc = 5)

  bm_fnames <- dir(v$bam_dir,".bam$", full.names = TRUE, recursive=T)  # re-define
  v$bm_fnames <- bm_fnames
  
  #Count mapped sequences
  mp_reads <- lapply(bm_fnames, function(file){
      cmd = paste0( "samtools view -F 2052 -c ", file)
      system(cmd, intern = TRUE)
  })
  
  #print("counted mapped seqs")
  increment_prog(progress, 15, "", n.inc = 5)
  
  l = length(bm_fnames)
  x <- basename(bm_fnames)
  #d_fnames <- unlist(lapply(x,function(x) paste0("/bam/", x )))
  lbl <- gsub("_s.bam$","",x)
  lbl <- sapply(strsplit(lbl,"__"), tail, 1)
  t.DF = data.frame(file.name = x, 
                    label = lbl, 
                    group = rep(1,l),
                    number.of.sqs = unlist(mp_reads),
                    stringsAsFactors = FALSE )
  
  return(t.DF)
}
