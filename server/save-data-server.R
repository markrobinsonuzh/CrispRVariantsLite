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
  
  
  observe({
    output$ui_save <- renderUI({
    if(!is.null(d$cset)){
      tags$div(
        downloadButton('downloadPlot', 'Download') 
      )
    }else{
      p("No Data to download - to create plot click on setting")
    }
    })
  })
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadPlot <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
   filename = function() { paste(simpletime(), '.pdf', sep='') },
    content = function(file) {
        ggsave(
            file, 
            plot = createCrispPlot(), 
            device = "pdf",
            width=8, 
            height=11, 
            dpi=100
            )
    }
  )
  