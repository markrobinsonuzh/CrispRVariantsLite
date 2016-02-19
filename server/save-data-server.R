   output$bams <- renderUI({
        if(!is.null(v$bm_fnames))
        {
            tags$div(
                p("BAM files stored")
            )
        }
        else
        {
            tags$div(
                fileInput('upload_bams', 'Upload Bams', multiple = F, width = "100%"),
                bsButton("select_FastQ", 'upload FastQ', style = "primary", block = TRUE),
                bsButton("select_AB1", 'upload AB1', style = "primary", block = TRUE),
                p()

            )
        }
    })
  
  observe({
    output$ui_save <- renderUI({
    if(!is.null(d$cset)){
      tags$div(
        splitLayout(p("Download plot"),
        downloadButton('downloadPlot', 'Download')),
        p(),
        splitLayout(p("Download BAM "),
        downloadButton('downloadBAM', 'Download')),
        p(),
        splitLayout(p("Download FASQ  "),
        downloadButton('downloadFASTQ', 'Download')),
        p(),
        splitLayout(p("Download Metadata "),
        downloadButton('downloadTable', 'Download')
        )
      )
    }else{
      p("No Data to download - to create plot click on setting")
    }
    })
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
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadFASTQ <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("FASTQ","zip", sep = ".")
    },
    
    # the argument 'file'.
    content = function(file){
      # Write to a file specified by the 'file' argument
        fs <- dir(v$fq_dir, pattern = ".fastq$", full.names = TRUE, recursive = TRUE)
        zip(zipfile = file, files = fs)
    },
     contentType = "application/zip"
  ) 
  
  # downloadHandler() takes two arguments, both functions.
# The content function is passed a filename as an argument, and
#   it should write out data to that filename.
output$downloadTable <- downloadHandler(
  # This function returns a string which tells the client
  # browser what name to use when saving the file.
  filename = function() {
    paste0("metadata", simpletime(), ".csv")
  },
  
  # This function should write data to a file given to it by
  # the argument 'file'.
  content = function(file) {
    # Write to a file specified by the 'file' argument
    write.table(t$DF, file, sep = ",", row.names = FALSE)
  }
)
  
  # downloadHandler() takes two arguments, both functions.
  # The content function is passed a filename as an argument, and
  #   it should write out data to that filename.
  output$downloadPlot <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
   filename = function() { paste(simpletime(), '.pdf', sep='') },
    content = function(file) {
        pdf(file, height = 5, useDingbats = FALSE)
        
        plotVariants(createCripSet(), txdb = d$txdb, 
     # col.wdth.ratio = c(4,2),
     # row.ht.ratio = c(1,6),
      gene.text.size = 8,
      left.plot.margin = ggplot2::unit(c(1,0,3,2), "lines"),
      
      
      
      plotAlignments.args = list(
        legend.cols =5,
        top.n = input$top.n,
        min.freq = input$min.freq,
        min.count = input$min.count,
        target.loc = d$t.loc,
        guide.loc = IRanges(
         start = d$seq.width + 1,
         end = end(d$guide) - start(d$guide) - d$seq.width + 1),
        axis.text.size = input$axis.text.size, 
        ins.size = input$ins.size,
        plot.text.size = input$plot.text.size, 
        legend.symbol.size = input$legend.symbol.size, 
        legend.text.size = input$legend.text.size
        ), 
        
      plotFreqHeatmap.args = list(
        top.n = input$top.n,
        min.freq = input$min.freq,
        min.count = input$min.count,
        type =  input$plot.type,
        x.size = input$x.size, 
        plot.text.size = input$plot.text.size, 
        legend.text.size = input$legend.text.size,
        x.angle = input$x.angle
        )
      )  
        
        dev.off()
    }
  )
  