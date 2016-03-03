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
                bsButton("select_AB1", 'Upload ZIP of AB1s', style = "primary", block = TRUE),
                bsButton("select_FastQ", 'Upload ZIP of FASTQs', style = "primary", block = TRUE),
                p(),
                fileInput('upload_bams', 'Upload ZIP of BAMs', multiple = F, width = "100%")

            )
        }
    })
  
  observe({
    output$ui_save <- renderUI({
    if(!is.null(d$cset)){
      if(!is.null(v$fq_fnames)){
           tags$div(
        
        splitLayout(p("Allele Variant Plot"),
        downloadButton('downloadPlot', 'Download PDF')),
        p(),
        splitLayout(p("Mapped BAM files"),
        downloadButton('downloadBAM', 'Download ZIP')),
        p(),
        splitLayout(p("FASTQ files"),
        downloadButton('downloadFASTQ', 'Download ZIP')),
        p(),
        splitLayout(p("Metadata table"),
        downloadButton('downloadTable', 'Download CSV')
        )
      ) 
        }else{
           tags$div(
        
        splitLayout(p("Allele Variant Plot"),
        downloadButton('downloadPlot', 'Download PDF')),
        p(),
        splitLayout(p("Mapped BAM files"),
        downloadButton('downloadBAM', 'Download ZIP')),
        p(),
        splitLayout(p("Metadata table"),
        downloadButton('downloadTable', 'Download CSV')
        )
      )  
        }
      
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
        
        pixelratio <- session$clientData$pixelratio
        
        cdata <- session$clientData
        cnames <- names(cdata)
        br_ht <- cdata[["output_crispplots_height"]]
        br_wd <- cdata[["output_crispplots_width"]]
          
        br_ht <- 600 # MAKE THESE GLOBAL VARIABLES
        res <- 72
        ht <- (br_ht*pixelratio)/res
        wd <- (br_wd*pixelratio)/res
                
        pdf(file, height = ht, width = wd, useDingbats = FALSE)
        
        group <- as.factor(t$DF$group)
        
         plotVariants(d$cset, txdb = d$txdb, 
#      col.wdth.ratio = c(5,2),
#      row.ht.ratio = c(1,6),
      gene.text.size = 8,
      left.plot.margin = ggplot2::unit(c(1,0,5,2), "lines"),
      
      plotAlignments.args = list(
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
        x.angle = input$x.angle,
        group = group
        )
      ) + theme(legend.position="none")  
    
        dev.off()
    }
  )
  
