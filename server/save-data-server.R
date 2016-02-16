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
  