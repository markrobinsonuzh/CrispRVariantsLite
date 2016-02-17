info_AB1 = "Convert AB1 format to FastQ and run the mapping"
info_FastQ = "Map the FastQ against a selected genome "
info_guides = "Set the target and reference sequence"
info_xls = "Set the metadata"
observe({
  if(input$info){
    #--------
    # AB1 information
    #--------
    
    addTooltip(session, "select_AB1", info_AB1, placement = "left", trigger = 'hover', options = list(container = "body"))
    
    #--------
    # FastQ information
    #--------
    
    addTooltip(session, "select_FastQ", info_FastQ, placement = "left", trigger = 'hover', options = list(container = "body"))
    
    #--------
    # AB1 information
    #--------
    
    addTooltip(session, "create_guides", info_guides, placement = "left", trigger = 'hover', options = list(container = "body"))
  
    #--------
    # AB1 information
    #--------
    
    addTooltip(session, "edit_xls", info_xls, placement = "left", trigger = 'hover', options = list(container = "body"))
    
  }
  else
  {
    removeTooltip(session, "select_AB1")
    removeTooltip(session, "select_FastQ")
    removeTooltip(session, "create_guides")
    removeTooltip(session, "edit_xls")
  }
  })