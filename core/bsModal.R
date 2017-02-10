# This code is taken from the shinyBS Github repo for shinyBS_0.62
# Cran currently supports v0.61.  We are including this code to allow
# installation of CrispRVariantsLite requirements from CRAN/Bioconductor.
# https://github.com/ebailey78/shinyBS/blob/shinyBS3/R/bsModal.R
addAttribs <- function(tag, ...) {
  a <- list(...)
  for(i in seq(length(a))) {
    tag$attribs[names(a)[i]] = a[[i]]
  }
  return(tag)
}


.bsModal <- function(id, title, trigger, ..., size, footer = NULL,
                     close.button = TRUE, width = NULL) {
  if(!missing(size)) {
    if(size == "large") {
      size = "modal-lg"
    } else if(size == "small") {
      size = "modal-sm"
    }
    size <- paste("modal-dialog", size)
    width = NULL
  } else {
    size <- "modal-dialog"
  }
  
  if(is.null(footer)) {
    footer <- tagList()
  } 
  
  if(close.button) {
    footer <- shiny::tagAppendChild(footer, tagList(shiny::tags$button(type = "button",
                     class = "btn btn-default", "data-dismiss" = "modal", "Close")))
  }
  
  bsTag <- shiny::tags$div(class = size, 
               shiny::tags$div(class = "modal-content", 
                   shiny::tags$div(class = "modal-header",
                          shiny::tags$button(type = "button", 
                              class = "close", "data-dismiss" = "modal", 
                              shiny::tags$span(shiny::HTML("&times;"))),
                              shiny::tags$h4(class = "modal-title", title)
                           ),
                           shiny::tags$div(class = "modal-body", list(...)),
                           shiny::tags$div(class = "modal-footer", footer
                   )
               )
  )
  
  if(!is.null(width)) {
     bsTag <- addAttribs(bsTag, style = paste0("width: ", width, " !important;"))
  }
  
  bsTag <- shiny::tags$div(class = "modal sbs-modal fade", id = id, tabindex = "-1",
                          "data-sbs-trigger" = trigger, bsTag)
}