#ui_save <- uiOutput("ui_save")
ui_options <- uiOutput("ui_options")

modal_save <- .bsModal("modal_save", "DOWNLOAD DATA : ",
  "save_data", size = "small", close.button = FALSE,
  bsAlert("alertAB1"),
  fluidRow(column(width = 12, ui_options))#,
  #footer = fluidRow(column(width = 3, ui_save))
)
