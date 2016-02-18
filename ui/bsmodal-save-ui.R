 ui_save <- uiOutput("ui_save")


modal_save <- bsModal(
  "modal_save", "DOWNLOAD DATA : ",
  "save_data", size = "small",
  bsAlert("alertAB1"),
  fluidRow(
    column(width = 12,
      ui_save
      )
  )
)
