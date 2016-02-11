modal_help <- bsModal(
  "modal_help", "Create the reference guide : ",
  "setting_btn", size = "large",
  bsAlert("alertSave"),
  fluidRow(
    column(width = 12,
    includeMarkdown("./help_file/CrispRVariants_reference_manual.Rmd")     
  ),
  fluidRow(
    column(width = 6),
    column(width = 6)
  )))
  
