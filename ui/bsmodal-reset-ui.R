################################################################################
# Reset button on the plot panel to restart analysis
################################################################################

confirm <- bsButton("confirm", 'confirm',  style = "danger", block = TRUE)
cancel <- bsButton("cancel", 'cancel',  style = "primary", block = TRUE)
modal_reset <- bsModal(
  "modal_reset", "INFORMATION",
  "reset", size = "small",
  fluidRow(
    column(width = 12,
      h4("If you continue all the current data will be lost")
  )),
  footer = fluidRow(
    column(width = 6,
    confirm),
    column(width = 6,
    cancel)
  )
)
