################################################################################
# Metadata table screen, which controls how samples are named
################################################################################

save_table <- bsButton("update_table", 'update', icon =  icon("list-alt"),
                       style = "primary", block = TRUE)
guide_from_table <- bsButton("guide_from_table", 'create guide',
                       icon =  icon("area-chart"), style = "success", block = TRUE)


modal_table <- .bsModal(
  "modal_table", "Edit Table : ",
  "update_table", size = "large", close.button = FALSE,
  # Metadata table
  fluidRow(
    column(width = 12,
      helpText("Changes to the table will be automatically saved to the source file."),
      htmlOutput("table")
    )
  ),
  # Footer containing action buttons - adjust table or next step (create guide)
  footer = fluidRow(
    column(width = 6,
      fluidRow(column(width = 6),column(width = 6))),
    column(width = 6,
      fluidRow(column(width = 6, guide_from_table),
               column(width = 6, save_table)))
  )
)
