################################################################################
# Metadata table screen, which controls how samples are named
################################################################################

save_table <- bsButton("update_table", 'update', icon =  icon("list-alt"), style = "primary", block = TRUE)
guide_from_table <- bsButton("guide_from_table", 'create guide', icon =  icon("area-chart"), style = "success", block = TRUE)

modal_table <- bsModal(
  "modal_table", "Edit Table : ",
  "update_table", size = "large",
  fluidRow(
    column(width = 12,
      helpText("Changes to the table will be automatically saved to the source file."),
      htmlOutput("table")
    )
  ),
  footer = fluidRow(
    column(width = 6,
      fluidRow(
        column(width = 6),
        column(width = 6)
      )),
    column(width = 6,
      fluidRow(
        column(width = 6, guide_from_table),
        column(width = 6, save_table)
        )
      )
    )
)
