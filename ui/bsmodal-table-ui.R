################################################################################
# Define  options modal Table 
################################################################################

#upload_table <- fileInput('upload_table', 'Apply table', multiple = TRUE, width = "100%")

#upload_metadata <- fileInput('upload_metadata', 'Upload metadata', multiple = F, width = "100%")
#select_genome <- selectInput("select_genome", "Select the genome", choices = genlist, width = "100%")

save_table <- bsButton("update_table", 'update', icon =  icon("list-alt"), style = "primary", block = TRUE)
guide_from_table <- bsButton("guide_from_table", 'create guide', icon =  icon("area-chart"), style = "success", block = TRUE)
modal_table <- .bsModal(
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
    ), close.button = F)